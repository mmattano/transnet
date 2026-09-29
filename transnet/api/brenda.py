"""BRENDA Enzyme Database API client.

BRENDA (https://www.brenda-enzymes.org) is the world's largest enzyme
information system.  This module provides a Python wrapper around the
BRENDA SOAP/WSDL API, returning results as pandas DataFrames and
caching queries to disk to avoid hammering the server.

Authentication
--------------
BRENDA requires a registered account.  Credentials can be supplied:

1. Explicitly::

       client = BrendaClient(email="you@example.com", password="secret")

2. Via environment variables::

       BRENDA_EMAIL=you@example.com
       BRENDA_PASSWORD=secret   # plain text — hashed automatically

3. Via interactive prompt (fallback when neither is provided).

The password is always stored as its SHA-256 hex digest, never in plain text.

Usage
-----
>>> from transnet.api.brenda import BrendaClient
>>> client = BrendaClient(email="you@example.com", password="secret")
>>> km = client.get_km_values("1.1.1.1", organism="Homo sapiens")
>>> km.head()

Or use the standalone convenience functions which manage client creation
internally:

    >>> from transnet.api.brenda import brenda_get_km_values
    >>> df = brenda_get_km_values("1.1.1.1", organism="Homo sapiens",
    ...                           email="you@example.com", password="secret")
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import logging
import os
import re
import time
from typing import Dict, List, Optional

import pandas as pd

logger = logging.getLogger(__name__)

__all__ = [
    "BrendaClient",
    "brenda_get_km_values",
    "brenda_get_kcat_values",
    "brenda_get_ki_values",
    "brenda_get_inhibitors",
    "brenda_get_activators",
    "brenda_get_substrates",
    "brenda_get_products",
    "brenda_get_cofactors",
    "brenda_enrich_proteins",
]

_BRENDA_WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
_REQUEST_DELAY = 0.5  # seconds between consecutive SOAP calls


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _hash_password(password: str) -> str:
    """Return the SHA-256 hex digest of *password*."""
    return hashlib.sha256(password.encode("utf-8")).hexdigest()


def _is_hashed(password: str) -> bool:
    """Return True if *password* looks like a SHA-256 hex digest."""
    return len(password) == 64 and all(c in "0123456789abcdef" for c in password.lower())


def _build_params(ec_number: str, organism: Optional[str] = None, **extra: str) -> str:
    """Build the BRENDA parameter string for a query.

    Format: ``ecNumber*{ec}#organism*{org}#field1*value1#...``
    """
    parts = [f"ecNumber*{ec_number}"]
    if organism:
        parts.append(f"organism*{organism}")
    for key, value in extra.items():
        if value:
            parts.append(f"{key}*{value}")
    return "#".join(parts)


def _serialize(raw) -> list:
    """Convert a zeep result to a plain list of dicts."""
    try:
        from zeep.helpers import serialize_object
        return [dict(item) for item in serialize_object(raw)] if raw else []
    except Exception:
        return list(raw) if raw else []


def _raw_to_df(raw, rename: Optional[Dict[str, str]] = None) -> pd.DataFrame:
    """Convert a BRENDA SOAP response to a DataFrame."""
    records = _serialize(raw)
    if not records:
        return pd.DataFrame()
    df = pd.DataFrame(records)
    if rename:
        df = df.rename(columns=rename)
    return df


# ---------------------------------------------------------------------------
# BrendaClient class
# ---------------------------------------------------------------------------

# XML 1.0 forbids most C0 control characters. BRENDA's free-text commentary
# occasionally contains them (\x05 has been seen), and zeep then rejects the
# whole response as invalid XML -- every activator or inhibitor for that EC is
# lost, and retrying returns the same bytes.
_ILLEGAL_XML = re.compile(rb"[\x00-\x08\x0b\x0c\x0e-\x1f]")


def _sanitising_transport():
    """A zeep transport that strips XML-illegal bytes from responses."""
    from zeep.transports import Transport

    class SanitisingTransport(Transport):
        def post(self, address, message, headers):
            response = super().post(address, message, headers)
            content = response.content
            if content and _ILLEGAL_XML.search(content):
                cleaned = _ILLEGAL_XML.sub(b"", content)
                logger.debug(
                    f"BRENDA: removed {len(content) - len(cleaned)} "
                    f"XML-illegal byte(s) from a response"
                )
                response._content = cleaned
            return response

    return SanitisingTransport()


class BrendaClient:
    """Authenticated BRENDA SOAP client.

    Parameters
    ----------
    email : str, optional
        BRENDA account e-mail.  Falls back to ``BRENDA_EMAIL`` env var,
        then interactive prompt.
    password : str, optional
        BRENDA account password (plain text or SHA-256 hex).  Falls back to
        ``BRENDA_PASSWORD`` env var, then interactive prompt.

    Attributes
    ----------
    email : str
    password_hash : str
        SHA-256 hex digest of the password — never stored in plain text.
    """

    def __init__(
        self,
        email: Optional[str] = None,
        password: Optional[str] = None,
    ) -> None:
        self.email = email or os.environ.get("BRENDA_EMAIL") or self._prompt_email()
        raw_pw = password or os.environ.get("BRENDA_PASSWORD") or self._prompt_password()
        self.password_hash = raw_pw if _is_hashed(raw_pw) else _hash_password(raw_pw)
        self._client = None  # lazy initialisation
        self._parameter_cache: Dict[str, List[str]] = {}

    # ------------------------------------------------------------------
    # Auth prompts
    # ------------------------------------------------------------------

    @staticmethod
    def _prompt_email() -> str:
        return input("BRENDA account e-mail: ").strip()

    @staticmethod
    def _prompt_password() -> str:
        import getpass
        return getpass.getpass("BRENDA password: ")

    # ------------------------------------------------------------------
    # SOAP client initialisation
    # ------------------------------------------------------------------

    def _get_client(self):
        """Return a cached zeep Client, creating it on first call."""
        if self._client is None:
            try:
                from zeep import Client
            except ImportError as exc:
                raise ImportError(
                    "zeep is required for BRENDA API access. "
                    "Install it with: pip install zeep"
                ) from exc
            logger.info("Initialising BRENDA SOAP client…")
            self._client = Client(_BRENDA_WSDL, transport=_sanitising_transport())
        return self._client

    # ------------------------------------------------------------------
    # Low-level SOAP call
    # ------------------------------------------------------------------

    def _call(self, method: str, ec_number: str, organism: Optional[str] = None, **extra) -> list:
        """Make a single SOAP call and return a list of serialized result dicts.

        BRENDA's SOAP API does not take ordinary named parameters. Every field
        after the credentials must be passed as a positional ``"field*value"``
        string -- ``"ecNumber*1.1.1.1"``, ``"organism*Mus musculus"`` -- and
        every field the method declares must be present, empty ones as
        ``"field*"``. Passing plain values instead is not an error: BRENDA
        simply returns an empty result, which is why a misuse looks exactly
        like "this enzyme has no data".
        """
        client = self._get_client()
        service_method = getattr(client.service, method, None)
        if service_method is None:
            raise AttributeError(f"BRENDA SOAP service has no method '{method}'")

        parameter_names = self._parameter_names(method)

        values = {"ecNumber": ec_number}
        if organism:
            values["organism"] = organism
        values.update(extra)

        # email and password go through raw; everything else is "field*value".
        arguments = [self.email, self.password_hash]
        for name in parameter_names[2:]:
            arguments.append(f"{name}*{values.get(name, '')}")

        logger.debug(f"BRENDA {method}({arguments[2:]})")
        try:
            from zeep.exceptions import TransportError, Fault
            try:
                raw = service_method(*arguments)
                time.sleep(_REQUEST_DELAY)
                return _serialize(raw)
            except Fault as exc:
                raise RuntimeError(f"BRENDA SOAP fault for {method}: {exc}") from exc
            except TransportError as exc:
                # zeep puts the entire response body in the message, which
                # turned each failure into kilobytes of SOAP in the build log.
                message = str(exc).split("Content:")[0].strip()
                raise RuntimeError(
                    f"BRENDA transport error for {method}: {message[:300]}"
                ) from exc
        except ImportError:
            raise

    def _parameter_names(self, method: str) -> List[str]:
        """Parameter names for a SOAP method, in the order the WSDL declares.

        Read from the WSDL rather than hard-coded, because the fields differ
        per method (``getSubstrate`` has ``reactionPartners`` and no
        ``literature``, for instance).
        """
        cached = self._parameter_cache.get(method)
        if cached is not None:
            return cached

        client = self._get_client()
        names: List[str] = []
        try:
            service = list(client.wsdl.services.values())[0]
            port = list(service.ports.values())[0]
            operation = port.binding._operations.get(method)
            # zeep exposes the request fields as the elements of the body's
            # complex type; `.parts` does not exist on this element.
            names = [name for name, _ in operation.input.body.type.elements]
        except Exception as exc:
            logger.debug(f"Could not read WSDL parameters for {method}: {exc}")
            names = []

        if not names:
            # Every BRENDA getter starts the same way; this keeps a WSDL
            # parsing change from breaking the common filters.
            names = ["email", "password", "ecNumber", "organism"]
            logger.debug(f"Falling back to default parameter names for {method}")

        self._parameter_cache[method] = names
        return names

    # ------------------------------------------------------------------
    # Public query methods
    # ------------------------------------------------------------------

    def get_km_values(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Michaelis constant (Km) values for an enzyme.

        Parameters
        ----------
        ec_number : str
            EC number (e.g. ``"1.1.1.1"``).
        organism : str, optional
            Filter by organism.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, substrate, kmValue, kmValueMaximum,
            commentary, organism, ligandStructureId.
        """
        records = self._call("getKmValue", ec_number, organism)
        if not records:
            return pd.DataFrame()
        df = pd.DataFrame(records)
        return df

    def get_kcat_values(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Catalytic constant (kcat / turnover number) values.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, substrate, turnoverNumber,
            turnoverNumberMaximum, commentary, organism.
        """
        records = self._call("getTurnoverNumber", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_ki_values(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Inhibition constant (Ki) values.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, inhibitor, kiValue, kiValueMaximum,
            commentary, organism.
        """
        records = self._call("getKiValue", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_inhibitors(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Inhibiting compounds for an enzyme.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, inhibitor, commentary, organism.
        """
        records = self._call("getInhibitors", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_activators(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Activating compounds for an enzyme.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, activatingCompound, commentary, organism.
        """
        records = self._call("getActivatingCompound", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_substrates(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Natural substrates for an enzyme.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, naturalSubstrate, commentary, organism.
        """
        records = self._call("getNaturalSubstrate", ec_number, organism)
        if not records:
            # Fall back to any substrate
            records = self._call("getSubstrate", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_products(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Natural products for an enzyme.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, naturalProduct, commentary, organism.
        """
        records = self._call("getNaturalProduct", ec_number, organism)
        if not records:
            records = self._call("getProduct", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_cofactors(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Cofactor requirements for an enzyme.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, cofactor, commentary, organism.
        """
        records = self._call("getCofactor", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_metals(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> pd.DataFrame:
        """Metal ion requirements for an enzyme.

        Returns
        -------
        pd.DataFrame
            Columns: ecNumber, metals, commentary, organism.
        """
        records = self._call("getMetals", ec_number, organism)
        if not records:
            return pd.DataFrame()
        return pd.DataFrame(records)

    def get_all_kinetics(
        self,
        ec_number: str,
        organism: Optional[str] = None,
    ) -> Dict[str, pd.DataFrame]:
        """Fetch all available kinetic data for an EC number.

        Queries Km, Kcat, Ki, inhibitors, activators, substrates,
        products, and cofactors in one call.

        Returns
        -------
        dict mapping query name → DataFrame
        """
        results: Dict[str, pd.DataFrame] = {}
        queries = [
            ("km_values", self.get_km_values),
            ("kcat_values", self.get_kcat_values),
            ("ki_values", self.get_ki_values),
            ("inhibitors", self.get_inhibitors),
            ("activators", self.get_activators),
            ("substrates", self.get_substrates),
            ("products", self.get_products),
            ("cofactors", self.get_cofactors),
        ]
        for name, method in queries:
            try:
                results[name] = method(ec_number, organism)
            except Exception as exc:
                logger.warning(f"BRENDA {name} query failed for {ec_number}: {exc}")
                results[name] = pd.DataFrame()
        return results


# ---------------------------------------------------------------------------
# Cached standalone functions (module-level, no client object required)
# ---------------------------------------------------------------------------

def _get_client_from_kwargs(kwargs: dict) -> BrendaClient:
    """Pop email/password from kwargs and return a BrendaClient."""
    email = kwargs.pop("email", None)
    password = kwargs.pop("password", None)
    return BrendaClient(email=email, password=password)


def brenda_get_km_values(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return Km values for *ec_number* from BRENDA.

    Parameters
    ----------
    ec_number : str
        EC number.
    organism : str, optional
        Restrict results to this organism.
    email, password : str, optional
        BRENDA credentials.  Fall back to env vars ``BRENDA_EMAIL`` /
        ``BRENDA_PASSWORD`` if not supplied.

    Returns
    -------
    pd.DataFrame
    """
    client = BrendaClient(email=email, password=password)
    return client.get_km_values(ec_number, organism)


def brenda_get_kcat_values(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return kcat (turnover number) values for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_kcat_values(ec_number, organism)


def brenda_get_ki_values(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return Ki (inhibition constant) values for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_ki_values(ec_number, organism)


def brenda_get_inhibitors(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return inhibitor compounds for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_inhibitors(ec_number, organism)


def brenda_get_activators(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return activating compounds for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_activators(ec_number, organism)


def brenda_get_substrates(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return natural substrates for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_substrates(ec_number, organism)


def brenda_get_products(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return natural products for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_products(ec_number, organism)


def brenda_get_cofactors(
    ec_number: str,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
) -> pd.DataFrame:
    """Return cofactor requirements for *ec_number* from BRENDA."""
    client = BrendaClient(email=email, password=password)
    return client.get_cofactors(ec_number, organism)


# ---------------------------------------------------------------------------
# Integration helper — enrich a list of Protein objects
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# On-disk cache
#
# BRENDA is queried one EC number at a time over SOAP, so enriching a proteome
# takes hours. Caching each answer makes the job resumable: interrupt it, run
# the same command again, and only the EC numbers not yet fetched are queried.
# ---------------------------------------------------------------------------

def brenda_cache_dir() -> Path:
    """Directory holding cached BRENDA answers.

    Override with the ``TRANSNET_BRENDA_CACHE`` environment variable.
    """
    configured = os.environ.get("TRANSNET_BRENDA_CACHE")
    if configured:
        return Path(configured)
    return Path.home() / ".cache" / "transnet" / "brenda"


def _cache_path(ec: str, field: str, organism: Optional[str]) -> Path:
    """One file per (EC, field, organism); the EC is slugged for the filename."""
    scope = (organism or "all").replace(" ", "_")
    slug = ec.replace(".", "_").replace("/", "_")
    return brenda_cache_dir() / scope / field / f"{slug}.json"


def _cache_read(ec: str, field: str, organism: Optional[str]):
    path = _cache_path(ec, field, organism)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        # A half-written file from an interrupted run: drop it and re-fetch.
        try:
            path.unlink()
        except OSError:
            pass
        return None


def _cache_write(ec: str, field: str, organism: Optional[str], names: list) -> None:
    path = _cache_path(ec, field, organism)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Write-then-rename so an interrupt cannot leave a truncated file behind.
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(names))
    temporary.replace(path)


def brenda_cache_stats() -> Dict[str, int]:
    """How many answers are cached, by field. Useful for reporting progress."""
    root = brenda_cache_dir()
    if not root.exists():
        return {}
    counts: Dict[str, int] = {}
    for path in root.glob("*/*/*.json"):
        counts[path.parent.name] = counts.get(path.parent.name, 0) + 1
    return dict(sorted(counts.items()))


def brenda_enrich_proteins(
    proteins: list,
    organism: Optional[str] = None,
    email: Optional[str] = None,
    password: Optional[str] = None,
    fields: Optional[List[str]] = None,
) -> None:
    """Enrich a list of Protein objects with BRENDA kinetic data in-place.

    For each unique EC number found across *proteins*, queries BRENDA for
    inhibitors, activators, substrates, and products, then assigns the
    compound names to the matching Protein objects.

    Parameters
    ----------
    proteins : list of Protein
        The proteins to enrich.  Each must have an ``ec_number`` list
        attribute and ``inhibitors``, ``activators``, ``substrates``,
        ``products`` list attributes (as defined in
        :class:`transnet.biology.elements.Protein`).
    organism : str, optional
        Restrict BRENDA queries to this organism name (e.g. "Homo sapiens").
        Leaving this as None returns data from all organisms — useful when
        the proteome contains proteins from multiple species.
    email, password : str, optional
        BRENDA credentials.
    fields : list of str, optional
        Subset of ``['inhibitors', 'activators', 'substrates', 'products',
        'cofactors']`` to fetch.  Fetches all by default.

    Notes
    -----
    Proteins without EC numbers are silently skipped.
    """
    if fields is None:
        fields = ["inhibitors", "activators", "substrates", "products", "cofactors"]

    client = BrendaClient(email=email, password=password)

    # Collect unique EC numbers and map each to the proteins that have it
    ec_to_proteins: Dict[str, list] = {}
    for protein in proteins:
        for ec in (protein.ec_number or []):
            ec = ec.strip()
            if not ec:
                continue
            ec_to_proteins.setdefault(ec, []).append(protein)

    n_ecs = len(ec_to_proteins)
    logger.info(f"Enriching {len(proteins)} proteins across {n_ecs} unique EC numbers from BRENDA")

    #: field -> (client method, candidate column names, protein attribute)
    field_spec = {
        "inhibitors": ("get_inhibitors",
                       ["inhibitor", "inhibitors"], "inhibitors"),
        "activators": ("get_activators",
                       ["activatingCompound", "activator", "activators"],
                       "activators"),
        "substrates": ("get_substrates",
                       ["naturalSubstrate", "substrate", "substrates"],
                       "substrates"),
        "products": ("get_products",
                     ["naturalProduct", "product", "products"], "products"),
        "cofactors": ("get_cofactors",
                      ["cofactor", "cofactors"], "cofactors"),
    }

    n_cached = n_fetched = 0

    for idx, (ec, ec_proteins) in enumerate(ec_to_proteins.items(), 1):
        for field in fields:
            spec = field_spec.get(field)
            if spec is None:
                continue
            method_name, columns, attribute = spec

            names = _cache_read(ec, field, organism)
            if names is None:
                try:
                    frame = getattr(client, method_name)(ec, organism)
                    column = _pick_col(frame, columns)
                    names = (
                        frame[column].dropna().unique().tolist() if column else []
                    )
                except Exception as exc:
                    logger.warning(f"BRENDA {field} failed for EC {ec}: {exc}")
                    # Not cached: a failure should be retried on the next run,
                    # not remembered as "this EC has none".
                    continue
                _cache_write(ec, field, organism, names)
                n_fetched += 1
            else:
                n_cached += 1

            if not names:
                continue
            for protein in ec_proteins:
                setattr(protein, attribute,
                        _merge_unique(getattr(protein, attribute, []) or [], names))

        if idx % 50 == 0 or idx == n_ecs:
            logger.info(
                f"BRENDA {idx}/{n_ecs} EC numbers "
                f"({n_cached} from cache, {n_fetched} fetched)"
            )

    logger.info(
        f"BRENDA enrichment complete: {n_cached} cached answers reused, "
        f"{n_fetched} fetched. Cache: {brenda_cache_dir()}"
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the first column name from *candidates* that exists in *df*."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _merge_unique(existing: list, new_items: list) -> list:
    """Return *existing* extended with items from *new_items* not already present."""
    existing_set = set(existing)
    for item in new_items:
        if item not in existing_set:
            existing.append(item)
            existing_set.add(item)
    return existing
