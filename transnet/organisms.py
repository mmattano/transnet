"""Which organisms a reference network can be built for, and what defines one.

A TransNet reference network is assembled from public databases, and every one
of those databases names organisms differently: KEGG uses a three-letter code,
UniProt a taxon id that is sometimes a strain rather than a species, Ensembl a
lowercase species name, ChIP-Atlas a genome assembly. An organism is therefore
defined by one entry in :data:`ORGANISMS` that says what to ask each database
for.

Five organisms are configured here and built into ``data/<organism>/latest``.
Adding a sixth means adding an entry, either by editing this module or, from
your own code, with :func:`register_organism`:

    >>> from transnet.organisms import register_organism
    >>> register_organism(
    ...     "zebrafish",
    ...     kegg_org="dre",
    ...     organism_full="Danio rerio",
    ...     ncbi_org="7955",
    ...     ensembl_org="danio_rerio",
    ...     ensembl_release=109,
    ...     genome_chip="danRer11",
    ... )

The keys
--------

``kegg_org`` : str
    KEGG's three-letter organism code, which decides the reactions, compounds,
    EC numbers and pathways. ``transnet.api.kegg.kegg_list_organisms()`` lists
    them.
``organism_full`` : str
    The species name as BRENDA writes it, used to restrict BRENDA's allosteric
    effectors to this organism.
``ncbi_org`` : str
    NCBI taxonomy id, used for STRING and as the default UniProt taxon.
``uniprot_org`` : str, optional
    UniProt taxon, for the case where UniProt files its reviewed entries under
    a strain rather than the species. Only needed when it differs from
    ``ncbi_org``, and getting it wrong is the most common way to build a
    network that looks complete and is not: see
    ``MIN_PLAUSIBLE_PROTEOME`` in ``maintenance/build_networks.py``.
``ensembl_org``, ``ensembl_release`` : str, int
    Ensembl species name and release, for the transcriptome.
``genome_chip`` : str
    ChIP-Atlas genome assembly, which decides the transcription-factor to
    target-gene edges. ChIP-Atlas carries a limited set of assemblies, and the
    one it carries is not always the current one.
``transcriptome_source`` : {"ensembl", "kegg"}, optional
    Where the transcriptome comes from. Defaults to Ensembl; set to ``"kegg"``
    for organisms Ensembl does not carry.
``expected_absent`` : list of str, optional
    Edge types that cannot exist for this organism, because no database
    supplies them. The build report names a missing edge type as a gap unless
    it is listed here, so this is how a known absence is distinguished from a
    failed download.
"""

from typing import Any, Dict, List

__all__ = ["ORGANISMS", "ORGANISM_CONFIG", "register_organism", "organism_config"]


#: Organism name -> the identifiers each database needs. See the module
#: docstring for what each key means.
ORGANISMS: Dict[str, Dict[str, Any]] = {
    "human": {
        "kegg_org": "hsa",
        "organism_full": "Homo sapiens",
        "ncbi_org": "9606",
        "ensembl_org": "homo_sapiens",
        "ensembl_release": 109,
        "genome_chip": "hg38",
    },
    "mouse": {
        "kegg_org": "mmu",
        "organism_full": "Mus musculus",
        "ncbi_org": "10090",
        "ensembl_org": "mus_musculus",
        "ensembl_release": 109,
        "genome_chip": "mm10",
    },
    "rat": {
        "kegg_org": "rno",
        "organism_full": "Rattus norvegicus",
        "ncbi_org": "10116",
        "ensembl_org": "rattus_norvegicus",
        "ensembl_release": 109,
        # rn6 is the only rat assembly ChIP-Atlas carries (52 factors), which is
        # why rat transcription-factor coverage is thin in every rat study.
        "genome_chip": "rn6",
    },
    "yeast": {
        "kegg_org": "sce",
        "organism_full": "Saccharomyces cerevisiae",
        "ncbi_org": "4932",
        # UniProt files reviewed yeast entries under strain S288c: 6,733
        # proteins, against 43 under the species id. STRING keeps 4932.
        "uniprot_org": "559292",
        "ensembl_org": "saccharomyces_cerevisiae",
        "ensembl_release": 109,
        "genome_chip": "sacCer3",
    },
    "ecoli": {
        "kegg_org": "eco",
        "organism_full": "Escherichia coli",
        "ncbi_org": "511145",
        # UniProt files reviewed E. coli entries under K-12: 4,531 proteins,
        # against 0 under MG1655 (511145, which holds only 8 at all).
        "uniprot_org": "83333",
        # Ensembl (and so pyensembl) does not carry bacteria; build the
        # transcriptome from KEGG genes instead.
        "transcriptome_source": "kegg",
        # ChIP-Atlas carries no bacterial genome, so this edge type cannot
        # exist for E. coli. Its absence is not a build failure.
        "expected_absent": ["transcriptional_regulation"],
        "ensembl_org": "escherichia_coli_str_k_12_substr_mg1655",
        "ensembl_release": 109,
        "genome_chip": "eschColi_K12",
    },
}

#: Kept so existing imports of the old name keep working.
ORGANISM_CONFIG = ORGANISMS

#: The keys every entry must carry. ``uniprot_org``, ``transcriptome_source``
#: and ``expected_absent`` are optional and defaulted by the builder.
REQUIRED_KEYS = (
    "kegg_org", "organism_full", "ncbi_org",
    "ensembl_org", "ensembl_release", "genome_chip",
)


def register_organism(name: str, replace: bool = False, **config: Any) -> Dict[str, Any]:
    """Add an organism to the registry.

    Parameters
    ----------
    name : str
        What to call it, as passed to ``build_networks.py --organisms``.
    replace : bool
        Overwrite an existing entry. Off by default, so a typo that collides
        with a built organism is an error rather than a silent redefinition of
        the network a study depends on.
    **config
        The keys above. The six in :data:`REQUIRED_KEYS` must all be given.

    Returns
    -------
    dict
        The stored entry.

    Raises
    ------
    ValueError
        If the name is already registered and ``replace`` is False, or if a
        required key is missing.
    """
    if name in ORGANISMS and not replace:
        raise ValueError(
            f"{name!r} is already registered; pass replace=True to overwrite it"
        )
    missing = [key for key in REQUIRED_KEYS if key not in config]
    if missing:
        raise ValueError(
            f"{name!r} is missing required key(s): {', '.join(missing)}. "
            f"See transnet.organisms for what each one is."
        )
    ORGANISMS[name] = dict(config)
    return ORGANISMS[name]


def organism_config(name: str) -> Dict[str, Any]:
    """The entry for one organism, with a message naming the alternatives."""
    try:
        return ORGANISMS[name]
    except KeyError:
        raise KeyError(
            f"Unknown organism {name!r}. Configured: "
            f"{', '.join(sorted(ORGANISMS))}. Add one with "
            f"transnet.organisms.register_organism()."
        ) from None


def available_organisms() -> List[str]:
    """Names in the registry, sorted."""
    return sorted(ORGANISMS)
