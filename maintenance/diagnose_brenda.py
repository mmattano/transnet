#!/usr/bin/env python3
"""Find out why BRENDA is returning nothing.

The mouse build made 3,334 BRENDA queries without error and got an empty answer
every time. That can mean several different things, and they need different
fixes, so this walks through them in order and prints a verdict.

Run it in the shell where BRENDA_EMAIL and BRENDA_PASSWORD are exported::

    python maintenance/diagnose_brenda.py

It prints no secrets -- the output is safe to paste.
"""

import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Enzymes chosen because BRENDA definitely holds inhibitor data for them.
PROBES = [
    ("1.1.1.1", "alcohol dehydrogenase"),
    ("2.7.1.1", "hexokinase"),
    ("1.1.1.27", "L-lactate dehydrogenase"),
]


def rule(title):
    print()
    print("=" * 70)
    print(title)
    print("=" * 70)


def main():
    rule("1. Credentials")
    email = os.environ.get("BRENDA_EMAIL")
    password = os.environ.get("BRENDA_PASSWORD")
    # Only whether it is set: this output lands in log files.
    print(f"  BRENDA_EMAIL    : {'set' if email else 'NOT SET'}")
    print(f"  BRENDA_PASSWORD : {'set, ' + str(len(password)) + ' chars' if password else 'NOT SET'}")
    if not (email and password):
        print("\n  -> Export both, then rerun. Without them the client would")
        print("     prompt, which is why an unattended build can hang or fail.")
        return 1

    rule("2. Connecting")
    try:
        from transnet.api.brenda import BrendaClient
        client = BrendaClient()
        print("  client constructed")
    except Exception:
        traceback.print_exc()
        print("\n  -> Could not construct the client. Is `zeep` installed?")
        return 1

    rule("3. A call with no organism filter")
    unfiltered = {}
    for ec, name in PROBES:
        try:
            frame = client.get_inhibitors(ec, None)
            unfiltered[ec] = len(frame)
            print(f"  EC {ec:<10} ({name:<24}) {len(frame):>5} rows  "
                  f"columns={list(frame.columns)[:5]}")
            if len(frame):
                first = frame.iloc[0].to_dict()
                shown = {k: str(v)[:28] for k, v in list(first.items())[:4]}
                print(f"      first row: {shown}")
        except Exception as exc:
            unfiltered[ec] = -1
            print(f"  EC {ec:<10} RAISED: {type(exc).__name__}: {str(exc)[:120]}")

    rule("4. The same calls filtered to Mus musculus")
    filtered = {}
    for ec, name in PROBES:
        try:
            frame = client.get_inhibitors(ec, "Mus musculus")
            filtered[ec] = len(frame)
            print(f"  EC {ec:<10} {len(frame):>5} rows")
        except Exception as exc:
            filtered[ec] = -1
            print(f"  EC {ec:<10} RAISED: {type(exc).__name__}: {str(exc)[:120]}")

    rule("5. What organism names BRENDA actually has for EC 1.1.1.1")
    try:
        frame = client.get_inhibitors("1.1.1.1", None)
        if len(frame) and "organism" in frame.columns:
            organisms = frame["organism"].dropna().unique()
            print(f"  {len(organisms)} distinct organisms; first 10:")
            for organism in list(organisms)[:10]:
                print(f"    {organism!r}")
            mouse = [o for o in organisms if "mus" in str(o).lower()]
            print(f"  entries mentioning 'mus': {mouse[:5] if mouse else 'NONE'}")
        else:
            print("  no rows, or no 'organism' column to inspect")
    except Exception as exc:
        print(f"  RAISED: {type(exc).__name__}: {str(exc)[:120]}")

    rule("6. The enrichment path end to end")
    from transnet.api.brenda import brenda_cache_dir, brenda_enrich_proteins
    from transnet.biology.elements import Protein

    probe = Protein(uniprot_id="PROBE", ec_number=["1.1.1.1"])
    try:
        brenda_enrich_proteins(
            [probe], organism=None, fields=["activators", "inhibitors"]
        )
        print(f"  inhibitors found: {len(probe.inhibitors)}")
        print(f"  activators found: {len(probe.activators)}")
        if probe.inhibitors:
            print(f"    e.g. {probe.inhibitors[:5]}")
        enrichment_ok = bool(probe.inhibitors or probe.activators)
    except Exception as exc:
        print(f"  RAISED: {type(exc).__name__}: {str(exc)[:150]}")
        enrichment_ok = False

    rule("Verdict")
    any_unfiltered = any(n > 0 for n in unfiltered.values())
    any_filtered = any(n > 0 for n in filtered.values())

    if enrichment_ok and any_unfiltered:
        print("  READY. BRENDA returns data and the enrichment path works.")
        if not any_filtered:
            print()
            print("  Note: the 'Mus musculus' filter returned nothing while")
            print("  unfiltered queries did. Section 5 shows the spellings")
            print("  BRENDA uses. The build passes an organism name, so if this")
            print("  persists you will get no allosteric edges -- tell me and I")
            print("  will drop the filter.")
        print()
        print("  Next:")
        print("    python maintenance/build_networks.py --organisms mouse --brenda")
        return 0

    if not any_unfiltered and all(n == 0 for n in unfiltered.values()):
        print("  NOT READY. Every call succeeded but returned zero rows.")
        print()
        print("  If the columns are empty too, the request shape is wrong:")
        print("  BRENDA needs each field as a positional 'field*value' string")
        print("  and answers a malformed call with an empty result rather than")
        print("  an error. That bug was fixed; if you still see this, the")
        print("  remaining candidate is the account itself -- check it is")
        print("  activated by following the link in BRENDA's signup email.")
    elif any_unfiltered and not any_filtered:
        print("  PARTLY WORKING. Unfiltered calls return data; filtering by")
        print("  'Mus musculus' returns nothing -- the organism filter is the")
        print("  problem. Section 5 lists the spellings BRENDA actually holds.")
    else:
        print("  Calls raised rather than returning empty -- see the errors above.")

    print()
    print("  (no credentials appear in this output; safe to paste)")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
