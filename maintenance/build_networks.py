"""Build the pre-built organism networks shipped in ``data/``.

Run from anywhere::

    python maintenance/build_networks.py --organisms mouse
    python maintenance/build_networks.py --organisms human mouse --brenda

``--brenda`` adds the allosteric edges, which is what makes the metabolite
regulation axis available; it needs BRENDA_EMAIL / BRENDA_PASSWORD.
"""

import os
import sys
import argparse
import logging
import traceback
import pickle
import shutil
from datetime import datetime

# Run from a clone without `pip install -e .` first. Python puts the *script's*
# directory on sys.path, not the working directory, so the repository root has
# to be added explicitly -- without this the import below fails wherever the
# package is not installed.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO_ROOT)

import pandas as pd

from transnet.biology.transnet import Transnet
from transnet.biology.layers import Pathways, Reactions, Proteome, Metabolome, Transcriptome

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("network_build.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Organism configuration lives in the package, so a network can be defined
# without editing this script: transnet.organisms.register_organism().
from transnet.organisms import ORGANISMS as ORGANISM_CONFIG  # noqa: E402


# ---------------------------------------------------------------------------
# Checkpointing
#
# A full build is hours of API calls, so it has to survive being interrupted.
# Each layer is written to disk as soon as it is finished; rerunning the same
# command reloads what is already there and continues from the first unfinished
# step. Nothing is ever partially written -- a checkpoint appears only after its
# step returns -- so Ctrl-C costs at most the step in flight.
# ---------------------------------------------------------------------------

def checkpoint_dir(output_dir, organism):
    return os.path.join(output_dir, ".checkpoints", organism)


def load_checkpoint(output_dir, organism, name):
    """Return a finished layer, or None if this step has not run yet."""
    path = os.path.join(checkpoint_dir(output_dir, organism), f"{name}.pkl")
    if not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as handle:
            layer = pickle.load(handle)
        logger.info(f"  [resume] {organism}/{name} loaded from checkpoint")
        return layer
    except Exception as exc:
        logger.warning(
            f"  [resume] {organism}/{name} checkpoint unreadable ({exc}); "
            f"rebuilding this step"
        )
        return None


def save_checkpoint(output_dir, organism, name, layer):
    """Persist a finished layer, atomically."""
    directory = checkpoint_dir(output_dir, organism)
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, f"{name}.pkl")
    temporary = path + ".tmp"
    with open(temporary, "wb") as handle:
        pickle.dump(layer, handle)
    os.replace(temporary, path)
    logger.info(f"  [checkpoint] {organism}/{name} saved")


def step(output_dir, organism, name, build, resume=True, chain=None):
    """Run one build step, or reload it if it already completed.

    ``chain`` links steps that each build on the previous one's result (the
    proteome: base -> STRING -> metabolites -> BRENDA -> ChIP-Atlas). Once a
    step in a chain actually runs, every later checkpoint in that chain was
    made from an older object and is stale: loading it would silently throw
    away what was just computed. Adding ``--brenda`` to an organism whose
    ChIP-Atlas step was already checkpointed did exactly that. So after the
    first rebuilt step, later steps in the chain rebuild too.
    """
    stale = chain is not None and chain.get("rebuilt")
    if resume and not stale:
        existing = load_checkpoint(output_dir, organism, name)
        if existing is not None:
            return existing
    elif resume and stale and os.path.exists(
            os.path.join(checkpoint_dir(output_dir, organism), f"{name}.pkl")):
        logger.info(f"  [resume] {organism}/{name} checkpoint is stale -- an earlier "
                    f"step it builds on was redone -- rebuilding")
    layer = build()
    save_checkpoint(output_dir, organism, name, layer)
    if chain is not None:
        chain["rebuilt"] = True
    return layer


def clear_checkpoints(output_dir, organism):
    directory = checkpoint_dir(output_dir, organism)
    if os.path.isdir(directory):
        shutil.rmtree(directory)
        logger.info(f"Cleared checkpoints for {organism}")


def generate_reactions_file(output_dir):
    """
    Generate a master reactions file to be used by all organisms.
    
    Parameters:
    -----------
    output_dir : str
        Directory to save the reactions file
    
    Returns:
    --------
    str
        Path to the reactions file
    """
    logger.info("Generating master reactions file")
    
    # Create reactions layer
    reactions = Reactions()
    
    # Populate from API
    reactions.populate(from_api=True)
    
    logger.info(f"Retrieved {len(reactions.reactions)} reactions from KEGG API")
    
    # Create directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Save reactions to file
    reactions_file = os.path.join(output_dir, "master_reactions.csv")
    
    def _to_semisep(val):
        """Safely convert a list/ndarray/str/None to a semicolon-separated string."""
        import numpy as np
        if val is None:
            return ""
        if isinstance(val, np.ndarray):
            items = [str(v) for v in val.flat if v is not None and str(v) != 'nan']
        elif isinstance(val, list):
            items = [str(v) for v in val if v is not None]
        else:
            items = [str(val)] if str(val) not in ('', 'nan') else []
        return ";".join(items)

    # Convert reactions to DataFrame
    reactions_data = []
    for reaction in reactions.reactions:
        stoich_subs = []
        if reaction.stoichiometry_substrates is not None:
            stoich_subs = [str(x) if x is not None else "None" for x in reaction.stoichiometry_substrates]
        stoich_prods = []
        if reaction.stoichiometry_products is not None:
            stoich_prods = [str(x) if x is not None else "None" for x in reaction.stoichiometry_products]

        reactions_data.append({
            'reaction': reaction.id,
            'name': reaction.name,
            'equation': reaction.equation,
            'definition': reaction.definition,
            'enzyme': _to_semisep(reaction.enzyme),
            'substrates': _to_semisep(reaction.substrates),
            'products': _to_semisep(reaction.products),
            'stoichiometry_substrates': ";".join(stoich_subs),
            'stoichiometry_products': ";".join(stoich_prods),
            'reversible': getattr(reaction, 'reversible', True),
        })
    
    reactions_df = pd.DataFrame(reactions_data)
    reactions_df.to_csv(reactions_file, index=False)
    
    logger.info(f"Saved {len(reactions_data)} reactions to {reactions_file}")
    
    return reactions_file

def load_reactions_from_file(file_path):
    """
    Load reactions from a CSV file.
    
    Parameters:
    -----------
    file_path : str
        Path to the reactions CSV file
    
    Returns:
    --------
    pd.DataFrame
        DataFrame containing reaction information
    """
    logger.info(f"Loading reactions from {file_path}")

    try:
        from transnet.datasets import read_master_reactions

        reactions_df = read_master_reactions(file_path)
        logger.info(f"Loaded {len(reactions_df)} reactions from file")

        # Reversibility is new; older reaction tables predate it and KEGG's
        # default arrow is reversible.
        if 'reversible' not in reactions_df.columns:
            reactions_df['reversible'] = True
        else:
            reactions_df['reversible'] = reactions_df['reversible'].fillna(True).astype(bool)

        return reactions_df
    except Exception as e:
        logger.error(f"Error loading reactions from {file_path}: {e}")
        logger.error(traceback.format_exc())
        return pd.DataFrame()


#: A reviewed proteome smaller than this means the UniProt query went wrong
#: (usually the taxon), not that the organism is small: E. coli, the smallest
#: organism built here, has 4,531 reviewed proteins.
MIN_PLAUSIBLE_PROTEOME = 1000


def _report_build_quality(transnet, organism, brenda, config=None):
    """Say what actually got built, and name anything missing.

    A build that finishes is not necessarily a build that worked: an
    unreachable database leaves a layer empty, and the network is then quietly
    missing a whole class of edge. Rather than logging "Success" and moving on,
    list what is there and warn about what is not.
    """
    from transnet.biology.schema import available_edge_types, available_layers

    graph = transnet.generate_graph()
    edge_types = available_edge_types(graph)

    logger.info(f"--- {organism} build report ---")
    logger.info(f"  {graph.number_of_nodes():,} nodes, {graph.number_of_edges():,} edges")
    logger.info(f"  layers: {available_layers(graph)}")
    for edge_type, count in edge_types.items():
        logger.info(f"    {edge_type:<30}{count:>9,}")

    problems = []
    counts = {
        "Transcriptome": len(getattr(transnet.transcriptome, "genes", []) or []),
        "Proteome": len(getattr(transnet.proteome, "proteins", []) or []),
        "Metabolome": len(getattr(transnet.metabolome, "metabolites", []) or []),
        "Reactions": len(getattr(transnet.reactions, "reactions", []) or []),
    }
    for layer, count in counts.items():
        if count == 0:
            problems.append(f"the {layer} layer is empty")

    # Presence is not enough. Yeast once built with 43 proteins -- the
    # species taxon instead of the strain UniProt files reviewed entries
    # under -- and every edge type was technically present, just tiny.
    n_proteins = counts["Proteome"]
    if 0 < n_proteins < MIN_PLAUSIBLE_PROTEOME:
        problems.append(
            f"the Proteome has only {n_proteins} proteins, implausibly few "
            f"for a reference network -- check the UniProt taxon "
            f"(config 'uniprot_org')"
        )

    expected = {
        "translation": "gene-protein links (needs both layers populated)",
        "catalysis": "enzyme-reaction links (needs EC numbers on proteins)",
        "substrate": "reaction inputs",
        "product": "reaction outputs",
        "protein_interaction": "STRING interactions",
        "transcriptional_regulation": "ChIP-Atlas TF targets",
    }
    if brenda:
        expected["allosteric_inhibition"] = "BRENDA allosteric regulation"

    expected_absent = set((config or {}).get("expected_absent", []))
    for edge_type, description in expected.items():
        if edge_type in edge_types:
            continue
        if edge_type in expected_absent:
            logger.info(
                f"  no {edge_type} edges -- expected for {organism}: "
                f"{description} has no data for this organism"
            )
            continue
        problems.append(f"no {edge_type} edges -- {description}")

    if problems:
        logger.warning(f"  {organism}: built with gaps")
        for problem in problems:
            logger.warning(f"    - {problem}")
        logger.warning(
            "  The network is usable, but analyses depending on the missing "
            "relationships will report weaker evidence."
        )
    else:
        logger.info(f"  {organism}: all expected layers and edge types present")
    return problems


def build_organism_network(organism, output_dir, reactions_file=None, debug=False,
                           brenda=False, resume=True, reviewed_only=True,
                           chip_score_threshold=100.0, signaling=False):
    """
    Build a network for a specific organism.
    
    Parameters:
    -----------
    organism : str
        Organism identifier (human, mouse, yeast, ecoli)
    output_dir : str
        Directory to save the network
    reactions_file : str, optional
        Path to the master reactions file
    debug : bool
        If True, run in debug mode with smaller network components
    resume : bool
        If True (default) reload any step that already completed, so an
        interrupted build continues instead of starting over. Pass False, or
        use --no-resume, to rebuild from scratch.
    signaling : bool
        If True, add a Signaling layer from KEGG signal-transduction KGML
        (kinase-substrate and kinase-TF relations). Needed by studies with
        phosphoproteomics; every analysis works without it.
    chip_score_threshold : float
        Minimum mean ChIP-Atlas binding score for a transcriptional
        regulation edge. See
        :meth:`Proteome.get_transcription_factor_targets`.
    reviewed_only : bool
        If True (default) build the proteome from reviewed Swiss-Prot entries
        only. Mouse drops from ~88,000 proteins to ~17,000; the remainder are
        unreviewed TrEMBL isoforms that make STRING time out and contribute
        little. Pass --all-proteins for the full set.
    brenda : bool
        If True, fetch allosteric activators and inhibitors from BRENDA.
        Without them the network has no allosteric edges and the metabolite
        regulation axis cannot be assessed. Credentials come from the
        BRENDA_EMAIL / BRENDA_PASSWORD environment variables.
    
    Returns:
    --------
    bool
        Whether the network build was successful
    """
    if organism not in ORGANISM_CONFIG:
        logger.error(f"Unknown organism: {organism}")
        return False
    
    config = ORGANISM_CONFIG[organism]
    logger.info(f"Building network for {organism}")
    
    try:
        # Create network layers
        # 1. Pathways layer
        def _build_pathways():
            logger.info("Building pathways layer")
            layer = Pathways()
            layer.kegg_organism = config["kegg_org"]
            layer.populate()
            layer.fill_pathways()
            logger.info(f"Created pathways layer with {len(layer.pathways)} pathways")
            return layer

        pathways = step(output_dir, organism, "pathways", _build_pathways, resume)
        
        # 2. Reactions layer
        def _build_reactions():
            logger.info("Building reactions layer")
            layer = Reactions()
            if reactions_file and os.path.exists(reactions_file):
                reactions_df = load_reactions_from_file(reactions_file)
                layer.populate(from_api=False, df=reactions_df)
                logger.info(f"Created reactions layer with {len(layer.reactions)} reactions from file")
            else:
                logger.warning("Reactions file not found, using API instead")
                layer.populate(from_api=True)
                logger.info(f"Created reactions layer with {len(layer.reactions)} reactions from API")
            return layer

        reactions = step(output_dir, organism, "reactions", _build_reactions, resume)
        
        # 3. Proteome layer -- the expensive one, so each external lookup is
        # checkpointed separately. STRING and BRENDA are hours apiece; losing
        # one should not cost the other.
        def _build_proteome():
            logger.info("Building proteome layer")
            layer = Proteome()
            layer.ncbi_organism = config["ncbi_org"]
            layer.kegg_organism = config["kegg_org"]
            layer.populate(
                uniprot=True, reviewed_only=reviewed_only,
                uniprot_organism=config.get("uniprot_org"),
            )
            if debug:
                layer.proteins = layer.proteins[:100]
            logger.info(f"Populated proteome with {len(layer.proteins)} proteins")
            return layer

        proteome_chain = {"rebuilt": False}
        proteome = step(output_dir, organism, "proteome_base",
                        _build_proteome, resume, chain=proteome_chain)

        def _add_string():
            logger.info("Getting protein-protein interactions (STRING)")
            ok = proteome.get_interaction_partners()
            if ok is False:
                # Raise so no checkpoint is written: a timeout must be retried
                # on the next run, not remembered as a completed step.
                raise RuntimeError(
                    "STRING lookup failed (see the error above). Rerun to "
                    "retry it; earlier steps are already checkpointed."
                )
            return proteome

        proteome = step(output_dir, organism, "proteome_string",
                        _add_string, resume, chain=proteome_chain)

        def _add_metabolites():
            logger.info("Getting protein-metabolite associations (KEGG)")
            proteome.get_metabolites()
            return proteome

        proteome = step(output_dir, organism, "proteome_metabolites",
                        _add_metabolites, resume, chain=proteome_chain)

        if brenda:
            def _add_brenda():
                logger.info("Getting allosteric regulators from BRENDA")
                logger.info(
                    "  answers are cached per EC number, so an interrupted "
                    "run resumes where it stopped"
                )
                try:
                    proteome.get_brenda_kinetics(
                        organism=config.get("organism_full"),
                        fields=["activators", "inhibitors"],
                    )
                except KeyboardInterrupt:
                    # Let the cache keep what it has and stop cleanly.
                    logger.warning("BRENDA interrupted; cached answers are kept")
                    raise
                except Exception as exc:
                    logger.error(f"BRENDA enrichment failed: {exc}")
                    logger.warning(
                        "Continuing without allosteric edges; the metabolite "
                        "regulation axis will be unavailable for this network"
                    )
                n_effectors = sum(
                    1 for protein in proteome.proteins
                    if protein.activators or protein.inhibitors
                )
                logger.info(
                    f"BRENDA: {n_effectors} enzymes with activators/inhibitors"
                )
                return proteome

            proteome = step(output_dir, organism, "proteome_brenda",
                            _add_brenda, resume, chain=proteome_chain)
        else:
            logger.info(
                "Skipping BRENDA (pass --brenda to include allosteric edges)"
            )

        logger.info(f"Completed proteome layer for {organism}")
        
        # 4. Metabolome layer
        def _build_metabolome():
            logger.info("Building metabolome layer")
            layer = Metabolome()
            layer.populate()
            logger.info(f"Created metabolome layer with {len(layer.metabolites)} metabolites")
            return layer

        metabolome = step(output_dir, organism, "metabolome",
                          _build_metabolome, resume)
        
        # 5. Transcriptome layer
        def _build_transcriptome():
            logger.info("Building transcriptome layer")
            layer = Transcriptome()
            layer.kegg_organism = config["kegg_org"]
            layer.organism_full = config["ensembl_org"]

            if debug:
                logger.info("Debug mode - skipping transcriptome population")
                layer.genes = []
                return layer

            try:
                if config.get("transcriptome_source") == "kegg":
                    # populate(kegg_api=True) resolves NCBI ids itself, so no
                    # fill_gene_info() pass is needed afterwards.
                    layer.populate(
                        kegg_organism=config["kegg_org"], kegg_api=True
                    )
                else:
                    layer.populate(
                        ensembl=True,
                        ensembl_release=config["ensembl_release"],
                    )
                    layer.fill_gene_info()
                logger.info(f"Created transcriptome layer with {len(layer.genes)} genes")
            except Exception as e:
                logger.error(f"Error populating transcriptome: {e}")
                logger.error(traceback.format_exc())
                logger.warning("Will continue with empty transcriptome layer")
                layer.genes = []
            return layer

        transcriptome = step(output_dir, organism, "transcriptome",
                             _build_transcriptome, resume)

        # 6. Signaling layer, from KEGG signal-transduction KGML. Optional: a
        # study without phosphoproteomics has nothing to map onto it, and every
        # analysis reads the layers that are present.
        signaling_layer = None
        if signaling and not debug:
            def _build_signaling():
                from transnet.api import kegg_signaling_relations
                from transnet.biology.layers import Signaling

                logger.info("Building signaling layer (KEGG KGML)")
                relations = kegg_signaling_relations(config["kegg_org"])
                if relations.empty:
                    logger.warning(
                        "KEGG returned no signaling relations; the Signaling "
                        "layer will be absent from this network"
                    )
                    return Signaling()

                layer = Signaling()
                layer.populate_from_table(
                    relations, source_column="source", target_column="target",
                    sign_column="sign", target_type_column="target_type",
                    evidence="KEGG signaling KGML",
                )
                subtypes = relations["subtype"].value_counts().to_dict()
                logger.info(
                    f"Signaling layer: {len(layer.nodes)} nodes from "
                    f"{len(relations)} relations ({subtypes})"
                )
                return layer

            signaling_layer = step(output_dir, organism, "signaling",
                                   _build_signaling, resume)
        elif signaling:
            logger.info("Debug mode - skipping signaling layer")

        # Create the integrated network
        logger.info("Creating integrated network")
        transnet = Transnet(
            name=f"{organism}_network",
            pathways=pathways,
            reactions=reactions,
            proteome=proteome,
            metabolome=metabolome,
            transcriptome=transcriptome
        )
        if signaling_layer is not None:
            transnet.signaling = signaling_layer
        
        # Transcription factor targets -- checkpointed because ChIP-Atlas is
        # another long download.
        if not debug:
            def _add_tf_targets():
                logger.info("Getting transcription factor targets (ChIP-Atlas)")
                ok = proteome.get_transcription_factor_targets(
                    genome_ChIP=config["genome_chip"],
                    distance_ChIP=5,
                    score_threshold=chip_score_threshold,
                )
                if ok is False:
                    # As for STRING: raise so no checkpoint is written and the
                    # next run retries instead of remembering a failed step.
                    raise RuntimeError(
                        "ChIP-Atlas lookup failed (see the error above). "
                        "Rerun to retry it; earlier steps are already "
                        "checkpointed."
                    )
                # Individual factors can still be lost to timeouts even though
                # the step as a whole succeeded. A few is normal; a lot means
                # the download was degraded, and checkpointing that would bake
                # the gap into the network permanently.
                lost = getattr(proteome, "chip_atlas_failed_factors", [])
                n_ok = sum(1 for protein in proteome.proteins
                           if protein.transcription_factor_targets)
                if lost and len(lost) > max(5, 0.1 * (n_ok + len(lost))):
                    raise RuntimeError(
                        f"ChIP-Atlas dropped {len(lost)} factor(s) to network "
                        f"errors ({lost[:5]}...), more than 10% of the "
                        f"{n_ok + len(lost)} attempted. Not checkpointing a "
                        f"degraded download -- rerun to retry."
                    )
                if lost:
                    logger.warning(
                        f"ChIP-Atlas: {len(lost)} factor(s) missing from this "
                        f"build after retries: {lost}"
                    )
                return proteome

            proteome = step(output_dir, organism, "proteome_chip_atlas",
                            _add_tf_targets, resume, chain=proteome_chain)
            transnet.proteome = proteome
        
        # Save the network
        timestamp = datetime.now().strftime("%Y%m%d")
        org_dir = os.path.join(output_dir, organism, timestamp)
        os.makedirs(org_dir, exist_ok=True)
        
        logger.info(f"Saving network to {org_dir}")
        transnet.save_network(org_dir)

        problems = _report_build_quality(transnet, organism, brenda, config)
        
        # Create a latest symlink
        latest_dir = os.path.join(output_dir, organism, "latest")
        if os.path.exists(latest_dir):
            if os.path.islink(latest_dir):
                logger.info(f"Removing existing symlink at {latest_dir}")
                os.unlink(latest_dir)
            else:
                logger.warning(f"{latest_dir} exists but is not a symlink. Removing...")
                import shutil
                shutil.rmtree(latest_dir)
        
        # Create relative symlink
        logger.info(f"Creating 'latest' symlink to {timestamp}")
        os.symlink(timestamp, latest_dir, target_is_directory=True)
        
        # Create a simple summary file
        summary = {
            "organism": organism,
            "build_date": timestamp,
            "pathways_count": len(pathways.pathways),
            "reactions_count": len(reactions.reactions),
            "proteins_count": len(proteome.proteins),
            "metabolites_count": len(metabolome.metabolites),
            "genes_count": len(transcriptome.genes) if hasattr(transcriptome, 'genes') else 0
        }
        
        with open(os.path.join(org_dir, "summary.txt"), "w") as f:
            for key, value in summary.items():
                f.write(f"{key}: {value}\n")
        
        if problems:
            logger.warning(f"Built network for {organism}, with gaps")
            return "gaps"
        logger.info(f"Successfully built network for {organism}")
        return True
    
    except Exception as e:
        logger.error(f"Error building network for {organism}: {e}")
        logger.error(traceback.format_exc())
        return False

def main():
    parser = argparse.ArgumentParser(description="Build integrated networks for model organisms")
    parser.add_argument("--organisms", nargs="+", default=["human", "mouse", "yeast", "ecoli"],
                        help="Organisms to build networks for")
    parser.add_argument("--output-dir", default="data",
                        help="Directory to save the networks")
    parser.add_argument("--debug", action="store_true",
                        help="Run in debug mode with smaller network components")
    parser.add_argument("--brenda", action="store_true",
                        help="fetch BRENDA allosteric regulators (needs "
                             "BRENDA_EMAIL / BRENDA_PASSWORD); slow but it is "
                             "what enables the metabolite regulation axis")
    parser.add_argument("--all-proteins", action="store_true",
                        help="include unreviewed TrEMBL entries (mouse: ~88,000 "
                             "proteins instead of ~17,000). Slower and less "
                             "reliable; reviewed-only is the default")
    parser.add_argument("--signaling", action="store_true",
                        help="add a Signaling layer from KEGG signal-transduction "
                             "KGML (kinase-substrate and kinase-TF relations). "
                             "Needed for studies with phosphoproteomics")
    parser.add_argument("--chip-score-threshold", type=float, default=100.0,
                        help="minimum mean ChIP-Atlas binding score for a "
                             "TF-target edge (default: 100). ChIP-Atlas lists "
                             "every gene with any binding, so 0 gives ~15,000 "
                             "targets per factor and transcriptional "
                             "regulation swamps the network")
    parser.add_argument("--no-resume", action="store_true",
                        help="ignore existing checkpoints and rebuild every "
                             "step from scratch")
    parser.add_argument("--clear-checkpoints", action="store_true",
                        help="delete the checkpoints for the selected organisms "
                             "and exit")
    parser.add_argument("--reactions-file", default=None,
                        help="reuse this KEGG reaction table instead of "
                             "re-fetching it (default: data/master_reactions.csv "
                             "if it exists)")
    parser.add_argument("--refresh-reactions", action="store_true",
                        help="re-fetch the KEGG reaction table even if one "
                             "exists; ~12,000 REST calls, so only when KEGG "
                             "reactions have actually changed")
    parser.add_argument("--skip-reactions-file", action="store_true",
                        help="do not use a reaction table at all and query the "
                             "API per reaction (slow; kept for compatibility)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.clear_checkpoints:
        for organism in args.organisms:
            clear_checkpoints(args.output_dir, organism)
        return

    # Resolve the KEGG reaction table. Fetching it is ~12,000 individual REST
    # calls and takes hours, while the table changes rarely -- so an existing
    # one is reused unless a refresh is asked for explicitly.
    reactions_file = None
    if args.skip_reactions_file:
        logger.info("Not using a reaction table; every reaction will be fetched")
    else:
        candidates = [
            args.reactions_file,
            os.path.join(args.output_dir, "master_reactions.csv"),
            os.path.join(REPO_ROOT, "data", "master_reactions.csv"),
        ]
        existing = next(
            (c for c in candidates if c and os.path.exists(c)), None
        )
        if existing and not args.refresh_reactions:
            logger.info(f"Reusing reaction table {existing}")
            logger.info("(pass --refresh-reactions to re-fetch it from KEGG)")
            reactions_file = existing
        else:
            if args.refresh_reactions:
                logger.info("Refreshing the reaction table from KEGG")
            else:
                logger.info("No reaction table found; fetching from KEGG")
            logger.info("This is ~12,000 REST calls and takes a while")
            reactions_file = generate_reactions_file(args.output_dir)
    
    results = {}
    for organism in args.organisms:
        logger.info(f"======= Starting build for {organism} =======")
        success = build_organism_network(
            organism, args.output_dir, reactions_file, args.debug,
            brenda=args.brenda, resume=not args.no_resume,
            reviewed_only=not args.all_proteins,
            chip_score_threshold=args.chip_score_threshold,
            signaling=args.signaling,
        )
        if success == "gaps":
            results[organism] = "Built with gaps (see warnings above)"
        else:
            results[organism] = "Success" if success else "Failed"
        logger.info(f"======= Completed build for {organism}: {results[organism]} =======")
    
    # Print summary
    logger.info("===== Build Summary =====")
    for organism, result in results.items():
        logger.info(f"{organism}: {result}")
    logger.info("========================")

    # A network with gaps used to exit 0 and read "Success", so anything
    # driving this script (run_all.sh) carried on with a broken network.
    if any(result != "Success" for result in results.values()):
        return 1
    return 0

if __name__ == "__main__":
    sys.exit(main())
