import argparse
import re
import logging
import multiprocessing as mp
import time
import sys
from os import sep
from pHinder.pHinder_7_0 import pHinder

def parse_chains(value: str):
    # Normalize separators: space, comma, or dot
    # Examples handled:
    # "A B C"   -> ["A", "B", "C"]
    # "A,B,C"   -> ["A", "B", "C"]
    # "A.B.C"   -> ["A", "B", "C"]
    # "ABC"     -> ["A", "B", "C"]
    # Also works without quotes if the shell doesn’t split: A,B,C  or A.B.C
    # But: unquoted `A B C` will already be 3 argv tokens, see note below.
    tokens = re.split(r"[,\.\s]+", value.strip())
    # If user gave "ABC" with no separators, tokens will be ["ABC"],
    # so split that into characters.
    if len(tokens) == 1:
        tokens = list(tokens[0])
    # Filter out any empty strings
    return [t for t in tokens if t]

def normalize_residues(raw_tokens):
    # raw_tokens is e.g. ["D", "E", "K", "R", "H"]
    # or ["DEKRH"] or ["D,E,K,R,H"] or ["D.E.K.R.H"] ...
    raw = " ".join(raw_tokens)
    # Split on comma, dot, or whitespace
    tokens = re.split(r"[,\.\s]+", raw.strip())
    if len(tokens) == 1:
        # Single token like "DEKRH" -> list of characters
        tokens = list(tokens[0])
    return [t for t in tokens if t]

def parse_arguments():
    parser = argparse.ArgumentParser(description="Run pHinder from the command line.")

    # Required inputs
    parser.add_argument("pdb_file", help="Path to input PDB or mmCIF file")
    # parser.add_argument("--chains", nargs="+", required=True, help="List of chain IDs to include")
    parser.add_argument(
        "--chains",
        required=True,
        type=parse_chains,
        help="Chains as ABC, A B C, A,B,C, A.B.C, etc."
    )
    parser.add_argument("--save_path", default=".", help="Directory to save results")

    # Residue selection options
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--residue-set", choices=[
        "allSet", "apolarSet", "polarSet", "ionizableSet",
        "ionizableSetNoCys", "acidicSet", "basicSet"
    ], default="ionizableSetNoCys", help="Use a predefined amino acid group (default: ionizableSetNoCys)")
    # group.add_argument("--residues", nargs="+", help="List of amino acids to use (e.g., D E K R H)")
    group.add_argument(
        "--residues",
        type=normalize_residues,
        help="AAs as D E K R H, DEKRH, D,E,K,R,H, D.E.K.R.H, etc.")

    # Calculation options
    parser.add_argument("--topology-calculation", action="store_true", default=True, help="Enable topology calculation (default: True)")
    parser.add_argument("--surface-calculation", action="store_true", default=False, help="Enable surface calculation (default: False)")
    parser.add_argument("--sidechain-classification", action="store_true", default=False, help="Enable sidechain classification (default: False)")
    parser.add_argument("--interface-classification", action="store_true", default=False, help="Enable interface classification (default: False)")
    parser.add_argument("--interface-distance-filter", type=float, default=8.0, help="Optional interface distance filter for interface classification")
    parser.add_argument("--virtual-screen-surfaces-calculation", action="store_true", default=False, help="Enable virtual screen surfaces calculation (default: False)")

    # Sidechain classification
    parser.add_argument("--core-cutoff", type=float, default=-3.0, help="Set CORE_CUTOFF (default: -3.0)")
    parser.add_argument("--margin-cutoff", type=float, default=-2.0, help="Set MARGIN_CUTOFF (default: -2.0)")
    parser.add_argument("--margin-cutoff-core-network", type=float, default=-2.0, help="Set MARGIN_CUTOFF_CORE_NETWORK (default: -2.0)")

    # Network options
    parser.add_argument("--max-network-edge-length", type=float, default=10.0, help="Set MAX_NETWORK_EDGE_LENGTH (default: 10.0)")
    parser.add_argument("--min-network-size", type=int, default=1, help="Set MIN_NETWORK_SIZE (default: 1)")
    parser.add_argument("--reduced-network-representation", type=int, default=1, help="Enable REDUCED_NETWORK_REPRESENTATION (default: 1)")
    parser.add_argument("--save-network-triangulation", type=int, default=1, help="Enable SAVE_NETWORK_TRIANGULATION (default: 1)")

    # Surface options
    parser.add_argument("--circumsphereradiuslimit", type=float, default=6.5, help="Set circumSphereRadiusLimit (default: 6.5)")
    parser.add_argument("--minarea", type=int, default=10, help="Set minArea (default: 10)")
    parser.add_argument("--high-resolution-surface", type=int, default=1, help="Enable HIGH_RESOLUTION_SURFACE (default: 1)")
    parser.add_argument("--save-surface", type=int, default=1, help="Enable SAVE_SURFACE (default: 1)")
    parser.add_argument("--allow-small-surfaces", type=int, default=0, help="Enable ALLOW_SMALL_SURFACES (default: 0)")
    parser.add_argument("--save-ligand-surfaces", type=int, default=0, help="Enable SAVE_LIGAND_SURFACES (default: 0)")
    parser.add_argument("--write-surface-creation-animation", type=int, default=0, help="Enable WRITE_SURFACE_CREATION_ANIMATION (default: 0)")

    # Virtual screening options
    parser.add_argument("--max-void-network-edge-length", type=float, default=2.0, help="Set MAX_VOID_NETWORK_EDGE_LENGTH (default: 2.0)")
    parser.add_argument("--min-void-network-size", type=int, default=10, help="Set MIN_VOID_NETWORK_SIZE (default: 10)")
    parser.add_argument("--virtual-clash-cutoff", type=float, default=2.5, help="Set VIRTUAL_CLASH_CUTOFF (default: 2.5)")
    parser.add_argument("--in-iterations", type=int, default=1, help="Set IN_ITERATIONS (default: 1)")
    parser.add_argument("--in-iterations-step-size", type=float, default=2.0, help="Set IN_ITERATIONS_STEP_SIZE (default: 2.0)")
    parser.add_argument("--out-iterations", type=int, default=1, help="Set OUT_ITERATIONS (default: 1)")
    parser.add_argument("--out-iterations-step-size", type=float, default=2.0, help="Set OUT_ITERATIONS_STEP_SIZE (default: 2.0)")

    # General options
    parser.add_argument("--group-chains", type=int, default=0, help="Group chains together for analysis (default: 0)")

    # Advanced
    parser.add_argument("--allow-cys-core-seeding", type=int, default=0, help="Enable ALLOW_CYS_CORE_SEEDING (default: 0)")
    parser.add_argument("--include-hydrogens", type=int, default=0, help="Enable INCLUDE_HYDROGENS (default: 0)")
    parser.add_argument("--include-water", type=int, default=0, help="Enable INCLUDE_WATER (default: 0)")
    parser.add_argument("--include-ions", type=int, default=0, help="Enable INCLUDE_IONS (default: 0)")
    parser.add_argument("--save-log-file", type=int, default=1, help="Enable SAVE_LOG_FILE (default: 1)")
    parser.add_argument("--python-recursion-limit", type=int, default=10000, help="Set PYTHON_RECURSION_LIMIT (default: 10000)")

    return parser.parse_args()

def write_phinder_log(p, log_path):
    log_lines = [
        "=======================================",
        "       pHinder Parameter Settings",
        "======================================="
    ]

    # Create a list of attributes you want to record
    attr_names = [
        "gui", "processes", "pdbFilePath", "pdbFileName", "outPath", "pdbFormat", "zip", "chains", "group_chains",
        "residueSet", "maxNetworkEdgeLength", "minNetworkSize", "reducedNetworkRepresentation",
        "saveNetworkTriangulation", "highResolutionSurface", "saveSurface", "allowSmallSurfaces",
        "saveLigandSurfaces", "writeSurfaceCreationAnimation", "coreCutoff", "marginCutoff",
        "marginCutoffCoreNetwork", "interface_distance_filter", "virtualClashCutoff",
        "inIterations", "inIterationsStepSize", "outIterations", "outIterationsStepSize",
        "allowCysCoreSeeding", "includeHydrogens", "includeWater", "includeIons"
    ]

    max_len = max(len(attr) for attr in attr_names)
    for attr in attr_names:
        value = getattr(p, attr, "<not set>")
        log_lines.append(f"{attr.ljust(max_len)} : {value}")

    log_lines.append("=======================================")

    # Write to file
    with open(log_path, "w") as f:
        f.write("\n".join(log_lines) + "\n")

    logging.info(f"pHinder attribute log saved to {log_path}")


def main():
    args = parse_arguments()

    logging.basicConfig(filename="phinder_run.log", level=logging.INFO)

    n_processes = max(1, mp.cpu_count() - 1)
    pdb_file_path = sep.join(args.pdb_file.split(sep)[:-1]) + sep
    pdb_file_name = args.pdb_file.split(sep)[-1]
    zip_status = 1 if ".gz" in pdb_file_name else 0
    if pdb_file_name.endswith('.cif') or pdb_file_name.endswith('.cif.gz'):
        pdb_format = "mmCIF"
    elif pdb_file_name.endswith('.pdb') or pdb_file_name.endswith('.pdb.gz'):
        pdb_format = "pdb"
    else:
        logging.warning(f"Unknown file extension for {pdb_file_name}, defaulting to PDB format")
        pdb_format = "pdb"
    save_path = args.save_path
    if not save_path or save_path.strip() == ".":
        save_path = pdb_file_path
    if not save_path.endswith(sep):
        save_path += sep
    logging.info(f"Resolved save_path: {save_path}")

    # Determine residue set
    if args.residues:
        residue_set = "customSet:" + ",".join([r.upper() for r in args.residues])
    else:
        residue_set = args.residue_set

    p = pHinder()
    p.gui = False
    p.processes = n_processes
    p.pdbFilePath = pdb_file_path
    p.pdbFileName = pdb_file_name
    p.outPath = save_path + "pHinderResults/"
    args.save_path = p.outPath
    p.pdbFormat = pdb_format
    p.zip = zip_status
    p.chains = args.chains
    p.group_chains = args.group_chains
    logging.info(f"Group chains option set to: {args.group_chains}")
    p.residueSet = residue_set

    # Set options from args
    p.maxNetworkEdgeLength = args.max_network_edge_length
    p.minNetworkSize = args.min_network_size
    p.reducedNetworkRepresentation = args.reduced_network_representation
    p.saveNetworkTriangulation = args.save_network_triangulation
    p.highResolutionSurface = args.high_resolution_surface
    p.saveSurface = args.save_surface
    p.allowSmallSurfaces = args.allow_small_surfaces
    p.saveLigandSurfaces = args.save_ligand_surfaces
    p.writeSurfaceCreationAnimation = args.write_surface_creation_animation
    p.coreCutoff = args.core_cutoff
    p.marginCutoff = args.margin_cutoff
    p.marginCutoffCoreNetwork = args.margin_cutoff_core_network
    p.interface_distance_filter = args.interface_distance_filter
    p.virtualClashCutoff = args.virtual_clash_cutoff
    p.inIterations = args.in_iterations
    p.inIterationsStepSize = args.in_iterations_step_size
    p.outIterations = args.out_iterations
    p.outIterationsStepSize = args.out_iterations_step_size
    p.allowCysCoreSeeding = args.allow_cys_core_seeding
    p.includeHydrogens = args.include_hydrogens
    p.includeWater = args.include_water
    p.includeIons = args.include_ions

    sys.setrecursionlimit(args.python_recursion_limit)

    start = time.time()
    p.setQuerySet()
    p.openPDBs(p.pdbFilePath, p.pdbFileName, zip_status=p.zip)
    p.hetLigand4D()
    p.hydrogens()
    p.makeAtomCollections()
    p.makeVertices4D()

    if args.topology_calculation:
        p.selectTscTriangulationAtoms()
        p.triangulateTscAtoms()
        p.writeTriangulation()
        p.pruneTriangulation()
        p.minimizePrunedTriangulation()
        p.identifyTightBonds()
        p.calculateNetworkParity()

    if args.surface_calculation:
        p.surface(args.circumsphereradiuslimit, args.minarea)
        p.writeSurface()
        p.surfaceLigands()
        p.writeLigandSurfaces()

    if args.sidechain_classification:
        if not args.surface_calculation:
            logging.info("Auto-running surface calculation for sidechain classification")
            p.surface(args.circumsphereradiuslimit, args.minarea)
            p.writeSurface()
        p.selectTscClassificationAtoms()
        p.classifySidechainLocation()
        p.identifyMissingTscAtoms()
        p.writeSidechainClassificationResults()

    if args.interface_classification:
        if not args.surface_calculation:
            logging.info("Auto-running surface calculation for interface classification")
            p.surface(args.circumsphereradiuslimit, args.minarea)
            p.writeSurface()

        if not args.sidechain_classification:
            logging.info("Auto-running sidechain classification for interface classification")
            p.selectTscClassificationAtoms()
            p.classifySidechainLocation()
            p.identifyMissingTscAtoms()
            p.writeSidechainClassificationResults()
        p.classifyInterfaceSidechains()

    if args.virtual_screen_surfaces_calculation:
        if not args.surface_calculation:
            logging.info("Auto-running surface calculation for interface classification")
            p.surface(args.circumsphereradiuslimit, args.minarea)
            p.writeSurface()

        p.makeSamplingGridUsingProteinSurface()
        p.filterSamplingPointsUsingClashes()
        p.triangulateRemainingGridPoints()
        p.identifyAndParseIndividualSamplingVoids(
            maxVoidNetworkEdgeLength=args.max_void_network_edge_length,
            minVoidNetworkEdgeLength=0.0,
            minVoidNetworkSize=args.min_void_network_size,
            psa=1
        )
        p.calculateSamplingVoidSurfaces(extend_sampling=True)

    print(f"pHinder completed in {time.time() - start:.2f} seconds.")

    # Write log file
    param_log_path = p.outPath + "pHinder_parameters.log"
    write_phinder_log(p, param_log_path)


if __name__ == "__main__":
    main()
