"""Hover-help content.

Every entry is grounded in how the value is actually consumed in
``pHinder_7_0.py`` and ``structure/pHinderSurface.py`` -- the comparison it
feeds, the function it is passed to -- rather than inferred from its name.
Where a parameter is set but never read, the help says so instead of inventing
a mechanism for it.

Format: key -> (title, body).
"""

# Depth convention, stated once and referenced by the cutoffs that use it:
# inLocalSurface() compares a sidechain's average distance to the surface,
# signed so that negative values are inside.
_DEPTH = ("Depth is signed against the molecular surface: negative is inside, "
          "positive is outside.")

CALCULATIONS = {
    "topologyCalculation": (
        "Residue network topology",
        "Triangulates the selected residues, prunes edges longer than the max "
        "network edge length, then analyses the networks that remain — tight "
        "bonds and network parity. Writes the triangulation and the pruned network."),
    "surfaceCalculation": (
        "Molecular surface",
        "Builds the pHinder molecular surface for the selected chains, and surfaces "
        "any ligands. Classification, interfaces and screening all need this surface "
        "and will build it themselves if you do not tick it here."),
    "sidechainClassification": (
        "Sidechain classification",
        "Measures how deep each selected sidechain sits relative to the surface and "
        "labels it core, margin, or exposed using the Classification cutoffs. "
        "Needs the triangulation and the surface."),
    "interfaceClassification": (
        "Interface classification",
        "Flags sidechains lying within the interface distance filter of another "
        "chain. Needs the classification, so it runs that first if you have not "
        "ticked it."),
    "virtualScreenSurfacesCalculation": (
        "Virtual screening surfaces",
        "Samples points from the surface — stepping inward and outward by the "
        "iteration settings — discards points that clash with atoms, then "
        "triangulates what is left and parses it into void volumes."),
}

PARAMETERS = {
    # --- Classification ----------------------------------------------------
    "CORE_CUTOFF": (
        "Core cutoff (Å)",
        f"A sidechain counts as core when its average distance to the surface is at "
        f"or below this value. {_DEPTH} The default of −3.0 means at least 3 Å inside. "
        "Make it more negative to call fewer residues core."),
    "MARGIN_CUTOFF": (
        "Margin cutoff (Å)",
        f"The boundary between margin and exposed, applied after the core test. "
        f"{_DEPTH} Anything deeper than this but shallower than the core cutoff is "
        "margin; anything shallower is exposed."),

    # --- Networks ----------------------------------------------------------
    "MAX_NETWORK_EDGE_LENGTH": (
        "Max network edge length (Å)",
        "Edges longer than this are pruned from the residue triangulation before "
        "the networks are analysed. Raising it keeps more distant residues connected "
        "and yields larger, sparser networks."),
    "MIN_NETWORK_SIZE": (
        "Min network size (residues)",
        "Networks smaller than this are dropped when the networks are ranked from "
        "largest to smallest. The default of 1 keeps every network, including "
        "isolated residues."),
    "SAVE_NETWORK_TRIANGULATION": (
        "Save network triangulation",
        "Write the triangulation itself alongside the network results, not just the "
        "pruned network."),

    # --- Surfaces ----------------------------------------------------------
    "circumSphereRadiusLimit": (
        "Circumsphere radius limit (Å)",
        "Upper bound on the circumsphere radius of a simplex kept when the surface "
        "is built — in effect, how tightly the surface is allowed to wrap the "
        "structure. Passed straight to the surface calculation."),
    "minArea": (
        "Minimum facet area (Å²)",
        "Facets smaller than this are not kept when the surface is assembled. "
        "Passed straight to the surface calculation."),
    "HIGH_RESOLUTION_SURFACE": (
        "High-resolution surface",
        "Build the surface at higher resolution. Slower, and worth it when you need "
        "fine detail in the depth measurements that classification depends on."),
    "SAVE_SURFACE": (
        "Save surface",
        "Write the computed molecular surface to the results directory."),
    "ALLOW_SMALL_SURFACES": (
        "Allow small surfaces",
        "Permit surfaces that would otherwise be rejected as too small. Useful for "
        "short chains, peptides, and fragments that the normal test discards."),
    "SAVE_LIGAND_SURFACES": (
        "Save ligand surfaces",
        "Write the surfaces computed for ligands as well as for the protein."),
    "WRITE_SURFACE_CREATION_ANIMATION": (
        "Write surface creation animation",
        "Emit the intermediate states of surface construction so the build can be "
        "played back. Diagnostic, and it writes a lot of files."),

    # --- Interfaces --------------------------------------------------------
    "INTERFACE_DISTANCE_FILTER": (
        "Interface distance filter (Å)",
        "A sidechain is treated as interface when its distance to another chain is "
        "strictly less than this. Only meaningful with more than one chain selected "
        "on the Input tab."),

    # --- Screening ---------------------------------------------------------
    "MAX_VOID_NETWORK_EDGE_LENGTH": (
        "Max void network edge length (Å)",
        "Edges longer than this are pruned when sampling points are triangulated "
        "into voids. Smaller values split the sampled volume into more, tighter "
        "voids; larger values merge them."),
    "MIN_VOID_NETWORK_SIZE": (
        "Min void network size (points)",
        "Voids built from fewer sampling points than this are discarded — the "
        "control that separates a real pocket from a handful of stray points."),
    "VIRTUAL_CLASH_CUTOFF": (
        "Virtual clash cutoff (Å)",
        "A sampling point is thrown away when it lies within this distance of any "
        "atom. Larger values carve more space away from the protein and leave "
        "smaller voids."),
    "IN_ITERATIONS": (
        "Inward iterations",
        "How many steps of sampling points are generated inward from each surface "
        "facet. More iterations reach deeper into the structure, at more cost."),
    "IN_ITERATIONS_STEP_SIZE": (
        "Inward step size (Å)",
        "How far each inward step moves. Step size times iterations sets how deep "
        "the sampling reaches."),
    "OUT_ITERATIONS": (
        "Outward iterations",
        "How many steps of sampling points are generated outward from each surface "
        "facet, into the solvent."),
    "OUT_ITERATIONS_STEP_SIZE": (
        "Outward step size (Å)",
        "How far each outward step moves. Step size times iterations sets how far "
        "from the surface the sampling extends."),

    # --- Advanced ----------------------------------------------------------
    "INCLUDE_HYDROGENS": (
        "Include hydrogens",
        "Keep hydrogen atoms from the structure instead of discarding them. Most "
        "crystal structures have none; models and NMR structures do."),
    "INCLUDE_WATER": (
        "Include water",
        "Treat crystallographic waters as part of the structure rather than "
        "stripping them."),
    "INCLUDE_IONS": (
        "Include ions",
        "Treat ions as part of the structure rather than stripping them."),
    "SAVE_LOG_FILE": (
        "Save log file",
        "Write pHinder_parameters.log into the results directory, recording every "
        "parameter this run used. Worth keeping — it is what makes a run reproducible."),
    "PYTHON_RECURSION_LIMIT": (
        "Python recursion limit",
        "The network traversal recurses, and large structures can exceed Python's "
        "default limit of 1000. Raised for the run and restored afterwards. Raise it "
        "further if a big structure stops with a recursion error."),
}

INPUT = {
    "file_path": (
        "Structure file",
        "The PDB or mmCIF file to analyse. Chain checkboxes appear once a file is "
        "chosen. Gzipped files are handled."),
    "save_path": (
        "Save location",
        "Results are written to a pHinderResults folder inside this directory. "
        "Left empty, they go beside the structure file."),
    "chains": (
        "Chains",
        "Which chains to include. Interface classification needs at least two."),
    "group_chains": (
        "Group chains",
        "Treat the selected chains as one structure rather than analysing each "
        "chain separately."),
    "residues": (
        "Residues",
        "Which sidechains are triangulated and classified. pHinder defaults to the "
        "ionizable set — D, E, K, R, H. Selecting a set that matches one of "
        "pHinder's named sets uses that name; any other combination becomes a "
        "custom set."),
}

ACTIONS = {
    "run": ("Run pHinder",
            "Runs the ticked calculations in order. Shared work — the triangulation, "
            "the surface, the classification — is done once and reused."),
    "stop": ("Stop",
             "Stops after the pHinder call in progress finishes. A single long step "
             "cannot be interrupted part-way, so this is not immediate."),
    "clear": ("Clear", "Empties the output pane. Does not affect a running calculation."),
}


# Which options are genuinely on/off. Decided per parameter, because the value
# alone cannot tell you: IN_ITERATIONS and MIN_NETWORK_SIZE both default to 1
# but are counts, and rendering them as checkboxes capped them at 1.
BOOLEAN_OPTIONS = {
    "SAVE_NETWORK_TRIANGULATION",
    "HIGH_RESOLUTION_SURFACE", "SAVE_SURFACE", "ALLOW_SMALL_SURFACES",
    "SAVE_LIGAND_SURFACES", "WRITE_SURFACE_CREATION_ANIMATION",
    "INCLUDE_HYDROGENS", "INCLUDE_WATER",
    "INCLUDE_INS_PLACEHOLDER", "INCLUDE_IONS", "SAVE_LOG_FILE",
}


def is_boolean(key):
    return key in BOOLEAN_OPTIONS


def for_option(key):
    """(title, body) for a parameter, or (None, None) if undocumented."""
    return PARAMETERS.get(key, (None, None))
