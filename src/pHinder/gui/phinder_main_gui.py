import tkinter as tk
from tkinter import ttk, messagebox, IntVar
from os.path import sep
from pHinder.gui.dynamic_option_widget import DynamicOptionWidget
from pHinder.gui.input_widgets import IntegerInputWidget, FloatInputWidget, StringInputWidget
from pHinder.gui.file_open import FilePathWidget
from pHinder.gui.dynamic_option_widget_amino_acid_selection import AminoAcidSelectionWidget
from pHinder.gui.collapsible_toggle import CollapsiblePane
from pHinder.gui.terminal_output_widget import TerminalOutputWidget  # Import terminal output widget
import io
import threading
import os
import logging

# Configure logging at the beginning of the script
logging.basicConfig(
    filename="phinder_run.log",  # Log file name
    level=logging.INFO,  # Logging level
    format="%(asctime)s - %(levelname)s - %(message)s"  # Log format
)

class TerminalBuffer:
    """
    A class to capture terminal output dynamically and display it in a tkinter widget.

    This class integrates with the TerminalOutputWidget by using a StringIO buffer to
    temporarily store output data. It periodically flushes the data into the widget to
    ensure that updates to the terminal are reflected in real-time.

    The buffer is reset after each flush to prevent redundant data from being displayed
    multiple times. Ensure flush is called frequently enough to avoid losing unprocessed
    data that may accumulate between updates.
    """
    def __init__(self, terminal_output_widget):
        self.buffer = io.StringIO()
        self.terminal_output_widget = terminal_output_widget

    def write(self, message):
        self.buffer.write(message)

    def flush(self):
        self.buffer.seek(0)
        data = self.buffer.read()
        if data.strip():  # Avoid blank updates
            self.terminal_output_widget.write(data)
        self.buffer = io.StringIO()  # Reset the buffer

# Function to retrieve PDB chains from a file
def get_pdb_chains(path):
    try:
        pdb_file_path = sep.join(path.split(sep)[:-1]) + sep
        pdb_file_name = path.split(sep)[-1]
        # Import appropriate PDB file handler
        if ".cif" in pdb_file_name:
            from pHinder._vendor.pdbFile import PDBfile 
        else:
            from pHinder._vendor.pdbFile import PDBfile
        # Check if the file is compressed
        if ".gz" in pdb_file_name:
            pdb = PDBfile(pdb_file_path, pdb_file_name, zip_status=1)
        else:
            pdb = PDBfile(pdb_file_path, pdb_file_name, zip_status=0)
        # Create IntVar variables for each chain
        chain_vars = {chain: IntVar() for chain in pdb.res_chains}
        return chain_vars
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load chains: {str(e)}")

# Helper function to align collapsible panes
def align_collapsible_pane(pane):
    pane.pack(pady=10, anchor="w")

# RedirectOutput class for capturing stdout and stderr
class RedirectOutput:
    def __init__(self, terminal_output):
        self.terminal_output = terminal_output

    def write(self, message):
        if message.strip():  # Avoid blank lines
            self.terminal_output.write(message)

    def flush(self):
        pass

# Default calculation options and variables
# These options control the calculation behavior
topologyCalculation = 0
surfaceCalculation = 0
sidechainClassification = 0
interfaceClassification = 0
virtualScreenSurfacesCalculation = 0

calculation_options = {
    "topologyCalculation": topologyCalculation,
    "surfaceCalculation": surfaceCalculation,
    "sidechainClassification": sidechainClassification,
    "interfaceClassification": interfaceClassification,
    "virtualScreenSurfacesCalculation": virtualScreenSurfacesCalculation
}

default_calculation_options = {
    "topologyCalculation": 0,
    "surfaceCalculation": 0,
    "sidechainClassification": 0,
    "interfaceClassification": 0,
    "virtualScreenSurfacesCalculation": 0
}

# Default parameters for sidechain classification
CORE_CUTOFF = -3.0
MARGIN_CUTOFF = -2.0
MARGIN_CUTOFF_CORE_NETWORK = -2.0

sidechain_classification_options = {
    "CORE_CUTOFF": CORE_CUTOFF,
    "MARGIN_CUTOFF": MARGIN_CUTOFF
}

default_sidechain_classification_options = {
    "CORE_CUTOFF": -3.0,
    "MARGIN_CUTOFF": -2.0
}

# Default network options
MAX_NETWORK_EDGE_LENGTH = 10.0
MIN_NETWORK_SIZE = 1
REDUCED_NETWORK_REPRESENTATION = 1
SAVE_NETWORK_TRIANGULATION = 1

network_options = {
    "MAX_NETWORK_EDGE_LENGTH": MAX_NETWORK_EDGE_LENGTH,
    "MIN_NETWORK_SIZE": MIN_NETWORK_SIZE,
    "SAVE_NETWORK_TRIANGULATION": SAVE_NETWORK_TRIANGULATION
}

default_network_options = {
    "MAX_NETWORK_EDGE_LENGTH": 10.0,
    "MIN_NETWORK_SIZE": 1,
    "SAVE_NETWORK_TRIANGULATION": 1
}

# Default surface parameters
circumSphereRadiusLimit = 6.5
minArea = 10
HIGH_RESOLUTION_SURFACE = 1
SAVE_SURFACE = 1
ALLOW_SMALL_SURFACES = 1
SAVE_LIGAND_SURFACES = 0
WRITE_SURFACE_CREATION_ANIMATION = 0

surface_options = {
    "circumSphereRadiusLimit": circumSphereRadiusLimit,
    "minArea": minArea,
    "HIGH_RESOLUTION_SURFACE": HIGH_RESOLUTION_SURFACE,
    "SAVE_SURFACE": SAVE_SURFACE,
    "ALLOW_SMALL_SURFACES": ALLOW_SMALL_SURFACES,
    "SAVE_LIGAND_SURFACES": SAVE_LIGAND_SURFACES,
    "WRITE_SURFACE_CREATION_ANIMATION": WRITE_SURFACE_CREATION_ANIMATION
}

default_surface_options = {
    "circumSphereRadiusLimit": 6.5,
    "minArea": 10,
    "HIGH_RESOLUTION_SURFACE": 1,
    "SAVE_SURFACE": 1,
    "ALLOW_SMALL_SURFACES": 1,
    "SAVE_LIGAND_SURFACES": 0,
    "WRITE_SURFACE_CREATION_ANIMATION": 0
}

# Default interface parameters
INTERFACE_DISTANCE_FILTER = 8.0

# Interface option
interface_options = {"INTERFACE_DISTANCE_FILTER": INTERFACE_DISTANCE_FILTER}

default_interface_options = {"INTERFACE_DISTANCE_FILTER":8.0}

# Virtual screening options
MAX_VOID_NETWORK_EDGE_LENGTH = 2.0
MIN_VOID_NETWORK_SIZE = 10
VIRTUAL_CLASH_CUTOFF = 2.5
IN_ITERATIONS = 1
IN_ITERATIONS_STEP_SIZE = 2.0
OUT_ITERATIONS = 1
OUT_ITERATIONS_STEP_SIZE = 2.0

virtual_screening_options = {
    "MAX_VOID_NETWORK_EDGE_LENGTH": MAX_VOID_NETWORK_EDGE_LENGTH,
    "MIN_VOID_NETWORK_SIZE": MIN_VOID_NETWORK_SIZE,
    "VIRTUAL_CLASH_CUTOFF": VIRTUAL_CLASH_CUTOFF,
    "IN_ITERATIONS": IN_ITERATIONS,
    "IN_ITERATIONS_STEP_SIZE": IN_ITERATIONS_STEP_SIZE,
    "OUT_ITERATIONS": OUT_ITERATIONS,
    "OUT_ITERATIONS_STEP_SIZE": OUT_ITERATIONS_STEP_SIZE
}

default_virtual_screening_options = {
    "MAX_VOID_NETWORK_EDGE_LENGTH": 2.0,
    "MIN_VOID_NETWORK_SIZE": 10,
    "VIRTUAL_CLASH_CUTOFF": 2.5,
    "IN_ITERATIONS": 1,
    "IN_ITERATIONS_STEP_SIZE": 2.0,
    "OUT_ITERATIONS": 1,
    "OUT_ITERATIONS_STEP_SIZE": 2.0
}

# Advanced options for calculations
ALLOW_CYS_CORE_SEEDING = 0
INCLUDE_HYDROGENS = 0
INCLUDE_WATER = 0
INCLUDE_IONS = 0
SAVE_LOG_FILE = 1
PYTHON_RECURSION_LIMIT = 10000

advanced_options = {
    "INCLUDE_HYDROGENS": INCLUDE_HYDROGENS,
    "INCLUDE_WATER": INCLUDE_WATER,
    "INCLUDE_IONS": INCLUDE_IONS,
    "SAVE_LOG_FILE": SAVE_LOG_FILE,
    "PYTHON_RECURSION_LIMIT": PYTHON_RECURSION_LIMIT
}

default_advanced_options = {
    "INCLUDE_HYDROGENS": 0,
    "INCLUDE_WATER": 0,
    "INCLUDE_IONS": 0,
    "SAVE_LOG_FILE": 1,
    "PYTHON_RECURSION_LIMIT": 10000
}

# Function to create a collapsible pane for calculation options
def create_calculation_options_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Calculation Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_calculation_options.items()):
        variable = tk.IntVar(value=default_value)
        option_variables[option] = variable
        checkbutton = ttk.Checkbutton(
            collapsible_pane.container, text=option, variable=variable
        )
        checkbutton.grid(row=row, column=0, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

# Function to create a collapsible pane for sidechain classification options
def create_sidechain_classification_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Sidechain Classification Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_sidechain_classification_options.items()):
        variable = tk.DoubleVar(value=default_value)
        option_variables[option] = variable
        label = ttk.Label(collapsible_pane.container, text=option)
        label.grid(row=row, column=0, padx=5, pady=5, sticky="w")
        entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
        entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

# Function to create a collapsible pane for network options
def create_network_options_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Network Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_network_options.items()):
        if option in ["REDUCED_NETWORK_REPRESENTATION", "SAVE_NETWORK_TRIANGULATION", "MIN_NETWORK_SIZE"]:
            variable = tk.IntVar(value=default_value)
            option_variables[option] = variable
            checkbutton = ttk.Checkbutton(
                collapsible_pane.container, text=option, variable=variable
            )
            checkbutton.grid(row=row, column=0, padx=5, pady=5, sticky="w")
        else:
            variable = tk.DoubleVar(value=default_value)
            option_variables[option] = variable
            label = ttk.Label(collapsible_pane.container, text=option)
            label.grid(row=row, column=0, padx=5, pady=5, sticky="w")
            entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
            entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

# Function to create a collapsible pane for surface options
def create_surface_options_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Surface Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_surface_options.items()):
        if default_value in [0, 1]:
            variable = tk.IntVar(value=default_value)
            option_variables[option] = variable
            checkbutton = ttk.Checkbutton(
                collapsible_pane.container, text=option, variable=variable
            )
            checkbutton.grid(row=row, column=0, padx=5, pady=5, sticky="w")
        else:
            variable = tk.DoubleVar(value=default_value)
            option_variables[option] = variable
            label = ttk.Label(collapsible_pane.container, text=option)
            label.grid(row=row, column=0, padx=5, pady=5, sticky="w")
            entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
            entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

# Function to create a collapsible pane for interface options
def create_interface_options_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Interface Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_interface_options.items()):
        variable = tk.DoubleVar(value=default_value)
        option_variables[option] = variable
        label = ttk.Label(collapsible_pane.container, text=option)
        label.grid(row=row, column=0, padx=5, pady=5, sticky="w")
        entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
        entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

# Function to create a collapsible pane for virtual screening options
def create_virtual_screening_options_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Virtual Screening Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_virtual_screening_options.items()):
        if option in ["IN_ITERATIONS", "OUT_ITERATIONS", "MIN_VOID_NETWORK_SIZE"]:
            variable = tk.IntVar(value=default_value)
            option_variables[option] = variable
            label = ttk.Label(collapsible_pane.container, text=option)
            label.grid(row=row, column=0, padx=5, pady=5, sticky="w")
            entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
            entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")
        else:
            variable = tk.DoubleVar(value=default_value)
            option_variables[option] = variable
            label = ttk.Label(collapsible_pane.container, text=option)
            label.grid(row=row, column=0, padx=5, pady=5, sticky="w")
            entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
            entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

# Function to create a collapsible pane for advanced options
def create_advanced_options_pane(parent):
    collapsible_pane = CollapsiblePane(parent, title="Advanced Options")
    option_variables = {}

    for row, (option, default_value) in enumerate(default_advanced_options.items()):
        if default_value in [0, 1]:
            variable = tk.IntVar(value=default_value)
            option_variables[option] = variable

            checkbutton = ttk.Checkbutton(
                collapsible_pane.container,
                text=option,
                variable=variable
            )
            checkbutton.grid(row=row, column=0, padx=5, pady=5, sticky="w")

        else:
            variable = tk.DoubleVar(value=default_value)
            option_variables[option] = variable

            label = ttk.Label(collapsible_pane.container, text=option)
            label.grid(row=row, column=0, padx=5, pady=5, sticky="w")

            entry = ttk.Entry(collapsible_pane.container, textvariable=variable)
            entry.grid(row=row, column=1, padx=5, pady=5, sticky="w")

    return collapsible_pane, option_variables

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

# Main function to build and display the GUI
def main():
    root = tk.Tk()
    root.title("pHinder by Isom Lab")

    # Set the initial width of the window
    root.geometry("1600x800")

    # Ensure the window opens at the top level and gains focus
    root.attributes("-topmost", True)  # Keep the window on top
    root.update_idletasks()  # Ensure changes take effect
    root.attributes("-topmost", False)  # Allow other windows to come on top later
    root.focus_force()  # Force focus on this window

    # Create main frame
    main_frame = tk.Frame(root)
    main_frame.pack(fill="both", expand=True)

    # # Create left frame for widgets
    # left_frame = tk.Frame(main_frame)
    # left_frame.pack(side="left", fill="y", padx=10, pady=10)

    # Create left frame for widgets with a fixed width
    left_frame = tk.Frame(main_frame, width=700)
    left_frame.pack(side="left", fill="y", padx=10, pady=10)
    left_frame.pack_propagate(False)  # Prevent the frame from resizing to its content

    # Create right frame for terminal output
    right_frame = tk.Frame(main_frame)
    right_frame.pack(side="right", fill="both", expand=True, padx=10, pady=10)

    # Add FilePathWidget
    extra_options = {"Group Chains":0}
    file_path_widget = FilePathWidget(left_frame, process_file=get_pdb_chains, options_label="Select Chains:", extra_options=extra_options)
    file_path_widget.frame.pack(pady=10, anchor="w")

    # Add AminoAcidSelectionWidget wrapped in a CollapsiblePane
    collapsible_amino_widget = CollapsiblePane(left_frame, title="Amino Acid Selection")
    collapsible_amino_widget.pack(pady=10, anchor="w")

    default_aa_selections = ["Aspartic Acid (D)", "Glutamic Acid (E)", "Lysine (K)", "Arginine (R)", "Histidine (H)"] # The pHinder set of ionizable groups
# [
#             "Alanine (A)", "Valine (V)", "Leucine (L)", "Isoleucine (I)", "Methionine (M)", "Phenylalanine (F)", "Tryptophan (W)", "Proline (P)",
#             "Serine (S)", "Threonine (T)", "Asparagine (N)", "Glutamine (Q)", "Tyrosine (Y)", "Cysteine (C)",
#             "Aspartic Acid (D)", "Glutamic Acid (E)", "Lysine (K)", "Arginine (R)", "Histidine (H)"
#         ]
    amino_acid_widget = AminoAcidSelectionWidget(collapsible_amino_widget.container, "Amino Acid Options", default_selections=default_aa_selections)
    amino_acid_widget.frame.pack(pady=5, anchor="w")

    # Add Calculation Options Pane
    calculation_options_pane, calculation_variables = create_calculation_options_pane(left_frame)
    calculation_options_pane.pack(pady=10, anchor="w")

    # Add Sidechain Classification Options Pane
    sidechain_classification_pane, sidechain_variables = create_sidechain_classification_pane(left_frame)
    sidechain_classification_pane.pack(pady=10, anchor="w")

    # Add Network Options Pane
    network_options_pane, network_variables = create_network_options_pane(left_frame)
    network_options_pane.pack(pady=10, anchor="w")

    # Add Surface Options Pane
    surface_options_pane, surface_variables = create_surface_options_pane(left_frame)
    surface_options_pane.pack(pady=10, anchor="w")

    # Add Interface Options Pane
    interface_options_pane, interface_variables = create_interface_options_pane(left_frame)
    interface_options_pane.pack(pady=10, anchor="w")

    # Add Virtual Screening Options Pane
    virtual_screening_options_pane, virtual_screening_variables = create_virtual_screening_options_pane(left_frame)
    virtual_screening_options_pane.pack(pady=10, anchor="w")

    # Add Advanced Options Pane
    advanced_options_pane, advanced_variables = create_advanced_options_pane(left_frame)
    advanced_options_pane.pack(pady=10, anchor="w")

    # Add Terminal Output Widget
    terminal_output = TerminalOutputWidget(right_frame, title="pHinder Output")
    terminal_output.pack(fill="both", expand=True)

    # Redirect sys.stdout and sys.stderr to the terminal buffer
    import sys
    terminal_buffer = TerminalBuffer(terminal_output)
    sys.stdout = terminal_buffer
    sys.stderr = terminal_buffer

    def update_terminal_output():
        """
        Periodically flush the terminal buffer and update the widget.

        The frequency of 10 ms ensures near real-time updates while balancing CPU usage.
        Adjust the interval based on application needs and performance testing.
        """
        terminal_buffer.flush()  # Ensures the widget reflects the latest data
        root.after(10, update_terminal_output)  # Schedule next update at 10ms intervals

    # Start the periodic terminal update
    update_terminal_output()

    def process_file():

        # Print to terminal output widget
        terminal_output.clear()  # Clear previous output
        terminal_output.write("Running with the following parameters:\n")
        # for category, options in results.items():
        #     terminal_output.write(f"{category}: {options}\n")

        ##################################################################################################################################

        def run_pHinder():

            import sys

            try:
                print("Starting pHinder process...")
                sys.stdout.flush()  # Explicit flush to update output

                # Dictionary to collect all widget parameter values
                results = {}

                # Make sure a PDB file path has been selected
                if not file_path_widget.get_file_path():
                    messagebox.showerror("Error", "File path is empty or no chains are selected.")
                    return

                # Get the file path and chains from the FilePathWidget
                file_path = file_path_widget.get_file_path()
                chains = file_path_widget.get_values()
                check_chains = []
                for c in chains:
                    if chains[c]:
                        check_chains.append(c)

                # Make sure a PDB file path and chains have been selected
                if not file_path or not check_chains:
                    messagebox.showerror("Error", "File path is empty or no chains are selected.")
                    return

                results["file_path"] = file_path
                results["chains"] = {chain: var for chain, var in chains.items()}

                # Get the save path
                save_path = file_path_widget.get_save_path()
                if save_path:
                    results["save_path"] = save_path
                else:
                    results["save_path"] = ""

                # Collect values from the Calculation Options Pane
                results["calculation_options"] = {
                    option: var.get() for option, var in calculation_variables.items()
                }

                # Collect values from the Sidechain Classification Options Pane
                results["sidechain_classification_options"] = {
                    option: var.get() for option, var in sidechain_variables.items()
                }

                # Collect values from the Network Options Pane
                results["network_options"] = {
                    option: var.get() for option, var in network_variables.items()
                }

                # Collect values from the Surface Options Pane
                results["surface_options"] = {
                    option: var.get() for option, var in surface_variables.items()
                }

                # Collect values from the Interface Options Pane
                results["interface_options"] = {
                    option: var.get() for option, var in interface_variables.items()
                }

                # Collect values from the Virtual Screening Options Pane
                results["virtual_screening_options"] = {
                    option: var.get() for option, var in virtual_screening_variables.items()
                }

                # Collect values from the Advanced Options Pane
                results["advanced_options"] = {
                    option: var.get() for option, var in advanced_variables.items()
                }

                # Collect amino acid selections
                results["amino_acid_selections"] = {
                    amino_acid: var for amino_acid, var in amino_acid_widget.get_values().items()
                }

                # Import pHinder
                from pHinder.pHinder_7_0 import pHinder

                # Portions of the virtual screening program are parallelized.
                # Set the number of processors to be used in the calculation.
                #############################################################
                import multiprocessing as mp
                n_processes = mp.cpu_count() - 1

                import time
                start_time = time.time()
                from os import sep
                print("Running pHinder for file", file_path.split(sep)[-1])
                print(80*"-")

                # Set the file and save paths and file name
                file_path = results["file_path"]
                pdb_file_path = sep.join(file_path.split(sep)[:-1]) + sep #.join(sep) + sep
                pdb_file_name = file_path.split(sep)[-1]
                save_path = results["save_path"] + sep

                # Set options for the appropriate PDB file handlers
                pdb_format = "pdb"
                if ".cif" in pdb_file_name:
                    pdb_format = "mmCIF"
                # Check if the file is compressed
                zip_status = 0
                if ".gz" in pdb_file_name:
                    zip_status = 1

                # Set PDB chain selection(s)
                chains = results["chains"]
                chain_selection = []
                for c in chains:
                    if chains[c]:
                        chain_selection.append(c)
                chain_selection.sort()
                chains = chain_selection

                # Set group chains options
                # If true, calculations will be done on all protein chains as one
                # If false, calculation will be done for each individual protein chain
                group_chains = results["chains"]["Group Chains"]

                # PHINDER INSTANTIATION AND FUNCTION CALLS.
                ###########################################

                # Create a pHinder instance.
                ############################
                try:
                    pHinderInstance = pHinder()
                except Exception as e:
                    print(f"Error occurred during pHinder execution: {e}")
                    sys.stdout.flush()
                    messagebox.showerror("Error", "File path is empty or no chains or amino acids are selected.")
                    return

                # Set the essential pHinder variables.
                ###################################### 

                # Set gui option
                pHinderInstance.gui = True

                # Set number of available processes
                pHinderInstance.processes = n_processes

                # Set file options
                pHinderInstance.pdbFilePath= pdb_file_path
                pHinderInstance.pdbFileName= pdb_file_name
                pHinderInstance.outPath= save_path + "pHinderResults/"
                pHinderInstance.pdbFormat = pdb_format
                pHinderInstance.zip=zip_status
                pHinderInstance.chains= chains
                pHinderInstance.group_chains=group_chains

                # Set the set of residues to be used in calculations
                aa_set = results["amino_acid_selections"]
                aa_set_keys, aa_set_selected = sorted(aa_set), []
                for aa_set_key in aa_set_keys:
                    if aa_set[aa_set_key]:
                        aa_set_selected.append(aa_set_key)
                aa_set_selected_key = tuple(aa_set_selected)

                aa_sets = {

                            ('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'): "allSet",
                            ('Ala', 'Ile', 'Leu', 'Met', 'Phe', 'Pro', 'Trp', 'Val'): "apolarSet",
                            ('Asn', 'Cys', 'Gln', 'Ser', 'Thr', 'Tyr'):"polarSet",
                            ('Arg', 'Asp', 'Cys', 'Glu', 'His', 'Lys'):"ionizableSet",
                            ('Arg', 'Asp', 'Glu', 'His', 'Lys'):"ionizableSetNoCys",
                            ('Asp', 'Glu'):"acidicSet",
                            ('Arg', 'His', 'Lys'):"basicSet"

                }

                # Standard aa set
                if aa_set_selected_key in aa_sets:
                    residue_set = aa_sets[aa_set_selected_key]
                    pHinderInstance.residueSet = residue_set
                    print("Standard Residue Set:", residue_set)
                # Custom aa set
                else:
                    custom_aa_set = ""
                    for custom_aa in aa_set_selected_key:
                        custom_aa_set += custom_aa.upper() + ","
                    custom_aa_set = "customSet:" + custom_aa_set[:-1]
                    residue_set = custom_aa_set
                    print("Custom Residue Set:", residue_set)
                    pHinderInstance.residueSet = residue_set

                # Set network options
                pHinderInstance.maxNetworkEdgeLength = results["network_options"]['MAX_NETWORK_EDGE_LENGTH']
                pHinderInstance.minNetworkSize = results["network_options"]['MIN_NETWORK_SIZE']
                pHinderInstance.reducedNetworkRepresentation = results["network_options"]['REDUCED_NETWORK_REPRESENTATION']
                pHinderInstance.saveNetworkTriangulation = results["network_options"]['SAVE_NETWORK_TRIANGULATION']

                # Set surface options
                circumSphereRadiusLimit = results["surface_options"]['circumSphereRadiusLimit']
                minArea = results["surface_options"]['minArea']
                pHinderInstance.highResolutionSurface = results["surface_options"]['HIGH_RESOLUTION_SURFACE']
                pHinderInstance.saveSurface = results["surface_options"]['SAVE_SURFACE']
                pHinderInstance.allowSmallSurfaces = results["surface_options"]['ALLOW_SMALL_SURFACES']
                pHinderInstance.saveLigandSurfaces = results["surface_options"]['SAVE_LIGAND_SURFACES']
                pHinderInstance.writeSurfaceCreationAnimation = results["surface_options"]['WRITE_SURFACE_CREATION_ANIMATION']

                # Set sidechain classification options
                pHinderInstance.coreCutoff = results["sidechain_classification_options"]['CORE_CUTOFF']
                pHinderInstance.marginCutoff = results["sidechain_classification_options"]['MARGIN_CUTOFF']
                pHinderInstance.marginCutoffCoreNetwork = results["sidechain_classification_options"]['MARGIN_CUTOFF_CORE_NETWORK']

                # Set interface options
                pHinderInstance.interface_distance_filter = results["interface_options"]["INTERFACE_DISTANCE_FILTER"]

                # Set virtual screening options
                maxVoidNetworkEdgeLength = results["virtual_screening_options"]['MAX_VOID_NETWORK_EDGE_LENGTH']
                minVoidNetworkSize = results["virtual_screening_options"]['MIN_VOID_NETWORK_SIZE']
                pHinderInstance.virtualClashCutoff = results["virtual_screening_options"]['VIRTUAL_CLASH_CUTOFF']
                pHinderInstance.inIterations = results["virtual_screening_options"]['IN_ITERATIONS']
                pHinderInstance.inIterationsStepSize = results["virtual_screening_options"]['IN_ITERATIONS_STEP_SIZE']
                pHinderInstance.outIterations = results["virtual_screening_options"]['OUT_ITERATIONS']
                pHinderInstance.outIterationsStepSize = results["virtual_screening_options"]['OUT_ITERATIONS_STEP_SIZE']

                # Set advanced options
                pHinderInstance.allowCysCoreSeeding = results["advanced_options"]['ALLOW_CYS_CORE_SEEDING']
                pHinderInstance.includeHydrogens = results["advanced_options"]['INCLUDE_HYDROGENS']
                pHinderInstance.includeWater = results["advanced_options"]['INCLUDE_WATER']
                pHinderInstance.includeIons = results["advanced_options"]['INCLUDE_IONS']

                # Set the calculation types
                topologyCalculation = results["calculation_options"]['topologyCalculation']
                surfaceCalculation = results["calculation_options"]['surfaceCalculation']
                sidechainClassification = results["calculation_options"]['sidechainClassification']
                interfaceClassification = results["calculation_options"]['interfaceClassification']
                virtualScreenSurfacesCalculation = results["calculation_options"]['virtualScreenSurfacesCalculation']

                # SET PYTHON RECURSION LIMIT.
                #############################
                # Python recursion limit (Default is 1000).
                # This may be needed to accommodate goFo recursion in larger calculations.
                ##########################################################################
                import sys
                print("Default recursion limit is:", sys.getrecursionlimit())
                python_recursion_limit = int(results["advanced_options"]['PYTHON_RECURSION_LIMIT'])
                sys.setrecursionlimit(python_recursion_limit)
                print("Recursion limit increased to:", sys.getrecursionlimit())

                # Essential/base pHinder function calls.
                ########################################
                pHinderInstance.setQuerySet()
                pHinderInstance.openPDBs(pdb_file_path, pdb_file_name, zip_status=pHinderInstance.zip)
                pHinderInstance.hetLigand4D()
                # pHinderInstance.hetWater4D()
                pHinderInstance.hydrogens()
                pHinderInstance.makeAtomCollections()
                pHinderInstance.makeVertices4D()

                write_phinder_log(pHinderInstance, pHinderInstance.outPath + "pHinder_parameters.log")

                # The list of pHinder calcuation options
                if topologyCalculation:

                    # Residue Network Topologies.
                    #############################

                    # Triangulate the user defined residue set.
                    ###########################################
                    pHinderInstance.selectTscTriangulationAtoms()
                    pHinderInstance.triangulateTscAtoms() # Done with chain adaptation...
                    pHinderInstance.writeTriangulation() # Done with chain adaptation...

                    # Calculate network topologies based on the user defined residue set.
                    #####################################################################
                    pHinderInstance.pruneTriangulation() # Done with chain adaptation...
                    pHinderInstance.minimizePrunedTriangulation() # Done with chain adaptation...

                    # Analyze network topologies based on the user defined residue set.
                    ###################################################################
                    pHinderInstance.identifyTightBonds() # Done with chain adaptation...
                    pHinderInstance.calculateNetworkParity() # Done with chain adaptation...

                if surfaceCalculation:

                    # Calculate the pHinder protein surface.
                    ########################################
                    pHinderInstance.surface(circumSphereRadiusLimit=circumSphereRadiusLimit, minArea=minArea) # Done with chain adaptation...
                    pHinderInstance.writeSurface() # Done with chain adaptation...
                    pHinderInstance.surfaceLigands() # Done with chain adaptation...
                    pHinderInstance.writeLigandSurfaces() # Done with chain adaptation...

                # Functions that require a protein surface.
                ###########################################

                if sidechainClassification:

                    if not topologyCalculation:

                        # Triangulate the user defined residue set.
                        ###########################################
                        pHinderInstance.selectTscTriangulationAtoms()
                        pHinderInstance.triangulateTscAtoms() # Done with chain adaptation...
                        pHinderInstance.writeTriangulation() # Done with chain adaptation...

                        # Calculate network topologies based on the user defined residue set.
                        #####################################################################
                        pHinderInstance.pruneTriangulation() # Done with chain adaptation...
                        pHinderInstance.minimizePrunedTriangulation() # Done with chain adaptation...

                        # Analyze network topologies based on the user defined residue set.
                        ###################################################################
                        pHinderInstance.identifyTightBonds() # Done with chain adaptation...
                        pHinderInstance.calculateNetworkParity() # Done with chain adaptation...


                    if not surfaceCalculation:

                        # Calculate the pHinder protein surface.
                        ########################################
                        pHinderInstance.surface(circumSphereRadiusLimit=circumSphereRadiusLimit, minArea=minArea) # Done with chain adaptation...
                        pHinderInstance.writeSurface() # Done with chain adaptation...
                        pHinderInstance.surfaceLigands() # Done with chain adaptation...
                        pHinderInstance.writeLigandSurfaces() # Done with chain adaptation...

                    # Classify sidechains.
                    ######################
                    # pHinderInstance.triangulateAllTscAtomsForClassfication()
                    pHinderInstance.selectTscClassificationAtoms() 
                    pHinderInstance.classifySidechainLocation() # Done with chain adaptation...
                    pHinderInstance.identifyMissingTscAtoms() # Done with chain adaptation...
                    pHinderInstance.writeSidechainClassificationResults() # Done with chain adaptation...

                if interfaceClassification:

                    # Classify sidechains at the interface of one or more protein chains.
                    #####################################################################

                    if not topologyCalculation:

                        # Triangulate the user defined residue set.
                        ###########################################
                        pHinderInstance.selectTscTriangulationAtoms()
                        pHinderInstance.triangulateTscAtoms() # Done with chain adaptation...
                        pHinderInstance.writeTriangulation() # Done with chain adaptation...

                        # Calculate network topologies based on the user defined residue set.
                        #####################################################################
                        pHinderInstance.pruneTriangulation() # Done with chain adaptation...
                        pHinderInstance.minimizePrunedTriangulation() # Done with chain adaptation...

                        # Analyze network topologies based on the user defined residue set.
                        ###################################################################
                        pHinderInstance.identifyTightBonds() # Done with chain adaptation...
                        pHinderInstance.calculateNetworkParity() # Done with chain adaptation...

                    if not surfaceCalculation:

                        # Calculate the pHinder protein surface.
                        ########################################
                        pHinderInstance.surface(circumSphereRadiusLimit=circumSphereRadiusLimit, minArea=minArea) # Done with chain adaptation...
                        pHinderInstance.writeSurface() # Done with chain adaptation...
                        pHinderInstance.surfaceLigands() # Done with chain adaptation...
                        pHinderInstance.writeLigandSurfaces() # Done with chain adaptation...

                    if not sidechainClassification:

                        # Classify sidechains.
                        ######################
                        # pHinderInstance.triangulateAllTscAtomsForClassfication()
                        pHinderInstance.selectTscClassificationAtoms() 
                        pHinderInstance.classifySidechainLocation() # Done with chain adaptation...
                        pHinderInstance.identifyMissingTscAtoms() # Done with chain adaptation...
                        pHinderInstance.writeSidechainClassificationResults() # Done with chain adaptation...

                    # Functions that require a protein surface.
                    ###########################################
                    pHinderInstance.classifyInterfaceSidechains()

                if virtualScreenSurfacesCalculation:

                    if not topologyCalculation:

                        # Triangulate the user defined residue set.
                        ###########################################
                        pHinderInstance.selectTscTriangulationAtoms()
                        pHinderInstance.triangulateTscAtoms() # Done with chain adaptation...
                        pHinderInstance.writeTriangulation() # Done with chain adaptation...

                        # Calculate network topologies based on the user defined residue set.
                        #####################################################################
                        pHinderInstance.pruneTriangulation() # Done with chain adaptation...
                        pHinderInstance.minimizePrunedTriangulation() # Done with chain adaptation...

                        # Analyze network topologies based on the user defined residue set.
                        ###################################################################
                        pHinderInstance.identifyTightBonds() # Done with chain adaptation...
                        pHinderInstance.calculateNetworkParity() # Done with chain adaptation...

                    if not surfaceCalculation:

                        # Calculate the pHinder protein surface.
                        ########################################
                        pHinderInstance.surface(circumSphereRadiusLimit=circumSphereRadiusLimit, minArea=minArea) # Done with chain adaptation...
                        pHinderInstance.writeSurface() # Done with chain adaptation...
                        pHinderInstance.surfaceLigands() # Done with chain adaptation...
                        pHinderInstance.writeLigandSurfaces() # Done with chain adaptation...

                    # Virtual screening: void volumes.
                    ##################################
                    pHinderInstance.makeSamplingGridUsingProteinSurface()
                    pHinderInstance.filterSamplingPointsUsingClashes()
                    pHinderInstance.triangulateRemainingGridPoints()
                    pHinderInstance.identifyAndParseIndividualSamplingVoids(maxVoidNetworkEdgeLength=maxVoidNetworkEdgeLength, minVoidNetworkEdgeLength=0.0, minVoidNetworkSize=minVoidNetworkSize, psa=1) # XXXX new
                    pHinderInstance.calculateSamplingVoidSurfaces(extend_sampling=True)

                duration = time.time() - start_time 
                print("Runtime (s):", duration)

                print("pHinder process completed successfully.")
                sys.stdout.flush()
            except Exception as e:
                print(f"Error occurred during pHinder execution: {e}")
                sys.stdout.flush()
                messagebox.showerror("Error", "File path is empty or no chains or amino acids are selected.")
                return

        # Run the pHinder process in a separate thread
        thread = threading.Thread(target=run_pHinder)
        thread.start()

        ##################################################################################################################################

    ttk.Button(left_frame, text="Run pHinder", command=process_file).pack(pady=10, anchor="w")

    root.mainloop()


if __name__ == "__main__":

    main()


