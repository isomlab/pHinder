# Set file paths.
#################
# where_to_save_output = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/SuperDarks/"
# where_the_pdb_file_is = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/SuperDarks/pdbFiles/"
# pdbCode = "AF-Q5FWE3-F1-model_v3--7vv5.R-aligned--199-0.69506"
# pdbChains = []

where_to_save_output = "/Users/danisom/Desktop/scott/"
where_the_pdb_file_is = "/Users/danisom/Desktop/scott/"
pdbCode = "8HQE"
pdbChains = []

# where_to_save_output = "/projectnb/isomlab/dan/in-and-out/input/alphaFold_v4/alphaFoldQueryStructures-Pegasus/rhod/"
# where_the_pdb_file_is = "/projectnb/isomlab/dan/in-and-out/input/alphaFold_v4/alphaFoldQueryStructures-Pegasus/rhod/"
# pdbCode = "RHO-1F88-1-348-A"
# pdbChains = []

# where_to_save_output = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/testing/"
# where_the_pdb_file_is = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/testing/"
# pdbCode = "GPR68-AF-Q15743-F1-model_v3--7wvu.R-aligned--235-0.8114-trimmed"
# pdbChains = [] 

# where_to_save_output = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/testing/"
# where_the_pdb_file_is = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/testing/"
# pdbCode = "4o9u-c-selection"
# pdbChains = ["C", "D"]

# where_to_save_output = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/cif_format/cif_files/proteome-tax_id-1630578-0_v3/"
# where_the_pdb_file_is = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/cif_format/cif_files/proteome-tax_id-1630578-0_v3/"
# pdbCode = "6ks0-aligned"
# pdbChains = []

# where_to_save_output = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/asics_macleod/pdb_files/"
# where_the_pdb_file_is = "/Users/danisom/BitTorrentSync/CalcSync-PhinderScratch/asics_macleod/pdb_files/"
# pdbCode = "2qts.A.B.C"
# pdbChains = ["A", "B", "C"]

# PDB format
# All pHinder output is in PDB format
# input options are "pdb" and "mmCIF"
#####################################
pdb_format = "mmCIF"

# Set zip status of the input PDB file.
#######################################
ZIP = 0

# If specific PDB chains are not specified...
# Should the chains be grouped or processed individually in the calculation?
############################################################################
group_chains = 0

# Set pHinder version.
######################
from pHinder.pHinder_7_0 import pHinder

# Select pHinder capabilities.
##############################
topologyCalculation = 1
surfaceCalculation = 0
sidechainClassification = 0
interfaceClassification = 0
virtualScreenSurfacesCalculation = 0

# Set surface carving and refinement parameters.
################################################
# This paramater controls the aggressiveness of surface carving. 
# A large value will lead to little or no carving.
# A small value will lead to overcarving and a "swiss cheese" like surface result.
# For proteins, a value of 6.5 consistently produces the best results.
# Tip: I never (or hardly ever) change this value for proteins.
# Tip: I do lower this value when carving surfaces for ligands or sampling voids.
circumSphereRadiusLimit = 6.5 
# Minimum facet area within the pHinder calculated surface.
# For proteins, a value of 10-20 Anstroms is ample for accurate side chain classification.
# The smaller the facet area, the more facets, the longer the surface calculation will take.
# For example, a value of 2 produces a very attractive surface, but takes longer to calculate.
# Tip: Use 10 unless you need a more "attractive" surface for a figure.
minArea = 10

# Select the amino acid residue set to be used in the calculation.
# Options:
# allSet: all residues
# ionizableSet: ionizable residues
# ionizableSetNoCys: ionizable residue, but excluding Cys
# acidicSet: acidic residues
# basicSet: basic residues 
# polarSet: polar residues
# apolarSet: apolar residues
# customSet: example for ASP and GLU only: "customSet:ASP,GLU"
####################################################################
RESIDUE_SET = "ionizableSetNoCys"    

# CORE PHINDER PARAMETERS.
##########################
# Tip: I hardly ever change the core parameters listed below.

# Set pHinder function parameters.
##################################
CORE_CUTOFF = -3.0 
MARGIN_CUTOFF = -2.0
MARGIN_CUTOFF_CORE_NETWORK = -2.0

# Set network parameters.
#########################
MAX_NETWORK_EDGE_LENGTH = 10.0 # Default = 10.0, but maybe it should be 12 moving forward?
MIN_NETWORK_SIZE = 1
REDUCED_NETWORK_REPRESENTATION = 1
SAVE_NETWORK_TRIANGULATION = 1

# Set pHinder-specific network paramater(s).
############################################
ALLOW_CYS_CORE_SEEDING = 1

# Special atom selections.
##########################
INCLUDE_HYDROGENS = 0
INCLUDE_WATER = 0
INCLUDE_IONS = 0

# Surface parameters.
#####################
HIGH_RESOLUTION_SURFACE = 1
SAVE_SURFACE = 1
ALLOW_SMALL_SURFACES = 1
SAVE_LIGAND_SURFACES = 0
WRITE_SURFACE_CREATION_ANIMATION = 0

# Python recursion limit (Default is 1000).
# This may be needed to accommodate goFo recursion in larger calculations.
##########################################################################
PYTHON_RECURSION_LIMIT = 10000

# VIRTUAL SCREENING PARAMETERS.
###############################

# Portions of the virtual screening program are parallelized.
# Set the number of processors to be used in the calculation.
#############################################################
import multiprocessing as mp
PROCESSES = mp.cpu_count() - 1

# Set sampling parameters for virtual screening.
################################################
MAX_VOID_NETWORK_EDGE_LENGTH = 2.0 # Default 2.0, 1.0 for olfactory; should to be <= IN_ITERATIONS_STEP_SIZE and OUT_ITERATIONS_STEP_SIZE
MIN_VOID_NETWORK_SIZE = 10 # Default 1

# Only relevant if using a cubic grid to generating sampling points, which is not my preferred approach.
########################################################################################################
GRID_INCREMENT = 3.0 # Default 2.0

# This must be at least 2.5 Angstroms (Default value is 2.5).
# If smaller, the surface voids that clash with the protein backbone,
# and coalesce into a void surface that snakes through the interior of the protein.
###################################################################################
# This parameter is critical for allowing voids. It has repeatedly proved to be a VERY sensitive parameter. Must be between 2 and 3.
####################################################################################################################################
VIRTUAL_CLASH_CUTOFF = 2.5 # Default 2.5 

# Below (in) and above (out) iteration parameters.
# These are used if generating sampling points using the proteins convex hull or protein surface.
#################################################################################################
IN_ITERATIONS = 1 # Default 1 ... 2022.04.14... was 7 for latest GPCR VS
IN_ITERATIONS_STEP_SIZE = 2.0 # Default 1
OUT_ITERATIONS = 1 #10 # Default 6 Pointless is proteins like Ras, 4LUC.A, where these points are remove by inHull...
OUT_ITERATIONS_STEP_SIZE = 2.0 # Default 1; 3.0 is typically too much for chemiforms

if __name__ == "__main__": 

	import time
	start_time = time.time()
	print("Running pHinder for PDB code", pdbCode)
	print(80*"-")

	# SET PYTHON RECURSION LIMIT.
	#############################
	# Python recursion limit (Default is 1000).
	# This may be needed to accommodate goFo recursion in larger calculations.
	##########################################################################
	import sys
	print("Default recursion limit is:", sys.getrecursionlimit())
	sys.setrecursionlimit(PYTHON_RECURSION_LIMIT)
	print("Recursion limit increased to:", sys.getrecursionlimit())

	# DIRECTORY STRUCTURE PARAMETERS.
	#################################
	# Set the input and output paths to the calculation data.
	#########################################################
	OUTPATH = where_to_save_output + "pHinderResults/"
	PATH = where_the_pdb_file_is


	# Set the list of PDB file paths (absolute paths).
	##################################################
	pdbFilePath, pdbFileName = PATH, ""
	if pdb_format == "mmCIF":
		if ZIP:
			pdbFileName = pdbCode + ".cif.gz"
		else:
			pdbFileName = pdbCode + ".cif"
	else:
		if ZIP:
			pdbFileName = pdbCode + ".pdb.gz"
		else:
			pdbFileName = pdbCode + ".pdb"

	# PHINDER INSTANTIATION AND FUNCTION CALLS.
	###########################################

	# Create a pHinder instance.
	############################
	pHinderInstance = pHinder()

	# Set the essential pHinder variables.
	###################################### 
	pHinderInstance.pdbFormat = pdb_format
	pHinderInstance.pdbFilePath=pdbFilePath
	pHinderInstance.pdbFileName=pdbFileName
	pHinderInstance.outPath=OUTPATH 
	pHinderInstance.chains=pdbChains 
	pHinderInstance.group_chains=group_chains
	pHinderInstance.maxNetworkEdgeLength=MAX_NETWORK_EDGE_LENGTH
	pHinderInstance.minNetworkSize=MIN_NETWORK_SIZE
	pHinderInstance.includeWater=INCLUDE_WATER
	pHinderInstance.includeIons=INCLUDE_IONS
	pHinderInstance.zip=ZIP
	pHinderInstance.reducedNetworkRepresentation=REDUCED_NETWORK_REPRESENTATION
	pHinderInstance.saveNetworkTriangulation=SAVE_NETWORK_TRIANGULATION
	pHinderInstance.highResolutionSurface=HIGH_RESOLUTION_SURFACE
	pHinderInstance.saveSurface=SAVE_SURFACE
	pHinderInstance.allowSmallSurfaces=ALLOW_SMALL_SURFACES
	pHinderInstance.saveLigandSurfaces=SAVE_LIGAND_SURFACES
	pHinderInstance.residueSet=RESIDUE_SET
	pHinderInstance.allowCysCoreSeeding=ALLOW_CYS_CORE_SEEDING
	pHinderInstance.writeSurfaceCreationAnimation=WRITE_SURFACE_CREATION_ANIMATION
	pHinderInstance.coreCutoff=CORE_CUTOFF
	pHinderInstance.marginCutoff=MARGIN_CUTOFF
	pHinderInstance.marginCutoffCoreNetwork=MARGIN_CUTOFF_CORE_NETWORK
	pHinderInstance.includeHydrogens=INCLUDE_HYDROGENS
	pHinderInstance.processes=PROCESSES
	pHinderInstance.gridIncrement=GRID_INCREMENT
	pHinderInstance.maxVoidNetworkEdgeLength=MAX_VOID_NETWORK_EDGE_LENGTH
	pHinderInstance.minVoidNetworkSize=MIN_VOID_NETWORK_SIZE
	pHinderInstance.virtualClashCutoff=VIRTUAL_CLASH_CUTOFF
	pHinderInstance.inIterations=IN_ITERATIONS
	pHinderInstance.inIterationsStepSize=IN_ITERATIONS_STEP_SIZE
	pHinderInstance.outIterations=OUT_ITERATIONS
	pHinderInstance.outIterationsStepSize=OUT_ITERATIONS_STEP_SIZE

	# Essential/base pHinder function calls.
	########################################
	pHinderInstance.setQuerySet()
	pHinderInstance.openPDBs(pdbFilePath, pdbFileName, zip_status=pHinderInstance.zip)
	pHinderInstance.hetLigand4D()
	pHinderInstance.hydrogens()
	pHinderInstance.makeAtomCollections()
	pHinderInstance.makeVertices4D()

	########################################################
	#
	# Optional pHinder function calls listed by group below.
	#
	########################################################

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

		# Protein Surface.
		##################

		# Calculate the pHinder surface.
		################################
		pHinderInstance.surface(circumSphereRadiusLimit=circumSphereRadiusLimit, minArea=minArea) # Done with chain adaptation...
		pHinderInstance.writeSurface() # Done with chain adaptation...
		pHinderInstance.surfaceLigands() # Done with chain adaptation...
		pHinderInstance.writeLigandSurfaces() # Done with chain adaptation...

	# Functions that require a protein surface.
	###########################################

	if sidechainClassification and surfaceCalculation:

		# Functions that require a protein surface.
		###########################################

		# Classify sidechains.
		######################
		# pHinderInstance.triangulateAllTscAtomsForClassfication()
		pHinderInstance.selectTscClassificationAtoms() 
		pHinderInstance.classifySidechainLocation() # Done with chain adaptation...
		pHinderInstance.identifyMissingTscAtoms() # Done with chain adaptation...
		pHinderInstance.writeSidechainClassificationResults() # Done with chain adaptation...

	if interfaceClassification and surfaceCalculation:

		# Classify sidechains at the interface of one or more protein chains.
		#####################################################################

		# Functions that require a protein surface.
		###########################################
		pHinderInstance.classifyInterfaceSidechains()


	if virtualScreenSurfacesCalculation and surfaceCalculation:

		# # For class...
		# pHinderInstance.makeSamplingGridCubic()
		# pHinderInstance.filterSamplingPointsUsingConvexHull3D()
		# pHinderInstance.filterSamplingPointsUsingProximity(proximityLimit=4.0)

		# Virtual screening: void volumes.
		##################################
		pHinderInstance.virtualScreen()
		pHinderInstance.inIterations=IN_ITERATIONS #+10 #4
		pHinderInstance.inIterationsStepSize=IN_ITERATIONS_STEP_SIZE #1 # was last 0.5
		pHinderInstance.outIterations=OUT_ITERATIONS #+6 #4
		pHinderInstance.outIterationsStepSize=OUT_ITERATIONS_STEP_SIZE #1 # was last 0.5
		pHinderInstance.makeSamplingGridUsingProteinSurface()
		pHinderInstance.filterSamplingPointsUsingClashes()
		#pHinderInstance.filterSamplingPointsUsingProximity(proximityLimit=0.9) # Work in progress; maybe not necessary.
		pHinderInstance.triangulateRemainingGridPoints()
		# pHinderInstance.samplingReduction(maxVoidNetworkEdgeLength=MAX_VOID_NETWORK_EDGE_LENGTH, minVoidNetworkEdgeLength=1.5)
		#pHinderInstance.samplingReduction(minVoidNetworkEdgeLength=0.9)
		pHinderInstance.identifyAndParseIndividualSamplingVoids(maxVoidNetworkEdgeLength=MAX_VOID_NETWORK_EDGE_LENGTH, minVoidNetworkEdgeLength=0.0, minVoidNetworkSize=MIN_VOID_NETWORK_SIZE, psa=1) # XXXX new
		pHinderInstance.calculateSamplingVoidSurfaces(extend_sampling=True)
		# pHinderInstance.reduceChemiforms(psa=1)

	duration = time.time() - start_time 
	print("Runtime (s):", duration)



