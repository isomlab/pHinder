global zero
zero = 1e-10

# Import dependencies here for enhanced performance.
####################################################
from pHinder._vendor.pdbFile import PseudoAtom
from pHinder._vendor.compGeometry import *
from math import sqrt
from copy import deepcopy
# from tempfile import TemporaryFile
# from xlwt import Workbook
from time import perf_counter
from os import mkdir, sep
from os.path import exists
from pHinder.geometry.convexHull4D import convexHull4D
import gzip
from pHinder.geometry.goFo import *

# Residue Sets (Gly is skipped).
#################################
ionizableSet = ["ASP","GLU","HIS","CYS","LYS","ARG"]
polarSet = ["ASN","GLN","SER","THR","TYR","TRP"]
apolarSet = ["ALA","ILE","LEU","MET","PHE","PRO","VAL"]
allSet = ionizableSet + polarSet + apolarSet

def s3StringOutput(v1,v2,v3):

	format = "%15.6f%15.6f%15.6f "
	string = format % (v1.x,v1.y,v1.z) + \
			format % (v2.x,v2.y,v2.z) + \
			format % (v3.x,v3.y,v3.z) + "\n"
	return string

# Non-essential function.
#########################
def writeEvaluationSurface(filePath, identifier, surfaceFacets):
	
	surfaceFacetString = ""
	for fKey in surfaceFacets:
		f = surfaceFacets[fKey]
		surfaceFacetString += s3StringOutput(f.v1, f.v2, f.v3)

		outFile = filePath + "evaluationSurface." + identifier + ".txt.gz"
		outputSurf = gzip.open(outFile,"wt")
		outputSurf.write(surfaceFacetString)
		outputSurf.close()

def writeSubTesselation(filePath, identifier, currentSurface, hull4D):
	
	outFile = filePath + "." + identifier + ".txt.gz"
	outputTesselation = gzip.open(outFile,"wt")

	surfaceFacetsString = ""
	for s3Key in currentSurface:
		
		triangleInstance = currentSurface[s3Key]
		t = triangleInstance

		format = "%15.6f%15.6f%15.6f "
		surfaceFacetsString += format % (t.v1.x, t.v1.y, t.v1.z) + \
								format % (t.v2.x, t.v2.y, t.v2.z) + \
								format % (t.v3.x, t.v3.y, t.v3.z) + "\n"

	outputTesselation.write(surfaceFacetsString)
	outputTesselation.close()

def getInitialSurface(subTess3D, s3CompleteLocal):


	# STEP II
	# Identify exposed s3 and s4 objects of the triangulation.
	# Also scan for, and remove, isolated s4 pairs. 
	###########################################################
	s3Exposed, s4Exposed, s4Isolated, s3RemoveDict = {}, {}, {}, {}
	for s3Key in s3CompleteLocal:

		# A) ISOLATED TETRAHEDRA.
		#
		# 2012.09.24 
		# Isolate surface tetrahedra:
		# Identify isolated s3 objects incident on exactly two s4 objects.
		# These objects comprise an isolated s4 pair that should be deleted
		# to avoid infinite while looping using the toTrim function output.
		#
		# This situation typically occurs in extended regions of structure, 
		# like N- or C-terminal peptides, where the surface algorithm 
		# rapidly erodes the thin surface contribution from the peptide
		# stretch.  In the end, such peptide are typically eroded to the 
		# extent that they are completely erased and make no contribution
		# to the pHinder surface. 
		#
		# The calculation is only relevant to s3 objects incident on two s4 objects.
		############################################################################
		if s3Key and len(s3CompleteLocal[s3Key]) == 2:
			s41, s42, s3List = s3CompleteLocal[s3Key][0], s3CompleteLocal[s3Key][1], []
			for s3Key2 in s41.subKeys:
				if s3Key2 != s3Key:
					s3List.append(s3Key2)
			for s3Key2 in s42.subKeys:
				if s3Key2 != s3Key: 
					s3List.append(s3Key2)

			# Isolated s4 pairs each have three unique s3 objects,
			# to give a total of six unique s3 objects per s4 pair.
			#######################################################
			t = 0
			for s3Test in s3List:
				if len(s3CompleteLocal[s3Test]) == 1:
					t += 1

			# If true, then this s4 pair is isolated.
			#########################################
			if t == 6:
				
				s4Isolated.update({s41:s41.subKeys})
				s4Isolated.update({s42:s42.subKeys})

				# Remove all of the s3 objects from subTess3D that are incident on the s4 pair.
				###############################################################################
				s3RemoveList = [s3Key,] + s3List
				for s3RemoveKey in s3RemoveList:
					if s3RemoveKey in subTess3D:
						del subTess3D[s3RemoveKey]
			
				# The exposed s3s of the isolated s4 pair need to be saved for removal from s3Exposed.
				######################################################################################
				for s3RemoveKey in s3RemoveList:
					s3RemoveDict.update({s3RemoveKey:None})
	
		# B) EXPOSED TETRAHEDRA
		# Identify s3Keys incident on only one s4.  These s3 objects are exposed.
		#########################################################################
		if len(s3CompleteLocal[s3Key]) == 1:
			s3Exposed.update({s3Key:s3CompleteLocal[s3Key]}) 
			s4 = s3CompleteLocal[s3Key][0]
			if s4 not in s4Exposed:
				s4Exposed.update({s4:[s3Key,]})
			elif s3Key not in s4Exposed[s4]:
				s4Exposed[s4].append(s3Key)

	# The exposed s3s of the isolated s4 pair need to be saved and removed from s3Exposed.
	#####################################################################################
	for s3Key in s3RemoveDict:
		if s3Key in s3Exposed:
			del s3Exposed[s3Key]


	# Remove any isolated s4 objects from the emerging pHinder surface.
	###################################################################
	for s4 in s4Isolated:
		if s4 in s4Exposed:
			del s4Exposed[s4]

	# This code links exposed s4 object to each other through their incidenct s3 facets.
	# The subset of s4s in s4Coincident are subsequently used in the function surfSizeTrim
	# to carve out the pHinder surface.
	#######################################################################################
	s4Coincident = {}
	for s4 in s4Exposed:
		for s3Key in s4.subKeys:
			if s3Key not in s4Exposed[s4]:
				if s4 not in s4Coincident:
					s4Coincident.update({s4:[s3Key]})
				elif s3Key not in s4Coincident[s4]:
					s4Coincident[s4].append(s3Key)

	return (s3Exposed, s4Exposed, s4Isolated, s4Coincident)

def identifySurface(subTess3D,tess3D,s3Complete={},checkSurfaceOn=0,noDelete=0):

	# STEP I 
	# Establish an empty dictionary for all s3 objects of subTess3D.
	################################################################

	# Go over full list: first time through.
	# This collects ALL of the s3 edges of the triangulation.
	# In subsequent runs, I just want to collect s3 edges within ?? what object
	###########################################################################

	# Smart search: after first time through.
	#########################################
	s3CompleteLocal = {}
	if not s3Complete:

		# Go over full list: first time through.
		# This collects ALL of the s3 edges of the triangulation.
		# In subsequent runs, I just want to collect s3 edges within ?? what object
		###########################################################################
		for s3Key in subTess3D:
			for s4 in subTess3D[s3Key]:
				for s3id in s4.subKeys:
					if s3id not in s3CompleteLocal:
						s3CompleteLocal.update({s3id:[s4]})
					elif s4 not in s3CompleteLocal[s3id]:
						s3CompleteLocal[s3id].append(s4)
	else:
		s3CompleteLocal = s3Complete

	# STEP II
	# Identify exposed s3 and s4 objects of the triangulation. Also scan for, and remove, isolated s4 pairs.
	########################################################################################################
	results = getInitialSurface(subTess3D, s3CompleteLocal)
	s3Exposed, s4Exposed, s4Isolated, s4Coincident = results[0], results[1], results[2], results[3]
	if checkSurfaceOn:
		
		# STEP III
		# Identify topology errors on the emerging pHinder surface.
		# Topology check requires a fully populated s3exposed facet subset calculated in STEP II.
		#########################################################################################

		# The function checkSurface() has the potential to alter s3Exposed and s4Exposed.
		# Iteratively call checkSurface() and identifySurface() until all topology errors
		# in the emerging pHinder surface have been removed.
		###################################################################################
		result = checkSurface(tess3D, subTess3D, s3CompleteLocal, s3Exposed, s4Exposed)
		topologyErrors, s3CompleteLocal, s3Exposed, s4Exposed = result[0], result[1], result[2], result[3]
		while topologyErrors:

			# Update the subset of exposed s3 and s4 objects.
			#################################################
			result = checkSurface(tess3D, subTess3D, s3CompleteLocal, s3Exposed, s4Exposed)
			topologyErrors, s3CompleteLocal, s3Exposed, s4Exposed = result[0], result[1], result[2], result[3]

		# Update s4Coincident after topology errors have been fixed.
		############################################################
		s4Coincident = {}
		for s4 in s4Exposed:
			for s3Key in s4.subKeys:
				if s3Key not in s4Exposed[s4]:
					if s4 not in s4Coincident:
						s4Coincident.update({s4:[s3Key]})
					elif s3Key not in s4Coincident[s4]:
						s4Coincident[s4].append(s3Key)
	
	# Return the surface.
	# There are no s4 objects to delete because the topology of the emerging pHinder surface is correct. 
	####################################################################################################
	return (s3Exposed, s4Exposed, s4Coincident, s3CompleteLocal)

def orientSurfaceTopology(s3Dict, s3Exposed, s4Exposed):#, s4Coincident):

	# This is where the magic happens.  This function uses the fact
	# that an exposed facet (3-simplex) is "above" the remaining
	# "reference" vertex that comprises the surface-exposed tetrahedral cell.
	# This topological relationship is used to orient the circuit of the
	# 3 facet vertices such that vertices "below" the plane of the facet
	# are "inside" the surface, and those "above" the plane of the facet
	# are "outside" the surface.
	#########################################################################


	# This may be more consistent using s4 centroid, but still have rotation error...
	orientedSurface = {}
	for s4 in s4Exposed:
		s4_centroid = centroid4D(s4.s4, PseudoAtom())
		for s3Key in s4Exposed[s4]:
			s3vList = s3Dict[s3Key].s3
			v1, v2, v3 = s3vList[0], s3vList[1], s3vList[2]
			orientedFacet = Triangle(v1, v2, v3, VertexR=s4_centroid)
			oV1 = orientedFacet.v1
			oV2 = orientedFacet.v2
			oV3 = orientedFacet.v3
			orientedFacet = Triangle(oV2, oV1, oV3, skip_orientation=1)
			orientedSurface.update({s3Key:orientedFacet})
	return orientedSurface 

	# orientedSurface = {}
	# for s4 in s4Exposed:
	# 	for s3Key in s4Exposed[s4]:

	# 		# May skip s3Keys already encountered to save time? XXXX

	# 		# 4-Simplex object.
	# 		###################
	# 		s4 = s3Exposed[s3Key][0]

	# 		s3vList = s3Dict[s3Key].s3
	# 		v1, v2, v3 = s3vList[0], s3vList[1], s3vList[2]
	# 		refV = None
	# 		for v in s4.s4:
	# 			if v not in [v1,v2,v3]:
	# 				refV = v
	# 				break

	# 		# XXXX searching for consistency, think this is the problem, may be the answer, no
	# 		# t = Tetrahedron(v1, v2, v3, refV)
	# 		# for e in t.edges:
	# 		#     if s3Key == t.edges[e].rf.id:
	# 		#         orientedSurface[s3Key] = t.edges[e].rf
	# 		#     elif s3Key == t.edges[e].lf.id:
	# 		#         orientedSurface[s3Key] = t.edges[e].lf

	# 		orientedFacet = Triangle(v1, v2, v3, VertexR=refV)
	# 		# The proper orientation is actually reversed for
	# 		# the surface.  I realized this property only by
	# 		# emperical inspection.
	# 		#################################################
	# 		oV1 = orientedFacet.v1
	# 		oV2 = orientedFacet.v2
	# 		oV3 = orientedFacet.v3
	# 		orientedFacet = Triangle(oV2, oV1, oV3, skip_orientation=1)
	# 		orientedSurface.update({s3Key:orientedFacet})

	# return orientedSurface 	

def checkSurface(tess3D, subTess3D, s3Complete, s3Exposed, s4Exposed):

	# Loop over all facets of the surface.
	######################################
	# checkedSurface, surfaceS3, deleteS3, s2Exposed = {}, {}, {}, {}
	checkedSurface, surfaceS3, deleteS3 = {}, {}, {} # XXXX 3/7
	for s3Key1 in s3Exposed:
		
		# !!!! Should not be necessary...likely remove.
		# Surface facets are incident on only one simplex-4.
		if len(s3Exposed[s3Key1]) != 1:
			continue

		# First facet information.
		##########################
		s31  = tess3D.s3Dict[s3Key1].getFacet()
		s211, s212, s213 = s31.e1, s31.e2, s31.e3
		# Get edges from 3-simplex 1.
		#############################
		eids1 = {s211.id:s211, s211.flipid:s211, 
				s212.id:s212, s212.flipid:s212, 
				s213.id:s213, s213.flipid:s213}

		# The exposed s3 as a Simplex3 object.
		######################################
		s31s3 = tess3D.s3Dict[s3Key1] # first facet

		# Within the Simplex3 object, access the s2 ids (subKeys).
		##########################################################
		topologyCheck = {}
		for s2id in s31s3.subKeys:

			# Collect all of the s3 objects that contain the s2.
			####################################################
			s2 = tess3D.s2Dict[s2id]
			for s2s3Key in s2.s3s:
				# Save s3Keys that belong to exposed s3 objects.
				# The number of s3Keys will be counted to identify surface topology errors.
				###########################################################################
				if s2s3Key in s3Exposed:
					topologyCheck.update({s2s3Key:None})

		j, edgeMatches = 0, {}
		# Using the list of s3Keys in topologyCheck reduces compute time.
		#################################################################
		for s3Key2 in topologyCheck.keys():

			if s3Key2 == s3Key1:
				continue

			# Surface facets are incident on only one simplex-4.
			####################################################
			if len(s3Exposed[s3Key2]) != 1:
				continue

			# Second facet information.
			###########################
			s32  = tess3D.s3Dict[s3Key2].getFacet()
			s221, s222, s223 = s32.e1, s32.e2, s32.e3
			# Get edges from 3-simplex 2.
			#############################
			eids2 = {s221.id:s221, s221.flipid:s221, 
					s222.id:s222, s222.flipid:s222, 
					s223.id:s223, s223.flipid:s223}

			# Identify if the first and second simplex-3 facets are coincident on a simplex-2.
			##################################################################################
			edgeMatch = None
			for eID in eids2:
				if eID in eids1:
					edgeMatch = eID

			# Count j to find surface s3 facets with abnormal topologies.	
			#
			# Every time there is an edge match, count it by incrementing j.
			# Normal topology for a surface s3 gives a j of 3.
			# Abnormal topologies for a surface s3 gives j > 3.
			#################################################################
			if edgeMatch:
				ei = eids1[edgeMatch]
				ei.lf, ei.rf = s32, s31
				edgeMatches.update({ei.id:ei})
				j += 1

		# The exposed s3 facet as a normal topology:
		# i.e it shares 3 edges (simplex-2 objects) with three neighboring s3 facets. 
		#############################################################################
		if j == 3:
			checkedSurface.update({s3Key1:s3Exposed[s3Key1]})
			#s2Exposed.update({s31:edgeMatches}) # XXXX 3/7

		# The exposed s3 has an abnormal topology:
		# i.e. it is shared by more than 2 surface exposed s3 facets. These facets, 
		# and their corresponding tetrhadral cells will cause rotation errors and need to be deleted. 
		##############################################################################################
		else:
			deleteS3.update({s3Key1:s3Exposed[s3Key1]})

	# This bit of code is key to the proper calculation of a fully 
	# rotatable surface. 3-simplices of s3Exposed that have more than
	# three neighbors will cause rotation errors.  This situation arises
	# when the removal of s4 objects from the working surface reveals
	# two or more new surface s4 objects that share the same 2-simplex edge.
	# This typically occurs in regions where holes and deep crevices 
	# begin to develop.
	#
	# The solution is to identify these newly revealed s4 objects and
	# delete them from the working sub-triangulation (i.e., subtess3D).
	#######################################################################

	# If deleteS3 is populated, there are topology errors. 
	######################################################
	if deleteS3:
		print("Fixing...", len(deleteS3), "...topology errors.")
		result = updateS3CompleteCheckSurface(deleteS3, subTess3D, s3Complete, s3Exposed, s4Exposed)
		return result
	else:
		print("No topology errors...")
		return (deleteS3, s3Complete, s3Exposed, s4Exposed)

def updateS3CompleteCheckSurface(deleteS3, subTess3D, s3Complete, s3Exposed, s4Exposed):

	# A) Identify the s4 objects to be removed.
	###########################################
	deleteS4 = {}
	for s3Key in deleteS3:
		deleteS4.update({deleteS3[s3Key][0]:deleteS3[s3Key][0]})

	# B) Update s4Exposed.
	######################
	for s4 in deleteS4:
		if s4 in s4Exposed:
			del s4Exposed[s4]

	# C) s3Complete is fully linked subTess3D !!!!!
	###############################################
	for s3Key in deleteS3:
		# For s4 with topology errors.
		##############################
		for s4 in s3Complete[s3Key]:
			# Access neighboring s4 objects.
			################################
			for s3id in s4.subKeys:
				# Skip s3 with topology error.
				##############################
				if s3id != s3Key:
					# Update s3Complete.
					####################
					s3Complete[s3id].remove(s4)

				# Empty s3Complete containing the topology error.
				#################################################
				else:
					# Update s3Complete.
					####################
					s3Complete[s3id].remove(s4)

	result = getInitialSurface(subTess3D, s3Complete)
	s3ExposedU, s4ExposedU, s4Isolated, s4Coincident = result[0], result[1], result[2], result[3]
	s3Exposed.update(s3ExposedU)
	s4Exposed.update(s4ExposedU)
		
	return (deleteS3, s3Complete, s3Exposed, s4Exposed)

# Essential function.
#####################
def surfaceGlobal(hull4D, circumSphereRadiusLimit, writeSurfaceCreationAnimation=""):
	
	# 2015.10.17 -- I made a major change to the surface algorithm.
	#
	# There is no longer an iterative trimming using progressively smaller circumsphere radii cutoffs.
	# Also, the iterative removal of surface tetrahedra no longer uses a convex hull topology evaluation.
	# Now, a single circumsphere radii cutoff is used to iteratively remove surface tetrahedra.
	# Now, newly exposed surface edges are identified in each iteration, and their coincident facets are
	# removed if their circumsphere radii are larger than circumSphereRadiusLimit. The while loop exits 
	# when no more facets can be removed.
	#
	# These changes were motivated by differences I observed in the surfaces of the potassium channel (4XDJ)
	# calculated using chains A.B.C.D versus A.B. The A.B. surface was "correct". However, the A.B.C.D was
	# not fully "carved-out". The was caused by the use of the convex hull as a topological discriminator.
	# Because the A.B.C.D structure contained two homodimers separated by a distance, the convex hull was large.
	# This caused the persistance of some large surface tetrahedra lying between the separate dimers to evaluate
	# as interior (because they were far beneath the convex hull based on the facet centroid evaluation).
	# The new algorithm does not use a convex hull, but instead "peels" large surface facets away iteratively.
	# This approach works fantastically, as the results of the A.B.C.D and A.B calculations are identical.

	
	# STEP I
	# Establish an empty dictionary for all s3 objects of subTess3D.
	################################################################
	
	# Smart search: after first time through.
	#########################################
	# Go over full list: first time through.
	# This collects ALL of the s3 edges of the triangulation.
	# In subsequent runs, I just want to collect s3 edges within s3Complete.
	########################################################################
	s3CompleteLocal, subTess3D = {}, hull4D.triangulation
	for s3Key in subTess3D:
		for s4 in subTess3D[s3Key]:
			for s3id in s4.subKeys:
				if s3id not in s3CompleteLocal:
					s3CompleteLocal.update({s3id:[s4]})
				elif s4 not in s3CompleteLocal[s3id]:
					s3CompleteLocal[s3id].append(s4)


	# If the creation of the surface is to be saved, make a new directory to store the creation files.
	##################################################################################################
	if writeSurfaceCreationAnimation:
		
		# Make the needed directory.
		############################
		if not exists(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/"):
			mkdir(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/")


	# STEP II
	# Identify exposed s3 and s4 objects of the triangulation. Also scan for, and remove, isolated s4 pairs.
	########################################################################################################
	results = getInitialSurface(subTess3D, s3CompleteLocal)
	s3Exposed, s4Exposed, s4Isolated, s4Coincident = results[0], results[1], results[2], results[3]

	# Write the first instance of the surface carving animation.
	############################################################
	# These facets correspond to the convex hull.
	#############################################
	if writeSurfaceCreationAnimation:
		iOrientedSurface = orientSurfaceTopology(hull4D.s3Dict, s3Exposed, s4Exposed)#,s4Coincident)
		currentSurface = iOrientedSurface
		writeSubTesselation(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/" + "SubTesselation", str(0), currentSurface, hull4D)


	print("\n\nInitial surface carving...")
	print(80*"-")
	j, facetDeletion, subTess3D, dEdges, edges = 1, 1, {}, {}, hull4D.triangulation
	while facetDeletion:
		
		print("Surface carving iteration...", j)
		j += 1
		
		subTess3D, dEdges, previousS3ExposedLength = {}, {}, len(s3Exposed)
		for s3Key in s3Exposed:
		
			# Access the pair of tetrahedral cells that share the common
			# edge and, for each tetrahedral cell, access its vertices.
			#############################################################
			if len(edges[s3Key]) != 1:
				s41 = edges[s3Key][0]
				s42 = edges[s3Key][1]
				cs1 = s41.getSphere(PseudoAtom())
				cs2 = s42.getSphere(PseudoAtom())
			
			else:
				s41, s42 = edges[s3Key][0], None
				cs1, cs2 = s41.getSphere(PseudoAtom()), None

			# s41 = edges[s3Key][0]
			# cs1 = s41.getSphere(PseudoAtom())
			# if cs1[0] > circumSphereRadiusLimit:
			#     dEdges.update({s3Key:edges[s3Key]})

			# The secret to the surface...
			#############################################################
			if cs1:
				# Only remove tetrahedral cells whose circumspheres...
				#####################################################
				if cs1[0] > circumSphereRadiusLimit:
					dEdges.update({s3Key:edges[s3Key]})
					continue
			if cs2:
				if cs2[0] > circumSphereRadiusLimit:
					dEdges.update({s3Key:edges[s3Key]})
					continue


		# Expose new facets by deleting s4s incident on dEdges.
		#######################################################
		for s3Key in dEdges:
			# s3Key has only one incident tetrahedron.
			for s4 in edges[s3Key]:
				for s3Key2 in s4.subKeys:
					# However, s3Key2(s) may have two incident tetrahedrons.
					if s3Key != s3Key2:
						s4Index = edges[s3Key2].index(s4)
						del edges[s3Key2][s4Index]

		# New 2015.12.02 test: fixes animation
		# Remove emptied surface edges.
		###############################
		for s3Key in edges:
			if not edges[s3Key]:
				dEdges.update({s3Key:edges[s3Key]})
	
		# Carving step: remove the surface tetrahedra for this iteration.
		#################################################################
		for s3Key in dEdges:
			del edges[s3Key]
		subTess3D = edges

		# Update the necessary loop variables.
		######################################
		s3CompleteLocal = {}
		for s3Key in subTess3D:
			for s4 in subTess3D[s3Key]:
				for s3id in s4.subKeys:
					if s3id not in s3CompleteLocal:
						s3CompleteLocal.update({s3id:[s4]})
					elif s4 not in s3CompleteLocal[s3id]:
						s3CompleteLocal[s3id].append(s4)
		results = getInitialSurface(subTess3D, s3CompleteLocal)
		s3Exposed, s4Exposed, s4Isolated, s4Coincident = results[0], results[1], results[2], results[3]

		# Determine if the interations should continue.
		###############################################
		if len(s3Exposed) == previousS3ExposedLength:
			facetDeletion = 0

		# Write the surface carving animation.
		######################################
		if writeSurfaceCreationAnimation:
			iOrientedSurface = orientSurfaceTopology(hull4D.s3Dict, s3Exposed, s4Exposed)#,s4Coincident)
			currentSurface = iOrientedSurface
			writeSubTesselation(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/" + "SubTesselation", str(j-1), currentSurface, hull4D)

	# Set s3Complete.
	#################
	s3Complete = s3CompleteLocal

	# Identify surface from iteration i and make topological corrections.
	#####################################################################
	iSurface = identifySurface(subTess3D,hull4D,s3Complete=s3Complete,checkSurfaceOn=1)
	s3Exposed, s4Exposed, s4Coincident, s3Complete = iSurface[0], iSurface[1], iSurface[2], iSurface[3]

	# After the last surface refinement, orient the facets of the surface. 
	# This is the clever part of the code.  It orients the circuit of each 
	# surface facet such that atoms "below" the plane of the facet register
	# as "inside" the surface, and aotms "above" the plane of the facet
	# register as outside the surface.  See comments within the orientSurfaceTopology()
	# function for more details. 
	###################################################################################
	iOrientedSurface = orientSurfaceTopology(hull4D.s3Dict, s3Exposed, s4Exposed)#,s4Coincident)
	currentSurface = iOrientedSurface

	# Write the surface carving animation.
	######################################
	if writeSurfaceCreationAnimation:
		writeSubTesselation(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/" + "SubTesselation", str(j), currentSurface, hull4D)

	# Surface output.
	#################
	modSurf = {}
	surfaceString = ""
	for s3Key in currentSurface:
		orientedFacet = currentSurface[s3Key]
		modSurf.update({s3Key:orientedFacet})
		surfaceString += s3StringOutput(orientedFacet.v1, orientedFacet.v2, orientedFacet.v3)

	return (modSurf, surfaceString, subTess3D)

##################################################################################
#
# Surface refinement functions.
#
##################################################################################

# Essential function.
#####################
def findSlivers(quadEdges, dCut=1.0, minAngleCut=90.0):

	# Rotate over each origin vertex to find extreme concavities.
	#############################################################
	eKeys = sorted(quadEdges) # XXXX more deterministic for chemiforms?
	# for eID in quadSurface:
	sliverRoots, originVertices = {}, {}
	for eID in eKeys:

		e = quadEdges[eID]
		# if e.origin_vertex.id in originVertices or e.destination_vertex in originVertices:
		# 	continue

		# Get the edges incident on the vertex vPop that will be "popped".
		##################################################################
		vPopEdges = rotateCwAboutTheOriginVertex(e, quadEdges)
		#print("-------------------")
		# pop_edges = 0
		for e_pop in vPopEdges:
			#print(e, e.Oprev)
			angle = dotProduct(e_pop.origin_vertex, e_pop.destination_vertex, quadEdges[e_pop.Oprev].destination_vertex, angle=1)
			#print(e.origin_vertex.id, e.destination_vertex.id, quadEdges[e.Oprev].destination_vertex.id, angle)
			if angle < 60:
				# pop_edges = 1
				# break
				sliverRoots[(e_pop.origin_vertex, e_pop.id)] = e_pop.id
				break

		# if pop_edges:
		# 	for e in vPopEdges:
		# 		sliverRoots[(e.origin_vertex, e.id)] = e.id

		# # originCommonEdges = rotateCwAboutTheOriginVertex(e, quadSurface)
		# originVertex = e.origin_vertex
		# originVertices[e.origin_vertex.id] = None

		# # Collect the repair or "rim" edges that will be used to build the new set of "popped" facets.
		# ##############################################################################################
		# rimEdges, rVertices, aveD, tinyFacet = [], {}, 0., 0
		# for cE in vPopEdges:
		# 	ridgeEdgeID = cE.Dnext 
		# 	rE = quadEdges[ridgeEdgeID]
		# 	rimEdges.append(rE)
		# 	rVertices.update({rE.origin_vertex:None, rE.destination_vertex:None})
		# 	d = distance(rE.origin_vertex, rE.destination_vertex)
		# 	# Skip popping tinyFacets that will create surface edges that violate the zero definition...
		# 	if d < zero: # 2020.05.27 fix...
		# 		tinyFacet = 1
		# 		break
		# 	aveD += distance(rE.origin_vertex, rE.destination_vertex)

		# # The sliver is too tiny to pop , so move on to next sliver...
		# # i.e., popping exceeds the distance limit set by zero...
		# if tinyFacet:
		# 	continue

		# aveD = aveD/float(len(rimEdges))

		# # Calculate the centroid of the rim edges...
		# psa = PseudoAtom()
		# nP = centroid4D(rVertices.keys(), psa)
		# d = distance(nP, originVertex)
		# if (d/aveD) >= dCut:
		# 	sliverRoots[(originVertex, eID)] = eID
		# 	# sliverRoots.append((originVertex, eID)) # XXXX Py2-3

	print("------------------------------->", 1, len(sliverRoots))

	unique_sliverRoots = {}
	for sliverRoot in sliverRoots:
		unique_sliverRoots[sliverRoot[0]] = sliverRoot
	print("------------------------------->", 2, len(unique_sliverRoots))
	sliverRootsFinal = {}
	for unique_sliverRoot in unique_sliverRoots:
		sliverRootsFinal[unique_sliverRoots[unique_sliverRoot]] = sliverRoots[unique_sliverRoots[unique_sliverRoot]]
	sliverRoots = sliverRootsFinal
	print("------------------------------->", 3, len(sliverRoots))
	return sliverRoots

# def findSliversSave(quadEdges, dCut=1.0, minAngleCut=90.0):

# 	# Rotate over each origin vertex to find extreme concavities.
# 	#############################################################
# 	eKeys = sorted(quadEdges) # XXXX more deterministic for chemiforms?
# 	# for eID in quadSurface:
# 	sliverRoots, originVertices = {}, {}
# 	for eID in eKeys:

# 		e = quadEdges[eID]
# 		if e.origin_vertex.id in originVertices or e.destination_vertex in originVertices:
# 			continue

# 		# Get the edges incident on the vertex vPop that will be "popped".
# 		##################################################################
# 		vPopEdges = rotateCwAboutTheOriginVertex(e, quadEdges)
# 		print("-------------------")
# 		for e in vPopEdges:
# 			print(e, e.Oprev)
# 			angle = dotProduct(e.origin_vertex, e.destination_vertex, quadEdges[e.Oprev].destination_vertex, angle=1)
# 			print(e.origin_vertex.id, e.destination_vertex.id, quadEdges[e.Oprev].destination_vertex.id, angle)

# 		# originCommonEdges = rotateCwAboutTheOriginVertex(e, quadSurface)
# 		originVertex = e.origin_vertex
# 		originVertices[e.origin_vertex.id] = None

# 		# Collect the repair or "rim" edges that will be used to build the new set of "popped" facets.
# 		##############################################################################################
# 		rimEdges, rVertices, aveD, tinyFacet = [], {}, 0., 0
# 		for cE in vPopEdges:
# 			ridgeEdgeID = cE.Dnext 
# 			rE = quadEdges[ridgeEdgeID]
# 			rimEdges.append(rE)
# 			rVertices.update({rE.origin_vertex:None, rE.destination_vertex:None})
# 			d = distance(rE.origin_vertex, rE.destination_vertex)
# 			# Skip popping tinyFacets that will create surface edges that violate the zero definition...
# 			if d < zero: # 2020.05.27 fix...
# 				tinyFacet = 1
# 				break
# 			aveD += distance(rE.origin_vertex, rE.destination_vertex)

# 		# The sliver is too tiny to pop , so move on to next sliver...
# 		# i.e., popping exceeds the distance limit set by zero...
# 		if tinyFacet:
# 			continue

# 		aveD = aveD/float(len(rimEdges))

# 		# Calculate the centroid of the rim edges...
# 		psa = PseudoAtom()
# 		nP = centroid4D(rVertices.keys(), psa)
# 		d = distance(nP, originVertex)
# 		if (d/aveD) >= dCut:
# 			sliverRoots[(originVertex, eID)] = eID
# 			# sliverRoots.append((originVertex, eID)) # XXXX Py2-3

# 	return sliverRoots

# Essential function.
#####################
def popDivit(veTuple, identifier, quadEdges):

	# XXXX 2020.05.15 changed...simplified...shortened...testing seems OK...

	# Get the vertex to be popped
	v = veTuple[0]	# vertex to be popped...
	e = quadEdges[veTuple[1]] # first edge...
	vPop = e.origin_vertex

	print(v.id, e)
	print(e.origin_vertex.id, e.origin_vertex.data)
	print(e.destination_vertex.id, e.destination_vertex.data)

	# Get the edges incident on the vertex vPop that will be "popped".
	##################################################################
	vPopEdges = rotateCwAboutTheOriginVertex(e, quadEdges)

	print("-------------------")
	for x in vPopEdges:
		print(x.id)

	# Collect the repair or "rim" edges that will be used to build the new set of "popped" facets.
	##############################################################################################
	# rimEdges, rVertices = [], {}
	# for cE in vPopEdges:
	# 	ridgeEdgeID = cE.Dnext 
	# 	rE = quadEdges[ridgeEdgeID]
	# 	rimEdges.append(rE)
	# 	#rVertices.update({rE.origin_vertex:None, rE.destination_vertex:None})
	# 	rVertices.update({rE.destination_vertex:None})
	rimEdges, rVertices = [], {}
	for rE in vPopEdges:
		# ridgeEdgeID = cE.Dnext 
		# rE = quadEdges[ridgeEdgeID]
		rimEdges.append(rE)
		#rVertices.update({rE.origin_vertex:None, rE.destination_vertex:None})
		rVertices.update({rE.destination_vertex:None})
	print("******")
	xsum, n = 0., 0.
	for x in rVertices:
		print(x.id, x.x, x.y, x.z)
		xsum += x.x
		n += 1
	print("Average x=", xsum/n)

	# Calculate the centroid of the vertices of the facet pair.
	###########################################################
	psa = PseudoAtom()
	psa.atom_serial = identifier
	nP = centroid4D(rVertices.keys(), psa)
	nP.id = identifier
	print(1, nP.x, nP.y, nP.z)
	print(2, nP.data)
	#nP.data = e.destination_vertex.data # Maintain link to surface vertex
	print(3, e.origin_vertex.data)
	print(4, e.destination_vertex.data)
	nP.reinitialize()
	print(5, nP.x, nP.y, nP.z)

	# XXXX 2020.05.16 deemed unecessary...
	# # Ensure general position of nP.
	# ################################
	# vGp = 0
	# while not vGp:
	#     for e in rimEdges:
	#         vGp = gp3D(e.rf.v1, e.rf.v2, e.rf.v3, nP)
	#         if not vGp:
	#             nP.jostle()
	#             break

	# # Replace vR with vPop: the actual "pop" step.
	# ##############################################
	# # The actual pop step is where the problem begins...
	# firstEdge = rimEdges[0]
	# # print(firstEdge.rf.vR.x, firstEdge.rf.vR.y, firstEdge.rf.vR.z, "first values...")
	# firstEdge.rf.vR.x = nP.x
	# firstEdge.rf.vR.y = nP.y
	# firstEdge.rf.vR.z = nP.z
	# firstEdge.rf.vR.u = nP.u  

	# Pop the divit...
	for rimEdge in rimEdges:
		print("popping divit")
		rimEdge.rf.vR = nP	
		quadEdges[rimEdge.id] = rimEdge

# def popDivitSave(veTuple, identifier, quadEdges):

# 	# XXXX 2020.05.15 changed...simplified...shortened...testing seems OK...

# 	# Should be get the vertex to be popped

# 	# Get the edge to be "popped".
# 	###############################
# 	v = veTuple[0]
# 	e = quadEdges[veTuple[1]]
# 	vPop = e.origin_vertex

# 	# Get the edges incident on the vertex vPop that will be "popped".
# 	##################################################################
# 	vPopEdges = rotateCwAboutTheOriginVertex(e, quadEdges)

# 	# Collect the repair or "rim" edges that will be used to build the new set of "popped" facets.
# 	##############################################################################################
# 	rimEdges, rVertices = [], {}
# 	for cE in vPopEdges:
# 		ridgeEdgeID = cE.Dnext 
# 		rE = quadEdges[ridgeEdgeID]
# 		rimEdges.append(rE)
# 		rVertices.update({rE.origin_vertex:None, rE.destination_vertex:None})
# 	print(rVertices)

# 	# Calculate the centroid of the vertices of the facet pair.
# 	###########################################################
# 	psa = PseudoAtom()
# 	psa.atom_serial = identifier
# 	nP = centroid4D(rVertices.keys(), psa)
# 	nP.id = identifier
# 	print(1, nP.x, nP.y, nP.z)
# 	nP.data = e.destination_vertex.data # Maintain link to surface vertex
# 	print(2, e.destination_vertex.data)
# 	nP.reinitialize()
# 	print(3, nP.x, nP.y, nP.z)

# 	# XXXX 2020.05.16 deemed unecessary...
# 	# # Ensure general position of nP.
# 	# ################################
# 	# vGp = 0
# 	# while not vGp:
# 	#     for e in rimEdges:
# 	#         vGp = gp3D(e.rf.v1, e.rf.v2, e.rf.v3, nP)
# 	#         if not vGp:
# 	#             nP.jostle()
# 	#             break

# 	# # Replace vR with vPop: the actual "pop" step.
# 	# ##############################################
# 	# # The actual pop step is where the problem begins...
# 	# firstEdge = rimEdges[0]
# 	# # print(firstEdge.rf.vR.x, firstEdge.rf.vR.y, firstEdge.rf.vR.z, "first values...")
# 	# firstEdge.rf.vR.x = nP.x
# 	# firstEdge.rf.vR.y = nP.y
# 	# firstEdge.rf.vR.z = nP.z
# 	# firstEdge.rf.vR.u = nP.u  

# 	# Pop the divit...
# 	for rimEdge in rimEdges:
# 		print("popping divit")
# 		rimEdge.rf.vR = nP	
# 		quadEdges[rimEdge.id] = rimEdge

# Essential function.
#####################
def dividePatch(quadEdge, identifier, quadEdges):

	e = quadEdge

	# Calculate the midpoint (i.e. centroid) of the divide path center, and create a new vertex4D object via centroid4D().
	# Using these 4 points results in a "smoother" surface and less jostle calls...
	######################################################################################################################
	psa = PseudoAtom()
	nP = centroid4D([quadEdge.destination_vertex, quadEdge.origin_vertex, quadEdge.lf.vR, quadEdge.rf.vR], psa)
	nP.id = identifier
	nP.data = quadEdge.destination_vertex.data # Maintain link to surface vertex
	#print("B", nP.id, nP.x, nP.y, nP.z, nP.data)
	nP.reinitialize()
	#print("A", nP.id, nP.x, nP.y, nP.z, nP.data)

	# XXXX 2020.05.16 deemed unecessary...
	# # Ensure general position of nP.
	# ################################
	# # May have helped...but still issues????? dividePatch is the problem...gp may be involved too...
	# vGp = 0
	# while not vGp:
	#     for f in [e.lf, e.rf]:
	#         vGp = gp3D(f.v1, f.v2, f.v3, nP)
	#         if not vGp:
	#             nP.jostle()
	#             break

	# Identify repair edges.
	########################
	rightRepairEdges = rotateCwOverTheRightFace(quadEdge, quadEdges)[1:]
	leftRepairEdges  = rotateCwOverTheLeftFace(quadEdge, quadEdges)[1:]
	repairEdges = rightRepairEdges + leftRepairEdges

	# 2020.05.16
	# Write now, I am skipping the test of general position because I think it is not needed...

	# # If nP is not in gP then dividePatch has met its natural terminal end.
	# # In the other words, the facet cannot be sub-divided further without 
	# # perpetually violating the algorithm's defined zero value.
	# ########################################################################
	# vGp = 1
	# for e in repairEdges:
	#     vGp = gp4D(e.rf.v1, e.rf.v2, e.rf.v3, nP, e.lf.v3)
	#     if not vGp:
	#         vGp = 0

	# # General position is NOT present and the sub-division process will terminate. 
	# ##############################################################################
	# if not vGp:
	#     print("quadEdge cannot be divided while maintaining general position...return as longEdgesToSkip...")
	#     return quadEdge.id

	# Quad-edge origin and destination vertices.
	############################################
	vO = quadEdge.origin_vertex
	vD = quadEdge.destination_vertex
	patchConeEdges = {}
	for e in repairEdges:

		# Algorithmically made to match popDivit.
		#########################################
		if e.lf.vR.id in [vO, vD]:
			pass
		else:
			e = quadEdges[e.flipid]

		# Update the edge and the quad surface.
		#######################################
		e.rf.vR = nP
		#e.rf.make_ltriangle()
		e.rf.update_triangle() # XXXX this fixes None vR problem by updating edge topology....
		e.lf = e.rf.ltriangle
		quadEdges[e.id] = e
		quadEdges[e.flipid] = e.rf.ltriangle.e1

		# Identify the coneFacet.
		#########################
		coneFacet = None
		if nP.id in e.rf.id:
			e.rf.e1 = e
			coneFacet = e.rf
		elif nP.id in e.lf.id:
			e.lf.e1 = e
			coneFacet = e.lf

		# Collect the edges of the cone facet, exclusive of the repair edge.
		####################################################################
		if coneFacet.e2.id not in patchConeEdges:
			patchConeEdges[coneFacet.e2.id] = {coneFacet.id:coneFacet}
		else:
			patchConeEdges[coneFacet.e2.id].update({coneFacet.id:coneFacet})

		if coneFacet.e2.flipid not in patchConeEdges:
			patchConeEdges[coneFacet.e2.flipid] = {coneFacet.id:coneFacet}
		else:
			patchConeEdges[coneFacet.e2.flipid].update({coneFacet.id:coneFacet})

		if coneFacet.e3.id not in patchConeEdges:
			patchConeEdges[coneFacet.e3.id] = {coneFacet.id:coneFacet}
		else:
			patchConeEdges[coneFacet.e3.id].update({coneFacet.id:coneFacet})

		if coneFacet.e3.flipid not in patchConeEdges:
			patchConeEdges[coneFacet.e3.flipid] = {coneFacet.id:coneFacet}
		else:
			patchConeEdges[coneFacet.e3.flipid].update({coneFacet.id:coneFacet})

	# Convert patchConeEdges into a format compatible with getQuadSurface.
	######################################################################
	matchedFacets = {}
	for e in patchConeEdges:
		matchedFacets.update({e:list(patchConeEdges[e].values())}) # XXXX added list()

	# Add updated repair edges to matchedFacets.
	############################################
	for e in repairEdges:
		e = quadEdges[e.id] 
		# Must be arranged as (lf, rf): topologically opposite to popDivit.
		###################################################################
		matchedFacets.update({e.id:(e.lf, e.rf)})

	# Update the quad-edge surface with the new patch cone quad edges. 
	##################################################################
	quadEdges.update(getQuadSurface(matchedFacets))

	# Delete the original quad edge of patch.
	#########################################
	del quadEdges[quadEdge.id]
	del quadEdges[quadEdge.flipid]

# Essential function.
#####################
def smoothen(quadEdge, quadEdges, already_smoothened_vertices):

	if quadEdge.origin_vertex in already_smoothened_vertices:
		return 0

	# Collect the cone edges coincident on the origin_vertex of the quadEdge...
	coneEdges = rotateCwAboutTheOriginVertex(quadEdge, quadEdges)
	
	# Collect the destination_vertex of each facet of the cone rim...
	coneRimVertices = {}
	for e in coneEdges:
		coneRimVertices.update({e.destination_vertex:None})

	# Calculate the centroid of the rim...
	psa = PseudoAtom()
	psa.atom_serial = quadEdge.origin_vertex.data.atom_serial
	nP = centroid4D(coneRimVertices.keys(), psa)

	# Update the position of the origin_vertex of the quadEdge...
	quadEdge.origin_vertex.x, quadEdge.origin_vertex.y, quadEdge.origin_vertex.z = nP.x, nP.y, nP.z 
	quadEdge.origin_vertex.reinitialize()

# def smoothenSave(quadEdge, quadEdges):

# 	coneEdges = rotateCwAboutTheOriginVertex(quadEdge, quadEdges)[1:]

# 	patchRingVertices = []
# 	for e in coneEdges:
# 		if e.destination_vertex not in patchRingVertices:
# 			patchRingVertices.append(e.destination_vertex)
# 	patchCentroid = centroid4D(patchRingVertices, PseudoAtom())
# 	patchCentroid.data = coneEdges[0].origin_vertex.data # Maintain link to surface vertex

# 	volumes, totalVolume = [], 0.
# 	for e in coneEdges:

# 		patchRingEdge = rotateCwOverTheRightFace(e, quadEdges)[1]
# 		patchTriangle = Triangle(patchRingEdge.origin_vertex, patchRingEdge.destination_vertex, patchCentroid, skip_orientation=1)
# 		# The patchTriangle cannot achieve gp because the zero definition will be violated...
# 		# Kill the attempt to smoothen this edge...
# 		# 2020.05.27...
# 		if not patchTriangle.gp:
# 			return 0
# 		volume = patchTriangle.beneath_beyond(quadEdge.origin_vertex, returnVolumeTuple=1)
# 		volumes.append(volume[0])
# 		totalVolume += volume[1]

# 	d = distance(quadEdge.origin_vertex, patchCentroid)

# 	# This condition and code improve the smoothness of the surface...
# 	# However, the addition of this codes slows the algorithm...
# 	# print("Smoothing surface...")
# 	# print(80*"-")
# 	#if 1 not in volumes and d > 1.0 and totalVolume > 0.1:
# 	# This condition prevents internal slivers...
# 	if d > 1.0 and totalVolume > 1e-3:

# 		quadEdge.origin_vertex.x, quadEdge.origin_vertex.y, quadEdge.origin_vertex.z = patchCentroid.x, patchCentroid.y, patchCentroid.z 
# 		quadEdge.origin_vertex.reinitialize()

# 		for e in coneEdges:
# 			patchRingEdge = rotateCwOverTheRightFace(e, quadEdges)[1]
# 			patchRingEdge.lf.vR.x, patchRingEdge.lf.vR.y, patchRingEdge.lf.vR.z = patchCentroid.x, patchCentroid.y, patchCentroid.z
# 			patchRingEdge.lf.vR.reinitialize()


# 		surfaceFacetsString = ""
# 		for e in coneEdges:

# 			format = "%15.6f%15.6f%15.6f "
# 			surfaceFacetsString += format % (e.lf.v1.x, e.lf.v1.y, e.lf.v1.z) + \
# 									format % (e.lf.v2.x, e.lf.v2.y, e.lf.v2.z) + \
# 									format % (e.lf.v3.x, e.lf.v3.y, e.lf.v3.z) + "\n"
# 			surfaceFacetsString += format % (e.rf.v1.x, e.rf.v1.y, e.rf.v1.z) + \
# 									format % (e.rf.v2.x, e.rf.v2.y, e.rf.v2.z) + \
# 									format % (e.rf.v3.x, e.rf.v3.y, e.rf.v3.z) + "\n"

# 		return surfaceFacetsString

# 	return 0

# Essential function.
#####################
def facetAreas(quadEdges, minArea=20.0):

	areaEdges = {}
	eKeys = sorted(quadEdges)
	for e in eKeys:
		qE = quadEdges[e]
		a1 = qE.lf.area()
		a2 = qE.rf.area()
		if a1 > minArea or a2 > minArea:	
			areaEdges[e] = None	
	return list(areaEdges)

# # Essential function.
# #####################
# def longEdgeCut(quadSurface, longCut, longEdgesToSkip={}):

# 	# I am thinking that on the back end of this function in dividePatch,
# 	# In situations in which there is large patch of long edges, the dividePatch
# 	# function may not be able to properly build and set the topology of the resultant surface, 
# 	# perhaps resulting in hung edges. 
# 	#
# 	# I am trying a solution in which only one edge of a triangle facet can be split at a time. 
# 	#
# 	# This situation really seems to be limited to small surfaces...but I cannot guarantee it hasn't
# 	# influenced cluster calculations in the past...I need to investigate.

# 	# longVids, longEdges, longEdgeFacets, longEdgeIDs = {}, {}, {}, []
# 	# eKeys = sorted(quadSurface) # XXXX more deterministic for chemiforms?
# 	# for eID in eKeys:

# 	# 	if eID in longEdgesToSkip:
# 	# 		continue

# 	# 	if eID in longEdges:
# 	# 		continue

# 	# 	# A surface edge...
# 	# 	e = quadSurface[eID]

# 	# 	# The left facet of the edge...
# 	# 	maxD, maxE = 0., None
# 	# 	for leTest in [e.lf.e1, e.lf.e2, e.lf.e3]:
# 	# 		d = distance(leTest.origin_vertex, leTest.destination_vertex)
# 	# 		if d > maxD:
# 	# 			maxD = d
# 	# 			maxE = leTest
# 	# 	if maxD > longCut:
# 	# 		#print("hello-------------------------------------------------------------->")
# 	# 		if maxE.origin_vertex.id not in longVids and maxE.destination_vertex.id not in longVids:
# 	# 			if maxE.id not in longEdges and maxE.flipid not in longEdges:
# 	# 				if e.lf.id not in longEdgeFacets:
# 	# 					#print("hello-------------------------------------------------------------->")
# 	# 					longEdgeIDs.append(maxE.id)
# 	# 					longVids[maxE.origin_vertex.id] = None
# 	# 					longVids[maxE.destination_vertex.id] = None
# 	# 					longEdges[maxE.id] = None
# 	# 					longEdges[maxE.flipid] = None
# 	# 					longEdgeFacets[e.lf.id] = None

# 	# 	# The right facet of the edge...
# 	# 	maxD, maxE = 0., None
# 	# 	for leTest in [e.rf.e1, e.rf.e2, e.rf.e3]:
# 	# 		d = distance(leTest.origin_vertex, leTest.destination_vertex)
# 	# 		if d > maxD:
# 	# 			maxD = d
# 	# 			maxE = leTest
# 	# 	if maxD > longCut:
# 	# 		if maxE.origin_vertex.id not in longVids and maxE.destination_vertex.id not in longVids:
# 	# 			if maxE.id not in longEdges and maxE.flipid not in longEdges:
# 	# 				if e.rf.id not in longEdgeFacets:
# 	# 					#print( "hello-------------------------------------------------------------->")
# 	# 					longEdgeIDs.append(maxE.id)
# 	# 					longVids[maxE.origin_vertex.id] = None
# 	# 					longVids[maxE.destination_vertex.id] = None
# 	# 					longEdges[maxE.id] = None
# 	# 					longEdges[maxE.flipid] = None
# 	# 					longEdgeFacets[e.rf.id] = None

# 	# # print(sorted(longEdgeIDs))
# 	# # print(sorted(longEdgeFacets.keys()))
# 	# return longEdgeIDs

# 	#2020.05.16 ... Maybe old version...nope, this version is faster and works better. 
# 	# Turns out maxID was the problem affecting all of the surface refinement stability and determinism...
# 	longEdgeIDs, alreadyScannedEdges = [], {}
# 	eKeys = sorted(quadSurface) # XXXX more deterministic for chemiforms?
# 	# for eID in quadSurface:
# 	for eID in eKeys:

# 		# Skip edges that were already scanned in a previous round.
# 		if eID in longEdgesToSkip or (eID[1], eID[0]) in longEdgesToSkip:
# 			continue 

# 		# Skip edges that were already scanned in this round.
# 		if eID in alreadyScannedEdges or (eID[1], eID[0]) in alreadyScannedEdges:
# 			continue 

# 		alreadyScannedEdges[eID] = None
# 		alreadyScannedEdges[(eID[1], eID[0])] = None

# 		# # Skip edges that were already scanned.
# 		# if eID in longEdgeIDs or (eID[1], eID[0]) in longEdgeIDs:
# 		# 	continue

# 		# # Skip any edges that should not be considered.
# 		# if eID in longEdgesToSkip or (eID[1], eID[0]) in longEdgesToSkip:
# 		# 	continue

# 		e = quadSurface[eID]
# 		# lf, rf = e.lf, e.rf
# 		d = distance(e.origin_vertex, e.destination_vertex)
# 		if d >= longCut:# and (eID[1], eID[0]) not in longEdgeIDs:
# 			#print(d, lf.area(), rf.area())
# 			# dl1 = distance(lf.e2.origin_vertex, lf.e2.destination_vertex)
# 			# dl2 = distance(lf.e3.origin_vertex, lf.e3.destination_vertex)
# 			# dr1 = distance(rf.e2.origin_vertex, rf.e2.destination_vertex)
# 			# dr2 = distance(rf.e3.origin_vertex, rf.e3.destination_vertex)
# 			# print(dl1, dl2, dr1, dr2)
# 			# from compGeometry import circumSphere
# 			# cs = circumSphere(lf.v1, lf.v2, lf.v3, lf.vR, PseudoAtom())
# 			# print(d, e, cs)
# 			longEdgeIDs.append(eID)

# 		# Simplified 2018.04.30.
# 		# if d < longCut:
# 		# 	continue
# 		# dl1 = distance(lf.e2.origin_vertex, lf.e2.destination_vertex)
# 		# dl2 = distance(lf.e3.origin_vertex, lf.e3.destination_vertex)
# 		# dr1 = distance(rf.e2.origin_vertex, rf.e2.destination_vertex)
# 		# dr2 = distance(rf.e3.origin_vertex, rf.e3.destination_vertex)
# 		# act = 0
# 		# for dTest in [dl1, dl2, dr1, dr2]:
# 		# 	if dTest >= d:
# 		# 		act = 1
# 		# 		break
# 		# if not act and (eID[1], eID[0]) not in longEdgeIDs: # XXXX (eID[1], eID[0]) not in longEdgeIDs # added to prevent duplicate calculations with symmetric edge.
# 		# 	longEdgeIDs.append(eID)	
# 	return longEdgeIDs

##################################################################################
# 
# Functions that operate on the pHinder surface.
#
##################################################################################

# Essential function.
#####################
def makeEdgeDict(orientedSurface):

	surfaceEdges = {}
	for s3Key1 in orientedSurface:

		# First facet.
		##############
		facet1 = orientedSurface[s3Key1]
		e11, e12, e13 = facet1.e1, facet1.e2, facet1.e3
		eids1 = {e11.id:e11, e11.flipid:e11, e12.id:e12, e12.flipid:e12, e13.id:e13, e13.flipid:e13}

		coincidentEdges = 0
		for s3Key2 in orientedSurface:
			if s3Key1 == s3Key2:
				continue

			# Second facet.
			###############
			facet2 = orientedSurface[s3Key2]
			e21, e22, e23 = facet2.e1, facet2.e2, facet2.e3
			eids2 = {e21.id:e21, e21.flipid:e21, e22.id:e22, e22.flipid:e22, e23.id:e23, e23.flipid:e23}

			# Identify if the first and second facets are coincident.
			#########################################################
			edgeMatch = None
			for e in eids2.keys():
				if e in eids1.keys():
					edgeMatch = eids1[e]
					coincidentEdges += 1

			if edgeMatch:
				if edgeMatch not in surfaceEdges:
					surfaceEdges.update({edgeMatch:[facet1,facet2]})
				# XXXX 3/7
				# elif facet1 not in surfaceEdges[edgeMatch]:
				# 	surfaceEdges[edgeMatch].append(facet1)
				# elif facet2 not in surfaceEdges[edgeMatch]:
				# 	surfaceEdges[edgeMatch].append(facet2)
				# Each facet can only gave 3 coincident edges...
				if coincidentEdges == 3: 
					break

	return surfaceEdges

# Essential function.
#####################
def getQuadSurface(surfaceEdges):

	####################################################
	# Note: This function is capable of altering vertex 
	# coordinates to achieve general position.
	####################################################

	quadSurface = {}
	for e in surfaceEdges:

		# # XXXX Sloppy:
		# # Sometimes I pass the edge object to this function,
		# # and sometimes I pass the edge key. Fix.  An example
		# # is in dividePatch.
		# #####################################################
		# eID, fID = None, None
		# try:
		# 	eID, fID = e.id, e.flipid
		# except:
		# 	eID, fID = e, (e[1],e[0])           

		# # BIG ALGORITHM TIME-SAVER.
		# # Skip edges already encountered to avoid unnecessary
		# # volume calculations, which are a huge drag on the 
		# # speed performance of the algorithm.
		# #####################################################
		# if eID in quadSurface or fID in quadSurface:
		# 	continue

		# Get the information that comprises the quad edge.
		#####################################################
		facet1 = surfaceEdges[e][0]
		facet2 = surfaceEdges[e][1]
		edges1 = (facet1.e1, facet1.e2, facet1.e3)
		edges2 = (facet2.e1, facet2.e2, facet2.e3)
		vertices1 = (facet1.v1, facet1.v2, facet1.v3)
		vertices2 = (facet2.v1, facet2.v2, facet2.v3)

		# Find the edge shared by the facets of the quad edge.
		######################################################
		eMatch = None
		for e1 in edges1:
			for e2 in edges2:
				if e1.id in [e2.id, e2.flipid]:
					eMatch = e1
					break

		# Find the first vertex that is not incident on the quad edge.
		##############################################################
		misMatchV1 = None
		for v in vertices1:
			if v not in [eMatch.origin_vertex, eMatch.destination_vertex]:
				misMatchV1 = v
				break

		# Find the second vertex that is not incident on the quad edge.
		###############################################################
		misMatchV2 = None
		for v in vertices2:
			if v not in [eMatch.origin_vertex, eMatch.destination_vertex]:
				misMatchV2 = v
				break

		# Make quad-edge with correctly oriented facets.
		################################################
		newTriangle1 = Triangle(eMatch.destination_vertex, eMatch.origin_vertex, misMatchV2, VertexR=misMatchV1, skip_orientation=1)
		newTriangle1.make_ltriangle()
		newTriangle1.set_edge_topology()
		quadSurface.update({newTriangle1.e1.id:newTriangle1.e1})
		quadSurface.update({newTriangle1.ltriangle.e1.id: newTriangle1.ltriangle.e1})

	return quadSurface

def rotateCwOverTheLeftFace(startEdge, edgeDict):

	# Need to flip the orientation of the start edge. This is unique to the left facet.  
	###################################################################################
	startEdge = edgeDict[(startEdge.id[1], startEdge.id[0])]

	rotate = [startEdge,]
	nextEdge = startEdge.Rprev
	n = 0
	while nextEdge != startEdge.id:

		rotate.append(edgeDict[nextEdge])
		nextEdge = edgeDict[nextEdge].Rprev

	return rotate  

def rotateCwOverTheRightFace(startEdge, edgeDict):

	rotate = [startEdge,]
	nextEdge = startEdge.Rprev
	n = 0
	while nextEdge != startEdge.id:

		rotate.append(edgeDict[nextEdge])
		nextEdge = edgeDict[nextEdge].Rprev

	return rotate   		

def rotateCwAboutTheOriginVertex(startEdge, edgeDict):

	# Rotate ccw over the edges that share the same origin vertex.
	##############################################################
	rotate = [startEdge,]
	next_edge = startEdge.Oprev
	while next_edge != startEdge.id:
		rotate.append(edgeDict[next_edge])
		next_edge = edgeDict[next_edge].Oprev 

	return rotate

def getClosestQuadEdges(vertexList, quadSurface):
	
	# Identify the quad edge closest to each test vertex.
	#####################################################
	eTouched, eMidpoints = {}, {}
	for eID in quadSurface:
		eTouched.update({eID:eID})
		if (eID[1], eID[0]) in eTouched:
			continue
		e = quadSurface[eID]
		eMidpoint = centroid3D([e.origin_vertex, e.destination_vertex], PseudoAtom())
		eMidpoints.update({eMidpoint:e})
	
	closestEdges = {}
	for vertex in vertexList:
		minD, closestEdge = 1000000., None
		for eMidpoint in eMidpoints:
			d = distance(vertex,eMidpoint)
			if d < minD:
				minD = d
				closestEdge = eMidpoints[eMidpoint]
		closestEdges.update({vertex:closestEdge})

	return closestEdges

def getSurroundingSurfaceCentroids(vertex, closestEdge, quadSurface):
	
	# First iteration over the left and right facet of the start edge.
	##################################################################
	localEdges = []
	rotate1L = rotateCwOverTheLeftFace(quadSurface[closestEdge.id], quadSurface)
	localEdges += rotate1L
	rotate1R = rotateCwOverTheRightFace(quadSurface[closestEdge.id], quadSurface)
	localEdges += rotate1R
	
	# Second iteration over the left facet.
	#######################################
	rotate2L, rotate2R = None, None
	for rE1L in rotate1L:
		rotate2L = rotateCwOverTheLeftFace(quadSurface[rE1L.id], quadSurface)
		localEdges += rotate2L
		rotate2R = rotateCwOverTheRightFace(quadSurface[rE1L.id], quadSurface)
		localEdges += rotate2R
	
	# Second iteration over the right facet.
	########################################
	for rE1R in rotate1R:
		rotate2R += rotateCwOverTheRightFace(quadSurface[rE1R.id], quadSurface)
		localEdges += rotate2R
		rotate2L += rotateCwOverTheLeftFace(quadSurface[rE1R.id], quadSurface)
		localEdges += rotate2L

#    # Third iteration over the left facet.
#    ######################################
#    for rE2L in rotate2L:
#        rotate3L = rotateCwOverTheLeftFace(quadSurface[rE2L.id], quadSurface)
#        localEdges += rotate3L
#        rotate3R = rotateCwOverTheRightFace(quadSurface[rE2L.id], quadSurface)
#        localEdges += rotate3R
#    
#    # Third iteration over the right facet.
#    #######################################
#    for rE2R in rotate2R:
#        rotate3R += rotateCwOverTheRightFace(quadSurface[rE2R.id], quadSurface)
#        localEdges += rotate3R
#        rotate3L += rotateCwOverTheLeftFace(quadSurface[rE2R.id], quadSurface)
#        localEdges += rotate3L

	# From the list of local quadEdges, make a unique list of facets.
	#################################################################
	localSurface = {}
	for e in localEdges:
		localSurface.update({e.lf.id:e.lf, e.rf.id:e.rf})
	
	# Calculate the centroid of each facet that comprises the local surface.
	########################################################################
	surroundingCentroids = {}
	for fKey in localSurface:

		# Access the left and right Simplex-3 facets stored by Simplex-3 edges of the quadSurface.
		##########################################################################################
		facet = localSurface[fKey]
			
		# Calculate the centroids in the plane of each Simplex-3 facet.
		###############################################################
		facetCentroid = centroid3D((facet.v1, facet.v2, facet.v3), PseudoAtom())
		if facet.beneath_beyond(vertex) < 0:
			facetSign = 1
		else:
			facetSign = -1
		surroundingCentroids.update({facet.id:(facetCentroid, facetSign)})

	# Return the local facet centroid set that will be used to assess the depth of side chain burial.
	#################################################################################################
	return surroundingCentroids

def generateVirtualScreeningVertices(quadSurface, inIterations=1, inIterationsStepSize=1, outIterations=1, outIterationsStepSize=1, perturbLevel=1):

	# A key to this function is the jostling of the sample vertices. 
	# By jostling the vertices when they are created, they are less likely
	# to trigger additional functions that ensure 4D general position.
	##################################################################

	# Get a unique dictionary of the surface facets to avoid duplicate points.
	##########################################################################
	vertices, facets = {}, {}
	for fKey in quadSurface:
		
		# Access the left and right Simplex-3 facets stored by Simplex-3 edges of the quadSurface.
		##########################################################################################
		leftFacet, rightFacet = quadSurface[fKey].lf, quadSurface[fKey].rf
		if leftFacet.id not in facets:
			facets[leftFacet.id] = leftFacet
			vertices[leftFacet.v1] = leftFacet.v1
			vertices[leftFacet.v2] = leftFacet.v2
			vertices[leftFacet.v3] = leftFacet.v3
		if rightFacet.id not in facets:
			facets[rightFacet.id] = rightFacet
			vertices[rightFacet.v1] = rightFacet.v1
			vertices[rightFacet.v2] = rightFacet.v2
			vertices[rightFacet.v3] = rightFacet.v3

	# Iteratively calculate the virtual screening vertices.
	#######################################################
	virtualScreeningVertices = []
	for fKey in facets:
		
		facet = facets[fKey]

		centroid = centroid4D((facet.v1, facet.v2, facet.v3), PseudoAtom())

		# # This appends the centroid facet...
		# virtualScreeningVertices.append(centroid)

		# i, inStep = 0, inIterationsStepSize
		# while i < inIterations:
		#     v = crossProductCoordinates(centroid, facet.v2, facet.v3, PseudoAtom(), multiplier=float(inStep), norm=1).v
		#     v.jostle()
		#     virtualScreeningVertices.append(v)
		#     inStep += inIterationsStepSize
		#     i += 1

		# i, outStep = 0, outIterationsStepSize
		# while i < outIterations:
		#     # Swap v2 and v3 so the orthogonal result of the cross product is outside of the facet.
		#     #######################################################################################
		#     v = crossProductCoordinates(centroid, facet.v3, facet.v2, PseudoAtom(), multiplier=float(outStep), norm=1).v
		#     v.jostle()
		#     virtualScreeningVertices.append(v)
		#     print(i, outStep, outIterations)
		#     outStep += outIterationsStepSize
		#     i += 1 

		# Swap v2 and v3 so the orthogonal result of the cross product is outside of the facet.
		#######################################################################################

		#XXXX apparently choice for chemiforms
		# outStep = outIterationsStepSize
		# v = crossProductCoordinates(centroid, facet.v3, facet.v2, PseudoAtom(), multiplier=float(outIterationsStepSize), norm=1).v
		# v.jostle()
		# virtualScreeningVertices.append(v)
		
		# # Original....
		# centroid = centroid4D((facet.v1, facet.v2, facet.v3), PseudoAtom())
		# virtualScreeningVertices.append(centroid)

		# i = inIterationsStepSize
		# while i <= inIterations:
		#     v = crossProductCoordinates(centroid, facet.v2, facet.v3, PseudoAtom(), multiplier=float(i), norm=1).v
		#     v.jostle()
		#     virtualScreeningVertices.append(v)
		#     i += inIterationsStepSize

		# i = outIterationsStepSize
		# while i <= outIterations:
		#     # Swap v2 and v3 so the orthogonal result of the cross product is outside of the facet.
		#     #######################################################################################
		#     v = crossProductCoordinates(centroid, facet.v3, facet.v2, PseudoAtom(), multiplier=float(i), norm=1).v
		#     v.jostle()
		#     virtualScreeningVertices.append(v)
		#     i += outIterationsStepSize 

		# XXXX
		if outIterations:
			i, step = 0, outIterationsStepSize
			while i <= outIterations:
				# Swap v2 and v3 so the orthogonal result of the cross product is outside of the facet.
				#######################################################################################
				v = crossProductCoordinates(centroid, facet.v3, facet.v2, PseudoAtom(), multiplier=float(step), norm=1).v
				# 2021.01.28 The perturbLevel argument is intended to provide sufficient coordinate variability
				# among sampling points generated along the facet normal using the crossProductCoordinates function.
				# 2020.10.19 Need perturn level 2 to reduce number of jostles necessary in the triangulation, 
				# otherwise go with 1. A perturb level of 2 is needed in the case of multiple iterations.
				v.jostle(perturbLevel=perturbLevel) 
				virtualScreeningVertices.append(v)
				step += outIterationsStepSize
				i += 1

		# XXXX
		if inIterations:
			i, step = 0, inIterationsStepSize
			while i <= inIterations:
				# Swap v2 and v3 so the orthogonal result of the cross product is outside of the facet.
				#######################################################################################
				v = crossProductCoordinates(centroid, facet.v3, facet.v2, PseudoAtom(), multiplier=-float(step), norm=1).v
				# 2021.01.28 The perturbLevel argument is intended to provide sufficient coordinate variability
				# among sampling points generated along the facet normal using the crossProductCoordinates function.
				# 2020.10.19 Need perturn level 2 to reduce number of jostles necessary in the triangulation, 
				# otherwise go with 1. A perturb level of 2 is needed in the case of multiple iterations.
				#v.jostle(perturbLevel=perturbLevel) 
				v.jostle(perturbLevel=perturbLevel)
				virtualScreeningVertices.append(v)
				step += inIterationsStepSize
				i += 1

	# # Extend the virtual screening vertices, calculated from the facet centroids, with the original surface vertices.
	# #################################################################################################################
	# virtualScreeningVertices.extend(vertices.values())  

	# Index virtual screening vertices by atom_serial number.
	#########################################################
	i, indexedVirtualScreeningVertices = 0, []
	for v in virtualScreeningVertices:
		v.id = i
		v.data.atom_serial = i
		v.data.reinitialize()
		indexedVirtualScreeningVertices.append(v)
		i += 1

	print("Initial number of virtual screening vertices:", len(indexedVirtualScreeningVertices))
	return indexedVirtualScreeningVertices

def inLocalSurface(vertexList, quadSurface, coreCutoff, marginCutoff):


	# print("Start getClosestQuadEdges....")
	closestQuadEdges = getClosestQuadEdges(vertexList, quadSurface)
	# print("End   getClosestQuadEdges....")

	classifications = []
	for v in vertexList:

		# Collect the surface centroids closest to v.
		#############################################
		surroundingSurfaceCentroids = getSurroundingSurfaceCentroids(v, closestQuadEdges[v], quadSurface)

		# Calculate the average distance from v the surface centroids.
		##############################################################
		minDistance, aveDistance, nDistances, distances = 1000000., 0., 0., []
		for facetKey in surroundingSurfaceCentroids:
			centroid = surroundingSurfaceCentroids[facetKey][0]
			sign = surroundingSurfaceCentroids[facetKey][1]
			d = distance(v, centroid)
			# if v.data.residue_sequence_number == 253:
			# 	print(d, sign*d)
			if d < minDistance:
				minDistance = d
			distances.append(sign*d)
			aveDistance += (sign*d)
			nDistances += 1.
		aveDistance = aveDistance/nDistances

		# print(v.data.residue_sequence_number)
		# if v.data.residue_sequence_number == 253:
		# 	print("--------------------------------------------> ave distance", aveDistance)

		if aveDistance <= coreCutoff:

			# # Original classification code as of 2023.03.24...

			# # Added 2019.04.23 to properly identify deepy buried side chains...
			# numberBelowCoreCutoff = 0.
			# for d in distances:
			# 	if d < coreCutoff:
			# 		numberBelowCoreCutoff += 1.0
			# coreScore = numberBelowCoreCutoff/float(len(distances))

			# # Condition added 2019.04.15 to properly identify exposed side chains that "appear" deeply buried...
			# # A very deeply buried side chain that is distant from surface facets...
			# if coreScore > 0.6:
			# 	classifications.append((v, -1, aveDistance))
			# 	#print("core", coreScore, (max(distances) - min(distances)), v.data)
			# # A side chain that is very exposed and distant from any surface facets...
			# elif (max(distances) - min(distances)) > 5.0:
			# 	# print(coreScore, distances)
			# 	# print("exposed", (max(distances) - min(distances)), v.data)
			# 	classifications.append((v, 1, aveDistance))
			# # A moderately buried side chain that is relatively close to the surface facets...
			# else:
			# 	classifications.append((v, -1, aveDistance))
			# 	# print(coreScore, "core", (max(distances) - min(distances)))


			# Added 2019.04.23 to properly identify deepy buried side chains...
			numberBelowCoreCutoff = 0.
			for d in distances:
				if d < coreCutoff:
					numberBelowCoreCutoff += 1.0
			coreScore = numberBelowCoreCutoff/float(len(distances))

			# Condition added 2019.04.15 to properly identify exposed side chains that "appear" deeply buried...
			# A very deeply buried side chain that is distant from surface facets...
			if coreScore > 0.6:
				classifications.append((v, -1, aveDistance))
				#print("core", coreScore, (max(distances) - min(distances)), v.data)
				# if v.data.residue_sequence_number == 326:
				# 	print("--------------> core 1", coreScore, (max(distances) - min(distances)), v.data)
			# A side chain that is very exposed and distant from any surface facets...
			elif (max(distances) - min(distances)) > 5.0:
				# print(coreScore, distances)
				# print("exposed", (max(distances) - min(distances)), v.data)
				classifications.append((v, 1, aveDistance))
				# if v.data.residue_sequence_number == 326:
				# 	print("--------------> exposed 1", coreScore, (max(distances) - min(distances)), v.data)
			# A moderately buried side chain that is relatively close to the surface facets...
			else:
				classifications.append((v, -1, aveDistance))
				# print(coreScore, "core", (max(distances) - min(distances)))
				# if v.data.residue_sequence_number == 326:
				# 	print("--------------> core 2", coreScore, (max(distances) - min(distances)), v.data)

		elif aveDistance <= marginCutoff:

			if aveDistance < 0. or minDistance <= marginCutoff:
				classifications.append((v, 0, aveDistance))
			else:
				# Some exposed residues have low aveDistance values by chance.
				# However, their minDistance values are larger, enabling them to be identified.
				###############################################################################
				classifications.append((v, 1, aveDistance))

		# Margin residues that are just above the marginCutoff when using aveDistance for classification.
		# For this subset of residues, it is more appropriate to use minDistance for a consistent margin classfication.
		###############################################################################################################
		elif minDistance <= marginCutoff:

			# Relative to previous pHinder publications, this condition extends the margin classification.
			##############################################################################################
			classifications.append((v, 0, minDistance))

		else:

			classifications.append((v, 1, aveDistance))

	return classifications

# def orientS2s(triangulation):
	
# 	# @triangulation really is an s1Dict of a triangulation
# 	for node in triangulation:
# 		orientedS2s = []
# 		for s2 in triangulation[node].s2s:
# 			if node == s2[0]:
# 				orientedS2s.append(s2)
# 			else:
# 				orientedS2s.append((s2[1],s2[0]))
# 		# Update the triangulation.
# 		triangulation[node].s2s = orientedS2s

def getSurroundingEdges(closestEdge, quadSurface, been):

	# Collect the surface to identify small independent surfaces that should be removed.
	####################################################################################
	closestEdge, localEdges = quadSurface[closestEdge.id], []
	
	# First iteration over the left and right facet of the start edge.
	##################################################################
	rotate1L = rotateCwOverTheLeftFace(quadSurface[closestEdge.id], quadSurface)
	localEdges += rotate1L
	rotate1R = rotateCwOverTheRightFace(quadSurface[closestEdge.id], quadSurface)
	localEdges += rotate1R

	uniqueLocalEdges = {}
	for e in localEdges:
		if e not in been:
			uniqueLocalEdges.update({e:None})

	# A version of the quadSurface for which all small-floating-isolated surfaces have been removed.
	################################################################################################
	return uniqueLocalEdges

def goFoSurface(firstStartEdge, quadSurface, allowSmallSurfaces=0):

	# XXXX 2015.10.14 This same dictionary strategy may make normal goFo faster.
	# This approach also avoids recursion.
	############################################################################
	
	startEdgeID, findingSurfaces, surfaces, floatingEdges, unfloatingEdges = firstStartEdge, 1, [], quadSurface, {}
	while floatingEdges:
	
		# Collect the surface to identify small independent surfaces that should be removed.
		####################################################################################
		been, beenDifference = {}, 1
		edgePursuits = [floatingEdges[startEdgeID],]
		while beenDifference:

			previousBeenLength, nextEdgePursuits, possiblePursuitsEdges = len(been), {}, {}
			for e in edgePursuits:
				possiblePursuitsEdges.update(getSurroundingEdges(e, quadSurface, been))
			for e in possiblePursuitsEdges:
				been.update({e:None})
				nextEdgePursuits.update({e:None})

			edgePursuits = nextEdgePursuits
			beenDifference = len(been) - previousBeenLength

		surfaces.append(been)
		unfloatingEdges.update(been)

		# Search for other independent surfaces.
		# quadSurface keys are edges IDs, and surface keys are edge objects.
		####################################################################
		floatingEdges = {}
		for e in quadSurface:
			if quadSurface[e] not in unfloatingEdges:
				if quadSurface[e] not in been:
					floatingEdges.update({e:quadSurface[e]})
		if floatingEdges:
			startEdgeID = sorted(floatingEdges)[0]

	# n = 1
	# for surface in surfaces:
	# 	print("Surface..........................", n, len(surface))
	# 	n += 1

	# Report the number of resultant surfaces.
	##########################################
	print("\n\nSurface assessment...")
	print(80*"-")
	i, surfaceLengths = 1, {}
	for surface in surfaces:
		surfaceLengths.update({(float(len(surface)), i):surface})
		i += 1
	maxLength = max(surfaceLengths.keys())

	i, keepSurfaces = 1, []
	for surfaceLength in surfaceLengths:
		#print("Surface", i, ":", surfaceLength[0]/maxLength[0], surfaceLength[0])
		smallSurfacePass = 0
		# Facet-count floor below which a surface piece is discarded as a
		# fragment. 1000 was far too high: trypsin, at 223 residues, produces a
		# 384-facet surface and lost it entirely, which silently classified every
		# side chain as exposed. 50 is the value pHinder already uses internally
		# for hetero and void surfaces, and it still culls genuine fragments.
		minLength = 50.0
		if allowSmallSurfaces:
			minLength = allowSmallSurfaces
		if surfaceLength[0] > minLength: #surfaceLength[0]/maxLength[0] > 0.2: Changed 2019.04.16 to make removal of small isolated surface more consistent/deterministic
			# smallSurface, a, n = surfaceLengths[surfaceLength], 0., 0.
			# for smallEdge in smallSurface:
			# 	la = smallEdge.lf.area()
			# 	ra = smallEdge.rf.area()
			# 	a += la
			# 	a += ra
			# 	n += 2.0
			# averageFacetArea = a/n
			# print("Surface", i, ":", surfaceLength[0]/maxLength[0], surfaceLength[0], averageFacetArea)
		# else:
		# 	smallSurfacePass = 1
		# if smallSurfacePass:
			print("Surface", i, ":", surfaceLength[0]/maxLength[0], surfaceLength[0])
			keepSurfaces.append(surfaceLengths[surfaceLength])
			i += 1

		else:
			print("Surface", i, " skipped because too small:", surfaceLength[0], "surface facets.")
			i += 1

	finalQuadEdges = {}
	for keepSurface in keepSurfaces:
		for e in keepSurface:
			finalQuadEdges.update({e.id:quadSurface[e.id]})

	return finalQuadEdges
	
# def goFoTermini(nStart,searchNet,selfNodes,targetNodes,terminalNodes,s2s,s1Dict,searchNets):

# 	# This bit may be not be necessary.
# 	if s1Dict[nStart].s1.data.residue_name == "PSA":
# 		return 0

# 	# Record network branches.
# 	# A terminal node is in only one s2.
# 	# A pair node is in only two s2s.
# 	# A branch node is in more than two s2s.

# 	for s2 in s2s:

# 		vGo = s2[1] # Only true if the triangulation was oriented using orientS2s()

# 		# Only search OUTSIDE of the network.
# 		# Keep the search from moving through the self subnetwork.
# 		if vGo not in selfNodes:
# 			# Encountering a node in the target sub-trianglution ends the recursion.
# 			if vGo not in targetNodes:
# 				# Ensures a forward move in the network (i.e., away from the self subnetwork).
# 				if vGo not in searchNet:
# 					# Termini of the ion-network stop the recursion.
# 					if vGo not in terminalNodes:

# 						# New in v_95, 2014.05.12.
# 						##########################
# 						searchNet.append(vGo)
# 						# Recursive function call.
# 						goFoTermini(vGo,searchNet,selfNodes,targetNodes,
# 							  terminalNodes,s1Dict[vGo].s2s,s1Dict,searchNets)

# 					# Within network branch, but at terminal.
# 					else:
# 						searchNet.append(vGo)
# 				# Nothing
# 				else: 
# 					pass
# 			# Hit the target. Record the hit.
# 			else:
# 				searchNet.append(vGo)

# def consolidateSubnetworks(subNetworksByChain):

# 	###################################################			
# 	#
# 	# Iteratively consolidate degenerate sub-networks.
# 	#
# 	###################################################
# 	loopConsolidationHasOccurred = 1
# 	count = 1
# 	while loopConsolidationHasOccurred:

# 		# There are no true bigs, thus, there are no networks in this chain.
# 		####################################################################
# 		if not subNetworksByChain:
# 			loopConsolidationHasOccurred = 0
# 			continue

# 		# Loop over the master networks.
# 		################################
# 		loopConsolidationHasOccurred, consolidatedNetworks = 0, {}
# 		for n1 in subNetworksByChain:

# 			masterKeys = []
# 			masterNetwork = n1
# 			for masterKey in masterNetwork:
# 				masterKeys.append(masterKey)
# 			masterKeys.sort()

# 			# Loop over the slave networks.
# 			###############################
# 			masterWasNotConsolidated = 1
# 			for n2 in subNetworksByChain:
# 				slaveKeys = []
# 				slaveNetwork = n2
# 				for slaveKey in slaveNetwork:
# 					slaveKeys.append(slaveKey)
# 				slaveKeys.sort()

# 				# Proceed only if the master and slave keys are different
# 				#########################################################
# 				if masterKeys != slaveKeys:

# 					# Determine if the slave is within the master.
# 					consolidate = 0
# 					for slaveKey in slaveKeys:
# 						if slaveKey in masterKeys:
# 							consolidate = 1
# 							break

# 					# If not consolidated already, check if master is in slave.
# 					###########################################################
# 					if not consolidate:
# 						for masterKey in masterKeys:
# 							if masterKey in slaveKeys:
# 								consolidate = 1
# 								break					

# 					# Consolidate the two networks.
# 					#############################################################################
# 					if consolidate:

# 						# The master and slave networks are consolidated here.
# 						#####################################################################
# 						consolidatedNetworkKeys = {}
# 						for consolidatedNetworkKey in masterKeys + slaveKeys:
# 							consolidatedNetworkKeys.update({consolidatedNetworkKey:None})
# 						consolidatedNetworkKey = consolidatedNetworkKeys.keys()
# 						consolidatedNetworkKey.sort()
# 						consolidatedNetworkKey = tuple(consolidatedNetworkKey)
# 						consolidatedNetwork = consolidatedNetworkKey

# 						# Because there was a consolidation made, tell the while loop to 
# 						# allow another interation.
# 						#####################################################################
# 						loopConsolidationHasOccurred = 1
# 						masterWasNotConsolidated = 0

# 						# Update the current iteration of the network set with the newly 
# 						# consolidated network.
# 						#####################################################################
# 						consolidatedNetworks.update({consolidatedNetwork:None})

# 			# The master network is already fully consolidated, or doesn't require consolidation.
# 			#####################################################################################
# 			if masterWasNotConsolidated:
# 				consolidatedNetworks.update({tuple(masterKeys):None})

# 		# Update the list of networks with the iterative result of the consolidation procedure.
# 		#######################################################################################
# 		if loopConsolidationHasOccurred:
# 			subNetworksByChain = deepcopy(consolidatedNetworks)

# 		# No consolidations whatsoever.
# 		###############################
# 		elif count == 1:
# 			subNetworksByChain = deepcopy(consolidatedNetworks)

# 		#print("Consolidation loop number:", count )

# 		count += 1

# 	return subNetworksByChain

# def printResNums(nodeList,t):
# 	nums = ""
# 	nodeList = list(nodeList)
# 	nodeList.sort()
# 	for node in nodeList:
# 		if t[node].s1.data.residue.name == "HOH": # HOH here
# 			nums += `t[node].s1.data.residue.num` + "-" + t[node].s1.data.residue.chn + "w "
# 		else:
# 			nums += `t[node].s1.data.residue.num` + "-" + t[node].s1.data.residue.chn + " "
# 	tuple(nodeList)
# 	return nums

def writeNetworkToFile(chains, networks, triangulation, networkName, reportName, displayDirectory, outputDirectory, pdbCode, maxNetworkEdgeLength, includeRes=0):

	#import gzip

	# Write network information to file.
	####################################
	fullReportStrings = ""

	# Write network information to file by chain.
	#############################################
	reportStrings = {}
	for chain in chains:
		reportStrings.update({chain:""})

	# Sort networks by size.
	########################
	rankedNetworks = {}
	for network in networks:
		rankedNetworks.update({(len(network),network):network})
	rankedNetworksKeys = rankedNetworks.keys()
	rankedNetworksKeys.sort()
	rankedNetworksKeys.reverse()

	# Prepare each network to be written to file.
	#############################################
	j = 1
	for key in rankedNetworksKeys:
		
		network = rankedNetworks[key]

		# Make network edges for writing and display in VMD.
		####################################################
		resNumList, resString, networkString, been, networkContactOrder, contactCount = {}, "", "", [], 0., 0.
		# For each node...
		###################
		for node in network:
			# Format the node-to-node coordinates.
			######################################
			for s2 in triangulation[node].s2s:
				node1, node2 = s2[0], s2[1]
				if node1 in network and node2 in network:

					# Access the node information.
					##############################
					v1 = triangulation[s2[0]].s1
					v1Key = (v1.data.residue.chn, v1.data.residue.num, v1.data.residue.name)
					v2 = triangulation[s2[1]].s1
					networkContactOrder += abs(v1.data.residue.num - v2.data.residue.num)
					contactCount += 1.
					v2Key = (v2.data.residue.chn, v2.data.residue.num, v2.data.residue.name)
					if v1.data.residue.num not in resNumList:
						resNumList.update({v1Key:v1.data.residue.num})
					if v2.data.residue.num not in resNumList:
						resNumList.update({v2Key:v2.data.residue.num})

					# Format the node information.
					##############################
					format = "%15.6f%15.6f%15.6f "
					if s2 in been:
						continue
					been.append(s2)

					# New in v_105.
					###############
					resCoords = format % (v1.x,v1.y,v1.z) + format % (v2.x,v2.y,v2.z) + "\n"
					networkString += format % (v1.x,v1.y,v1.z) + format % (v2.x,v2.y,v2.z) + "\n"
					resStringFormat = "%s-%i-%s"
					resString += "%13s %13s | " % (resStringFormat % (v1.data.residue.chn, v1.data.residue.num, v1.data.residue.name), 
													resStringFormat % (v2.data.residue.chn, v2.data.residue.num, v2.data.residue.name))
					resString += resCoords

		# Calculate the contact order for the network.
		# Skip bridge networks and disulfide bonds.
		##############################################
		if contactCount:
			networkContactOrder = networkContactOrder/contactCount
		else:
			networkContactOrder = 0.

		# New in v_105.
		###############
		if resString:
			output = gzip.open(displayDirectory+pdbCode+"."+networkName+"Res."+str(j)+".txt.gz","wt")
			output.write(resString)
			output.close()			

		# Write network node files. v_104. 
		##################################
		if networkString:
			output = gzip.open(displayDirectory+pdbCode+"."+networkName+"."+str(j)+".txt.gz","wt")
			output.write(networkString)
			output.close()

			# Format networks to be written to file.
			#########################################
			fullResString = ""
			for vKey in resNumList:
				# Added v_104.
				##############
				if vKey[2] == "HOH":
					fullResString += str(vKey[1]) + "w-" + vKey[0] + " "
				else:
					fullResString += str(vKey[1]) + "-" + vKey[0] + " "

			# Skip bridge networks and disulfide bonds.
			###########################################
			if networkName in ["bridgeNetwork", "disulfideBonds"]:
				format = "%4s %4i %4i :  %s: \n"
			else:
				format = "%4s %4i %4i :  %s: %.0f\n"
			if networkName == "clusterNetwork":
				fullReportStrings += format % ("C",j,len(network),fullResString,networkContactOrder)
			elif networkName == "ionNetwork":
				fullReportStrings += format % ("I",j,len(network),fullResString,networkContactOrder)
			elif networkName == "tightAntiPairs":
				fullReportStrings += format % ("AP",j,len(network),fullResString,networkContactOrder)
			elif networkName == "tightIonPairs":
				fullReportStrings += format % ("IP",j,len(network),fullResString,networkContactOrder)
			elif networkName == "disulfideBonds":
				fullReportStrings += format % ("DB",j,len(network),fullResString)
			elif networkName == "bridgeNetwork":
				fullReportStrings += format % ("BN",j,len(network),fullResString)
			elif networkName == "acidWeakPointNetwork":
				fullReportStrings += format % ("AWP",j,len(network),fullResString,networkContactOrder)
			elif networkName == "weakPointNetworkHCO":
				fullReportStrings += format % ("WPHCO",j,len(network),fullResString,networkContactOrder)

			# Format networks to be written to file by chain.
			###################################################
			for chain in chains:

				# Added v_104.
				# Make and output string from the network nodes.
				################################################
				resNumListByChain = []
				for vKey in resNumList:
					if vKey[0] == chain:
						numString = ""
						if vKey[2] == "HOH":
							numString = str(vKey[1]) + "w" 
						else:
							numString = str(vKey[1])
						resNumListByChain.append(numString)
				resNumListByChain.sort()
				resString = ""
				for numString in resNumListByChain:
					resString += numString + " "

				# Update networkReportStrings: type, number, size, residue numbers.
				###################################################################
				format = "%4s %4i %4i :  %s\n"
				if networkName == "clusterNetwork":
					reportStrings[chain] += format % ("C",j,len(network),resString)
				if networkName == "ionNetwork":
					reportStrings[chain] += format % ("I",j,len(network),resString)
				if networkName == "tightAntiPairs":
					reportStrings[chain] += format % ("AP",j,len(network),resString)
				if networkName == "tightIonPairs":
					reportStrings[chain] += format % ("IP",j,len(network),resString)
				if networkName == "disulfideBonds":
					reportStrings[chain] += format % ("DB",j,len(network),resString)
				if networkName == "bridgeNetwork":
					reportStrings[chain] += format % ("BN",j,len(network),resString)
				if networkName == "acidWeakPointNetwork":
					reportStrings[chain] += format % ("AWP",j,len(network),resString)
				if networkName == "weakPointNetworkHCO":
					reportStrings[chain] += format % ("WPHCO",j,len(network),resString)
			j += 1


	# Write report strings to file.
	###############################
	if fullReportStrings:

		# Covert to Microsoft Excel format.
		####################################
		classBook = openpyxl.Workbook()
		sheet = classBook["Sheet"]
		classBook.remove(sheet)
		sheet1 = classBook.create_sheet(networkName + " Report")
		cell = sheet1.cell(row=1, column=1)
		cell.value = "Network Type"
		cell = sheet1.cell(row=1, column=2)
		cell.value = "Network Number"
		cell = sheet1.cell(row=1, column=3)
		cell.value = "Network Length"
		cell = sheet1.cell(row=1, column=4)
		cell.value = "Network Nodes"
		# book = Workbook()
		# sheet1 = book.add_sheet(networkName + " Report")
		# sheet1.write(0, 0, "Network Type")
		# sheet1.write(0, 1, "Network Number")
		# sheet1.write(0, 2, "Network Length")
		# sheet1.write(0, 3, "Network Nodes")
		if networkName not in ["bridgeNetwork", "disulfideBonds", "tightIonPairs"]:
			cell = sheet1.cell(row=1, column=5)
			cell.value = "Network Contact Order"			
			# sheet1.write(0, 4, "Network Contact Order")
		
		allStrings = fullReportStrings.split("\n")
		row = 1
		for x in allStrings:
			splitX = x.split(":")
			if len(splitX) == 1:
				continue
			fullSplit = []
			for y in splitX[0].split():
				fullSplit.append(y)
			fullSplit.append(splitX[1])
			fullSplit.append(splitX[2])
			for col in range(len(fullSplit)):
				if col in (1, 2):
					cell = sheet1.cell(row=row, column=col+1)
					cell.value = int(fullSplit[col])
					# sheet1.write(row, col, int(fullSplit[col]))
				else:
					cell = sheet1.cell(row=row, column=col+1)
					cell.value = fullSplit[col]
					# sheet1.write(row, col, fullSplit[col])
			row += 1
		if networkName == "clusterNetwork":
			book.save(outputDirectory+pdbCode+".ClusterNetworkReport.xls")
			# book.save(TemporaryFile())
		elif networkName == "ionNetwork":
			book.save(outputDirectory+pdbCode+".IonNetworkReport.xls")
			# book.save(TemporaryFile())
		elif networkName == "bridgeNetwork":
			book.save(outputDirectory+pdbCode+".BridgeNetworkReport.xls")
			# book.save(TemporaryFile())
		elif networkName == "acidWeakPointNetwork":
			book.save(outputDirectory+pdbCode+".AcidWeakPointReport.xls")
			# book.save(TemporaryFile())
		elif networkName == "weakPointNetworkHCO":
			book.save(outputDirectory+pdbCode+".WeakPointHCOReport.xls")
			# book.save(TemporaryFile())

	for chain in chains:
		if reportStrings[chain]:
			report = gzip.open(displayDirectory+pdbCode+"."+reportName+"."+chain+".txt.gz","wt")
			report.write(reportStrings[chain])
			report.close()

def calculateSurface(surfaceVertexSet, circumSphereRadiusLimit, minArea=10, highResolutionSurface=1, allowSmallSurfaces=0, writeSurfaceCreationAnimation=""):
	

	# Not enough informations (i.e. points) to calculate surface.
	##############################################################
	if len(surfaceVertexSet) < 5:
		return (None, None)
	
	# Surface calculation start time.
	#################################
	sTime = perf_counter()

	# Calculate the global convex hull and triangulation
	########################################################################
	print("\n\nTriangulating....", len(surfaceVertexSet), "vertices for pHinder surface calculation.")
	print(80*"-")
	globalHull4D = convexHull4D(surfaceVertexSet)
	print("Done triangulating....", len(surfaceVertexSet), "vertices for pHinder surface calculation. t =", (perf_counter()-sTime)/60., "min")
	
	maxID = len(surfaceVertexSet) + 1

	# Write the original triangulation to file.
	# writeSurfaceCreationAnimation must be set to the desired file output path.
	############################################################################
	if writeSurfaceCreationAnimation:
		caTriangulationS2s = globalHull4D.getS2s()
		caTriangulationTrimmed = {}
		for s2 in caTriangulationS2s:
			s2 = list(s2)
			s2.sort()
			caTriangulationTrimmed.update({tuple(s2):s2})
		caTriangulationString = ""
		caTriangulationFormat = "%15.6f%15.6f%15.6f "
		for s2 in caTriangulationTrimmed:
			v1 = globalHull4D.s1Dict[s2[0]].s1
			v2 = globalHull4D.s1Dict[s2[1]].s1
			caTriangulationString += (caTriangulationFormat % (v1.x,v1.y,v1.z) + caTriangulationFormat % (v2.x,v2.y,v2.z) + "\n")
		if not exists(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/"):
			mkdir(writeSurfaceCreationAnimation + "SurfaceCreationAnimation/")
		outFile = writeSurfaceCreationAnimation + "SurfaceCreationAnimation/" + "triangulation.CA.txt.gz"
		outputTriangulation = gzip.open(outFile,"wt")
		outputTriangulation.write(caTriangulationString)
		outputTriangulation.close()

	# Calculate the surface from .
	########################################################################
	globalSurface = surfaceGlobal(globalHull4D, circumSphereRadiusLimit, writeSurfaceCreationAnimation=writeSurfaceCreationAnimation)
	globalOrientedSurface = globalSurface[0]
	print("\n\nStart pHinder surface calculation and refinement...")
	print(80*"-")
	print("Course surface complete... t =", (perf_counter()-sTime)/60., "min")
	
	# Make edge dictionary and quad edge surface representation
	########################################################################
	globalEdgeDict = makeEdgeDict(globalOrientedSurface)
	print("MakeEdgeDict complete... t =", (perf_counter()-sTime)/60., "min")
	quadSurface = getQuadSurface(globalEdgeDict) 
	print("QuadSurface complete... t =", (perf_counter()-sTime)/60., "min")

	# Calculate the high-resolution surface.
	########################################
	if highResolutionSurface:

		print("Minimum area:", minArea)
		
		# 1) Regularize the surface by facet area...
		areaEdges, area_iteration = facetAreas(quadSurface, minArea=minArea), 1
		while areaEdges:
			print("Surface: facet refinement iteration...", area_iteration)
			for e in areaEdges:
				if e in quadSurface:
					dividePatch(quadSurface[e], maxID, quadSurface)
					maxID += 1
			areaEdges = facetAreas(quadSurface, minArea=minArea)
			area_iteration += 1
		print("Surface: facet area refinement complete... t =", (perf_counter()-sTime)/60., "min")

		
		# 2) Regularize the surface by smoothing facet divots...
		print("Surface: Smoothing surface...")
		# writeString = ""
		eKeys, already_smoothened_vertices = sorted(quadSurface), {}
		for eID in eKeys: 
			smoothen(quadSurface[eID], quadSurface, already_smoothened_vertices)
			already_smoothened_vertices[quadSurface[eID].origin_vertex] = None
			# newString = smoothen(quadSurface[eID], quadSurface)
			# if newString:
			#     writeString += newString
		# outputTesselation.write(writeString)
		# outputTesselation.close()

		print("Surface: Smoothing surface complete... t = ", (perf_counter()-sTime)/60., "min")

	# Return the surface results.
	#############################
	quadSurfaceKeys = list(quadSurface.keys())
	# There is a surface...
	if quadSurfaceKeys:
		finalQuadSurface = goFoSurface(quadSurfaceKeys[0], quadSurface, allowSmallSurfaces=allowSmallSurfaces)
	# There is no surface...
	else:
		return (None, None)

	if not finalQuadSurface:
		return (None, None)

	return (globalHull4D, finalQuadSurface)

def writeSurface_as_pymol_script(quadSurface, pdbCode, outputDirectory):

	coordinateTriples = []
	surfaceFacets = {}
	for e in quadSurface:
		surfaceFacets.update({quadSurface[e].lf.id:quadSurface[e].lf})
		surfaceFacets.update({quadSurface[e].rf.id:quadSurface[e].rf})
	surfaceFacetString = ""
	for fKey in surfaceFacets:
		f = surfaceFacets[fKey]
		t1 = f.v1.coordinate_tuple[0:-1]
		t2 = f.v2.coordinate_tuple[0:-1]
		t3 = f.v3.coordinate_tuple[0:-1]
		coordinateTriples.append((t1, t2, t3))

	script_body = "surfaceName = " + "\"surface." + pdbCode + "\"\n"
	script_body += "coordinateTriples = " + str(coordinateTriples) + "\n"
	script_body += "from pymol.cgo import *\n"
	script_body += "from pymol import cmd\n"
	script_body += "background_color = \"white\"\n"
	script_body += "cmd.bg_color(color=background_color)\n"

	script_body += "cgo_line_width = 4.0\n"
	script_body += "colorR, colorG, colorB = 0.55, 0.57, 0.67\n"
	script_body += "obj = [ BEGIN, LINES,\n"
	script_body += "		COLOR, float(colorR), float(colorG), float(colorB),]\n"
	script_body += "for coordinateTriple in coordinateTriples:\n"
	script_body += "	c = coordinateTriple\n"
	script_body += "	x1, y1, z1 = float(c[0][0]), float(c[0][1]), float(c[0][2])\n"
	script_body += "	x2, y2, z2 = float(c[1][0]), float(c[1][1]), float(c[1][2])\n"
	script_body += "	x3, y3, z3 = float(c[2][0]), float(c[2][1]), float(c[2][2])\n"
	script_body += "	obj.extend([VERTEX, x1, y1, z1,VERTEX, x2, y2, z2,])\n"
	script_body += "	obj.extend([VERTEX, x1, y1, z1,VERTEX, x3, y3, z3,])\n"
	script_body += "	obj.extend([VERTEX, x2, y2, z2,VERTEX, x3, y3, z3,])\n"
	script_body += "obj.extend([END,])\n"
	# script_body += "cmd.set('cgo_line_width', cgo_line_width)\n"
	# script_body += "cmd.load_cgo(obj, surfaceName)\n"

	script_body += "obj.extend([ BEGIN, TRIANGLES,\n"
	script_body += "		COLOR, float(colorR), float(colorG), float(colorB),])\n"
	script_body += "obj.extend([ALPHA, 0.3,])\n"
	script_body += "for coordinateTriple in coordinateTriples:\n"
	script_body += "	c = coordinateTriple\n"
	script_body += "	x1, y1, z1 = float(c[0][0]), float(c[0][1]), float(c[0][2])\n"
	script_body += "	x2, y2, z2 = float(c[1][0]), float(c[1][1]), float(c[1][2])\n"
	script_body += "	x3, y3, z3 = float(c[2][0]), float(c[2][1]), float(c[2][2])\n"
	script_body += "	obj.extend([VERTEX, x1, y1, z1,])\n"
	script_body += "	obj.extend([VERTEX, x2, y2, z2,])\n"
	script_body += "	obj.extend([VERTEX, x3, y3, z3,])\n"
	script_body += "obj.extend([END,])\n"
	script_body += "cmd.load_cgo(obj, surfaceName)\n"

	script_body += "cmd.set('ambient', 0.25)\n"
	script_body += "cmd.set('cgo_line_width', cgo_line_width)\n"

	outFile = open(outputDirectory + pdbCode + "-surface.py", "w")
	outFile.write(script_body)
	outFile.close()

	print("Done writing surface Pymol script...")

