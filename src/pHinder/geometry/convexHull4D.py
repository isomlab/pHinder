# 2013.11.05
# A major speed improvement using horizon ridge trick (starting with version 15).
# Version 16 uses this trick, but also checks for the general position of the test
# vertex against the emerging horizon ridge.
################################################################################## 

# Set the global value of zero for volume and hyperplane calculations.
# The size of zero effects how many moves it takes to jostle points into general position.
# On 2014.02.20 I changed zero from 1e-3 to 1e-5 to accommodate composite networks derived from aligned structures.
# The jostling was too extreme and added too much noise to the reduced network representation.
###################################################################################################################
# On 2015-05-05 I reduced the zero value to 1e-10 to avoid excessive attempts to achieve general position.
# Thus far, there have been no negative consequences associated with this change.
###################################################################################################################
zero = 1e-12

#print("Using version convexHull4D_1_22")
#print("Zero is set to", zero)

# Import dependencies here for enhanced performance.
####################################################
from pHinder._vendor.compGeometry import distance, Triangle, circumCircle, circumSphere, centroid3D, centroid4D, Vertex4D
from pHinder._vendor.compGeometry import geom_tol, planeCoefficients3D, planeCoefficients4D, gp2D, gp3D, gp4D
from pHinder.geometry.general_position import JostleBudget
from pHinder._vendor.determinants import det3x3, det4x4
from pHinder._vendor.pdbFile import PseudoAtom
# from convexHull3D_3_0 import convexHull3D # not used
from random import choice

############################################
# Classes for fundamental geometry objects.
############################################

class Simplex1:

	def __init__(self,v,s1Dict,bottom=0):
		
		self.s1 = v
		self.id = v.id
		# True if coincident on a bottom facet in the 4D hull.
		self.bottom = bottom

		s1Dict.update({self.id:self})
		
		# Super: 2-Simplex 
		self.s2s = []

	def copy(self):

		newCopy = Simplex1(self.s1, {})
		newCopy.bottom = self.bottom
		newCopy.s2s = self.s2s
		return newCopy

class Simplex2:

	def __init__(self,v1,v2,s1Dict,s2Dict,bottom=0):

		self.vDict = {v1.id:v1, v2.id:v2}
		self.s2 = (v1,v2)
		self.id = [v1.id,v2.id]
		self.id.sort()
		self.id = tuple(self.id)
		# True if coincident on a bottom facet in the 4D hull.
		self.bottom = bottom

		# Subs: 1-Simplex
		self.subKeys = self.vDict.keys()

		for key in self.vDict:
			va = self.vDict[key]
			if key in s1Dict:
				if key not in s1Dict[key].s2s:
					s1Dict[key].s2s.append(self.id)
			else:
				if self.bottom:
					newSimplex = Simplex1(va,s1Dict,bottom=self.bottom)
				else:
					newSimplex = Simplex1(va,s1Dict)
				s1Dict.update({key:newSimplex})   
				s1Dict[key].s2s.append(self.id)

		self.subKeys = self.vDict.keys()

		s2Dict.update({self.id:self})

		# Super: 3-simplex
		self.s3s = []
		
class Simplex3:
	
	def __init__(self,v1,v2,v3,s1Dict,s2Dict,s3Dict,bottom=0):

		self.vDict = {v1.id:v1, v2.id:v2, v3.id:v3}
		self.s3 = (v1,v2,v3)
		self.id = [v1.id,v2.id,v3.id]
		self.id.sort()
		self.id = tuple(self.id)
		# True if coincident on a bottom facet in the 4D hull.
		self.bottom = bottom

		# Subs: 2-Simplex
		s12 = [v1.id,v2.id]
		s12.sort()
		s13 = [v1.id,v3.id]
		s13.sort()
		s23 = [v2.id,v3.id]
		s23.sort()

		# Clean - up
		s2s_of_s3 = {tuple(s12):[[v1,v2],], 
						tuple(s13):[[v1,v3],], 
						tuple(s23):[[v2,v3],]}

		for key in s2s_of_s3:
			va, vb = self.vDict[key[0]], self.vDict[key[1]]
			if key in s2Dict:
				if s2s_of_s3[key][0] not in s2Dict[key].s3s:
					s2Dict[key].s3s.append(self.id)
			else:
				if self.bottom:
					newSimplex = Simplex2(va,vb,s1Dict,s2Dict,bottom=self.bottom)
				else:
					newSimplex = Simplex2(va,vb,s1Dict,s2Dict)
				s2Dict.update({key:newSimplex})
				s2Dict[key].s3s.append(self.id)

		self.subKeys = s2s_of_s3.keys()
		s3Dict.update({self.id:self})
		# Super: 3-Simplex
		self.s4s = []
		self.facet = None

	def getFacet(self):
		if self.facet:
			return self.facet
		else:
			self.facet = Triangle(self.s3[0],self.s3[1],self.s3[2],skip_orientation=1)
			return self.facet

	def __repr__(self):

		format = "%15.6f%15.6f%15.6f "
		v1, v2, v3 =  self.s3[0], self.s3[1], self.s3[2]
		string = format % (v1.x,v1.y,v1.z) + \
					format % (v2.x,v2.y,v2.z) + \
					format % (v3.x,v3.y,v3.z) + "\n"
		return string

class Simplex4:
	
	def __init__(self,v1,v2,v3,v4,s1Dict,s2Dict,s3Dict,s4Dict,bottom=0):
		
		self.vDict = {v1.id:v1, v2.id:v2, v3.id:v3, v4.id:v4}
		self.s4 = (v1,v2,v3,v4)
		self.id = [v1.id,v2.id,v3.id,v4.id]
		self.id.sort()
		self.id = tuple(self.id)
		# True if coincident on a bottom facet in the 4D hull.
		self.bottom = bottom

		# Subs: 3-Simplex
		s123 = [v1.id,v2.id,v3.id]
		s123.sort()
		s124 = [v1.id,v2.id,v4.id]
		s124.sort()
		s134 = [v1.id,v3.id,v4.id]
		s134.sort()
		s234 = [v2.id,v3.id,v4.id]
		s234.sort()

		# Clean - up
		self.s3s_of_s4 = {tuple(s123):[[v1,v2,v3,v4],],
							tuple(s124):[[v1,v2,v3,v4],],
							tuple(s134):[[v1,v2,v3,v4],],
							tuple(s234):[[v1,v2,v3,v4],]}

		for key in self.s3s_of_s4:
			va = self.vDict[key[0]]
			vb = self.vDict[key[1]]
			vc = self.vDict[key[2]]

			if key in s3Dict:
				if self.s3s_of_s4[key][0] not in s3Dict[key].s4s:
					if self.id not in s3Dict[key].s4s:
						s3Dict[key].s4s.append(self.id)
			else:
				if self.bottom:
					newSimplex = Simplex3(va,vb,vc,s1Dict,s2Dict,s3Dict,bottom=self.bottom)
				else:
					newSimplex = Simplex3(va,vb,vc,s1Dict,s2Dict,s3Dict)
				s3Dict.update({key:newSimplex})
				if self.id not in s3Dict[key].s4s:
					s3Dict[key].s4s.append(self.id)

		self.subKeys = self.s3s_of_s4.keys()
		s4Dict.update({self.id:self})
		self.circumSphere = None
		self.centroid = None

	def getCentroid(self, psa):
		if self.centroid:
			return self.centroid
		else:
			self.centroid = centroid4D(self.s4, psa)
			return self.centroid

	def getSphere(self, psa, returnSquareDistance=0):

		self.circumSphere = circumSphere(self.s4[0],self.s4[1],self.s4[2],self.s4[3],psa,
											returnSquareDistance=returnSquareDistance)
		return self.circumSphere

	def setBottom(self,value,s1Dict,s2Dict,s3Dict):

		self.bottom = value
		for a in self.subKeys:
			s3Dict[a].bottom = 1
			for b in s3Dict[a].subKeys:
				s2Dict[b].bottom = 1
				for c in s2Dict[b].subKeys:
					s1Dict[c].bottom = 1

##############################################################
# Convenient functions.
##############################################################
# def planeCoefficients3D(v1,v2,v3):

#     v1x, v1y, v1z = v1.x, v1.y, v1.z
#     v2x, v2y, v2z = v2.x, v2.y, v2.z
#     v3x, v3y, v3z = v3.x, v3.y, v3.z

#     A = det3x3([1., 1., 1., v1y, v2y, v3y, v1z, v2z, v3z])
#     B = det3x3([v1x, v2x, v3x, 1., 1., 1., v1z, v2z, v3z])
#     C = det3x3([v1x, v2x, v3x, v1y, v2y, v3y, 1., 1., 1.])
#     D = -det3x3([v1x, v2x, v3x, v1y, v2y, v3y, v1z, v2z, v3z])

#     return (A,B,C,D) 

# def planeCoefficients4D(v1,v2,v3,v4):

#     v1x, v1y, v1z, v1u = v1.x, v1.y, v1.z, v1.u
#     v2x, v2y, v2z, v2u = v2.x, v2.y, v2.z, v2.u
#     v3x, v3y, v3z, v3u = v3.x, v3.y, v3.z, v3.u
#     v4x, v4y, v4z, v4u = v4.x, v4.y, v4.z, v4.u

#     A = det4x4([1., 1., 1., 1., v1y, v2y, v3y, v4y, v1z, v2z, v3z, v4z, v1u, v2u, v3u, v4u])
#     B = det4x4([v1x, v2x, v3x, v4x, 1., 1., 1., 1., v1z, v2z, v3z, v4z, v1u, v2u, v3u, v4u])
#     C = det4x4([v1x, v2x, v3x, v4x, v1y, v2y, v3y, v4y, 1., 1., 1., 1., v1u, v2u, v3u, v4u])
#     D = det4x4([v1x, v2x, v3x, v4x, v1y, v2y, v3y, v4y, v1z, v2z, v3z, v4z, 1., 1., 1., 1.])
#     E = -det4x4([v1x, v2x, v3x, v4x, v1y, v2y, v3y, v4y, v1z, v2z, v3z, v4z, v1u, v2u, v3u, v4u])

#     return (A,B,C,D,E) 

def orientTetrahedronVerticesHP(v1,v2,v3,v4,vR):

	c  = planeCoefficients4D(v1,v2,v3,v4)
	aboveBelow = c[0]*vR.x + c[1]*vR.y + c[2]*vR.z + c[3]*vR.u + c[4]
	if abs(aboveBelow) < zero:
		print("Warning possible non-degeneracy in 4D....Jostling to remedy.")
		v1.jostle()
		v2.jostle()
		v3.jostle()
		v4.jostle()
		c  = planeCoefficients4D(v1,v2,v3,v4)
		aboveBelow = c[0]*vR.x + c[1]*vR.y + c[2]*vR.z + c[3]*vR.u + c[4]

	"""By convention, v1, v2, v3 and v4 must be arranged such that 
		vR evaluated to below the hyper-plane."""
	if aboveBelow > zero:
		
		# Swapping the first two vertices standardizes the orientation.
		temp = v1
		v1 = v2
		v2 = temp

		# Update plane coefficients.
		c  = planeCoefficients4D(v1,v2,v3,v4)

	return ((v1,v2,v3,v4),(c[0],c[1],c[2],c[3],c[4]))

def orientTetrahedronVertices(v1,v2,v3,v4,vR):

	c  = planeCoefficients4D(v1,v2,v3,v4)
	A, B, C, D, E = c[0], c[1], c[2], c[3], c[4]
	aboveBelow = A*vR.x + B*vR.y + C*vR.z + D*vR.u + E
	if abs(aboveBelow) < zero:
		print("Warning possible non-degeneracy in 4D....Jostling to remedy.")
		v1.jostle()
		v2.jostle()
		v3.jostle()
		v4.jostle()
		c  = planeCoefficients4D(v1,v2,v3,v4)
		A, B, C, D, E = c[0], c[1], c[2], c[3], c[4]
		aboveBelow = A*vR.x + B*vR.y + C*vR.z + D*vR.u + E

	"""By convention, v1, v2, v3 and v4 must be arranged such that 
		vR evaluated to below the hyper-plane."""
	if aboveBelow > zero:
		temp = v1
		v1 = v2
		v2 = temp
		
	return (v1,v2,v3,v4)

def simplex3_of_simplex4(v1,v2,v3,v4):

	vID = [v1.id, v2.id, v3.id, v4.id]
	vID.sort()
	v1id, v2id, v3id, v4id = vID[0], vID[1], vID[2], vID[3]

	return {(v1id,v2id,v3id):[[v1,v2,v3,v4],], 
			(v1id,v2id,v4id):[[v1,v2,v3,v4],], 
			(v1id,v3id,v4id):[[v1,v2,v3,v4],], 
			(v2id,v3id,v4id):[[v1,v2,v3,v4],]}

def updateHull(edgeDict,newEdges):

	for key in newEdges.keys():
		if key in edgeDict.keys():
			edgeDict[key] += newEdges[key]
		else:
			edgeDict.update({key:newEdges[key]})

def buildSimplex5(sortedVertices):

	print("Length of sorted vertices...", len(sortedVertices))

	# First, test for general position.
	#####################################

	# 3 Vertices.
	#############
	iterate = 1
	v1, v2, v3, v4, v5 = None, None, None, None, None

	# Test for general position d1.
	###############################
	v1, v2 = sortedVertices[0], sortedVertices[1]
	if distance(v1, v2) <= zero:
		budget = JostleBudget("the first two vertices of the seed simplex are coincident")
		while distance(v1, v2) <= zero:
			budget.jostle(v2)

	# Test for general position d2.
	###############################
	vGp = 0
	v3 = sortedVertices[2]
	budget = JostleBudget("the first three vertices of the seed simplex are collinear")
	while not vGp:
		vGp = gp2D(v1, v2, v3)
		if not vGp:
			budget.jostle(v3)
			continue

	# Test for general position d3.
	###############################
	vGp = 0
	v4 = sortedVertices[3]
	budget = JostleBudget("the first four vertices of the seed simplex are coplanar")
	while not vGp:
		vGp = gp3D(v1, v2, v3, v4)
		if not vGp:
			budget.jostle(v4)
			continue

	# Test for general position d4.
	###############################
	vGp = 0
	v5 = sortedVertices[4]
	budget = JostleBudget("the first five vertices of the seed simplex are cospherical")
	while not vGp:
		vGp = gp4D(v1, v2, v3, v4, v5)
		if not vGp:
			budget.jostle(v5)
			continue

	# Delete used vertices.
	#######################
	for v in [v1,v2,v3,v4,v5]:
		sortedVertices.remove(v)

	# Now that general position has been established,
	# build the initial hull in dimension 4. 
	#################################################

	# General position vertices.
	############################
	a, b, c, d, e = v1, v2, v3, v4, v5

	psa = PseudoAtom()
	vR = centroid4D([a,b,c,d,e], psa)

	t1 = orientTetrahedronVertices(a,b,c,d,vR)
	t2 = orientTetrahedronVertices(a,b,c,e,vR)
	t3 = orientTetrahedronVertices(a,b,d,e,vR)
	t4 = orientTetrahedronVertices(a,c,d,e,vR)
	t5 = orientTetrahedronVertices(b,c,d,e,vR)

	# Link facets. 
	hull4D = {}
	# t1
	hull4D.update(simplex3_of_simplex4(t1[0],t1[1],t1[2],t1[3]))
	# t2
	t2_simplex3 = simplex3_of_simplex4(t2[0],t2[1],t2[2],t2[3])
	updateHull(hull4D,t2_simplex3)
	# t3
	t3_simplex3 = simplex3_of_simplex4(t3[0],t3[1],t3[2],t3[3])
	updateHull(hull4D,t3_simplex3)
	# t4
	t4_simplex3 = simplex3_of_simplex4(t4[0],t4[1],t4[2],t4[3])
	updateHull(hull4D,t4_simplex3)
	# t5
	t5_simplex3 = simplex3_of_simplex4(t5[0],t5[1],t5[2],t5[3])
	updateHull(hull4D,t5_simplex3)

	# Check orientation
	newEdgeDict = {}
	for e in hull4D:

		f1 = hull4D[e][0]
		f2 = hull4D[e][1]
		vT1, vs1 = None, []
		for vID in e:
			for v in (f1[0],f1[1],f1[2],f1[3]):
				if vID == v.id:
					vs1.append(v)
				elif v.id not in e:
					vT1 = v
		vT2, vs2 = None, []
		for vID in e:
			for v in (f2[0],f2[1],f2[2],f2[3]):
				if vID == v.id:
					vs2.append(v)
				elif v.id not in e:
					vT2 = v

		# Now all edges and their tetrahera are oriented consistently
		v1, v2, v3, v4 = vT1, vs1[0], vs1[1], vs2[2]
		nF = orientTetrahedronVertices(v1,v2,v3,v4,vT2)
		nF_simplex3 = simplex3_of_simplex4(nF[0],nF[1],nF[2],nF[3])
		if e not in newEdgeDict:
			newEdgeDict.update({e:nF_simplex3[e]})
		else:
			newEdgeDict[e] += nF_simplex3[e]

		v1, v2, v3, v4 = vT2, vs2[0], vs2[1], vs2[2]
		nF = orientTetrahedronVertices(v1,v2,v3,v4,vT1)
		nF_simplex3 = simplex3_of_simplex4(nF[0],nF[1],nF[2],nF[3])
		if e not in newEdgeDict:
			newEdgeDict.update({e:nF_simplex3[e]})
		else:
			newEdgeDict[e] += nF_simplex3[e]

	hull4D = newEdgeDict

	return hull4D # This is a 5-simplex.

def getSigns(vT, t1, t2, hyperPlanes):

	# Now fully consistent, i.e. identicle, with gp4D

	# Get the tetrahedron incident on the startEdge.
	################################################
	c = hyperPlanes[tuple(t1)]
	v1 = (c[0]*vT.x + c[1]*vT.y + c[2]*vT.z + c[3]*vT.u + c[4])/24.
	c = hyperPlanes[tuple(t2)]
	v2 = (c[0]*vT.x + c[1]*vT.y + c[2]*vT.z + c[3]*vT.u + c[4])/24.

	if abs(v1) <= zero or abs(v2) <= zero:
		# Not in general position.
		return 0
	else:
		# Convert 4D hyper-volumes to binary score.
		return (v1/abs(v1), v2/abs(v2))

	# Convert 4D hyper-volumes to binary score.
	###########################################
	return (v1/aV1, v2/aV2)

def getEmergingEdges(startEdge, hull4D):

	# This version is a bit faster because it uses dictionary look-ups
	# which are faster than list look-ups.  
	##################################################################

	step = 1
	nextEdges, iEdges, hullEdges  = [startEdge,], [startEdge,], {}
	history = {}
	while nextEdges:

		nextEdges = {}
		for e in iEdges:

			step += 1

			# Update hull edges.
			####################
			hullEdges.update({e:hull4D[e]})        

			# Get the tetrahedron incident on the startEdge.
			################################################
			t1, t2 = hull4D[e][0], hull4D[e][1]

			# Get the triangle facet "edges" of tetrahedral facet 1.
			########################################################
			edgeKeys = simplex3_of_simplex4(t1[0], t1[1], t1[2], t1[3])

			# Get the triangle facet "edges" of tetrahedral facet 2.
			########################################################
			edgeKeys.update(simplex3_of_simplex4(t2[0], t2[1], t2[2], t2[3]))

			# Targeted, recursive search from horizon edges.
			################################################
			for pursuantEdge in edgeKeys:
				if pursuantEdge in history:
					continue
				else:
					history.update({pursuantEdge:None}) #may not be correct
					nextEdges.update({pursuantEdge:None})
		
		nextEdges = nextEdges.keys()

		iEdges = nextEdges
		
	return hullEdges

def findHorizonEdge(vT, startEdge, hull4D, hyperPlanes):

	nextEdges, iEdges, horizonEdge  = [startEdge,], [startEdge,], None
	history = {}
	while nextEdges:

		nextEdges = {}
		for e in iEdges:

			# Get the facets of the edge, e.
			################################
			t1, t2 = hull4D[e][0], hull4D[e][1]

			signs = getSigns(vT, t1, t2, hyperPlanes)
			# vT is not in general position relative to iEdge.
			##################################################
			if not signs:
				return -1
			# vT is in general position relative to iEdge.
			##################################################
			v1, v2 = signs[0], signs[1]

			# STOP at horizon edge.
			#######################
			if v1 != v2:
				horizonEdge = e
				nextEdges = {}
				break

			# Get the triangle facet "edges" of tetrahedral facet 1.
			########################################################
			edgeKeys = simplex3_of_simplex4(t1[0], t1[1], t1[2], t1[3])

			# Get the triangle facet "edges" of tetrahedral facet 2.
			########################################################
			edgeKeys.update(simplex3_of_simplex4(t2[0], t2[1], t2[2], t2[3]))

			# Targeted, recursive search from horizon edges.
			################################################
			for pursuantEdge in edgeKeys:
				if pursuantEdge in history:
					continue
				else:
					history.update({pursuantEdge:None}) #may not be correct
					nextEdges.update({pursuantEdge:None})
		
		nextEdges = nextEdges.keys()

		iEdges = nextEdges

	return horizonEdge

def getHorizonRidge(vT, startEdge, hull4D, hyperPlanes, tiers=2):
	
	# !!!! This functions assumes that startEdge is a horizon edge.
	# Have to set startEdge using findHorizonEdge prior to calling this function.
	#############################################################################
	
	nextEdges, iEdges, horizonEdges  = [startEdge,], [startEdge,], {}
	history = {} # testing
	while nextEdges:

		nextEdges = {}
		for e in iEdges:

			#history.update({e:hull4D[e]}) # added

			if e not in horizonEdges and e not in history:

				# Get the facets of the edge, e.
				################################
				t1, t2 = hull4D[e][0], hull4D[e][1]

				signs = getSigns(vT, t1, t2, hyperPlanes)
				if not signs:
					return -1
				v1, v2 = signs[0], signs[1]

				# Condition for horizon edge.
				##############################
				if v1 != v2:

					t_duplicates = {tuple(t1):None, tuple(t2):None}

					horizonEdges.update({e:hull4D[e]})
					history.update({e:hull4D[e]}) # added

					allTierEdges = {}

					# Always need tier 1.
					#####################
					tier1edges = simplex3_of_simplex4(t1[0], t1[1], t1[2], t1[3])
					tier1edges.update(simplex3_of_simplex4(t2[0], t2[1], t2[2], t2[3]))
					allTierEdges.update(tier1edges)

					# Always need tier 2.
					#####################
					tier2edges = {}
					for pairEdges in tier1edges:
						t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
						if tuple(t12) not in t_duplicates:
							tier2edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							t_duplicates[tuple(t12)] = None
						if tuple(t22) not in t_duplicates:
							tier2edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3])) 
							t_duplicates[tuple(t22)] = None
					allTierEdges.update(tier2edges)

					# Add tier 3.
					#############
					if tiers == 3:

						# In some cases need tier 3.
						############################
						tier3edges = {}
						for pairEdges in tier2edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if tuple(t12) not in t_duplicates:
								tier3edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
								t_duplicates[tuple(t12)] = None
							if tuple(t22) not in t_duplicates:
								tier3edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3])) 
								t_duplicates[tuple(t22)] = None
						allTierEdges.update(tier3edges)

					# Add tier 4.
					if tiers == 4:

						# In some cases need tier 4.
						tier3edges = {}
						for pairEdges in tier2edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if tuple(t12) not in t_duplicates:
								tier3edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
								t_duplicates[tuple(t12)] = None
							if tuple(t22) not in t_duplicates:
								tier3edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
								t_duplicates[tuple(t22)] = None 
						allTierEdges.update(tier3edges)

						# In some cases need tier 4.
						tier4edges = {}
						for pairEdges in tier3edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if tuple(t12) not in t_duplicates:
								tier4edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
								t_duplicates[tuple(t12)] = None
							if tuple(t22) not in t_duplicates:
								tier4edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
								t_duplicates[tuple(t22)] = None 
						allTierEdges.update(tier4edges)
				
					if tiers == 5:
						
						# In some cases need tier 5.
						tier3edges = {}
						for pairEdges in tier2edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier3edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier3edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier3edges)
						
						# In some cases need tier 5.
						tier4edges = {}
						for pairEdges in tier3edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier4edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier4edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier4edges)
			
						# In some cases need tier 5.
						tier5edges = {}
						for pairEdges in tier4edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier5edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier5edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier5edges)
	
					if tiers == 6:

						# In some cases need tier 6.
						tier3edges = {}
						for pairEdges in tier2edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier3edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier3edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier3edges)
						
						# In some cases need tier 6.
						tier4edges = {}
						for pairEdges in tier3edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier4edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier4edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier4edges)
						
						# In some cases need tier 6.
						tier5edges = {}
						for pairEdges in tier4edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier5edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier5edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier5edges)
							
						# In some cases need tier 6.
						tier6edges = {}
						for pairEdges in tier5edges:
							t12, t22 = hull4D[pairEdges][0], hull4D[pairEdges][1]
							if t12 not in [t1, t2]:
								tier6edges.update(simplex3_of_simplex4(t12[0], t12[1], t12[2], t12[3]))
							if t22 not in [t1, t2]:
								tier6edges.update(simplex3_of_simplex4(t22[0], t22[1], t22[2], t22[3]))
						allTierEdges.update(tier6edges)

					# Targeted, recursive search from horizon edges.
					################################################
					for pursuantEdge in allTierEdges:
						# try:
						# 	history[pursuantEdge]
						# 	continue
						# except:
						# 	history.update({pursuantEdge:None})
						nextEdges.update({pursuantEdge:None})

				else:
					history.update({e:hull4D[e]})
		
		nextEdges = nextEdges.keys()

		iEdges = nextEdges

	return horizonEdges

##############################################################
# Class for the convex hull of a set of points in dimension-4.
##############################################################

class convexHull4D:

	def __init__(self, vList, startID=1, skipSort=0, walkDistanceLimit=0.5, gpFailTimesLimit=10, skip_jostle=1):

		"""Argument description.
			vList -- a list of 4D vertex objects.
					vertex4D is located in the compGeometry module.
		"""
		
		vList, self.empty_hull = list(vList), False
		if not vList:

			print("The vertex list is empty...exiting calculation.")
			self.empty_hull = True

		else:

			# If vList does not contain at least 5 vertices, 
			# it cannot be triangulated without adding dummy vertices.
			##########################################################
			i = 1
			if len(vList) < 5:

				# Collect the possible chain names for the necessary PseudoAtom objects.
				########################################################################
				vChains = {}
				for v in vList:
					vChains.update({v.data.residue.chn:None})

				# Generate possible PseudoAtom objects.
				#######################################
				templateTetrahedron = [(1.,0.,0.), (-1.,0.,0.), (0.,1.,0), (0.,0.,1.)]
				for c in templateTetrahedron:
					x, y, z, u = c[0], c[1], c[2], c[0]**2 + c[1]**2 + c[2]**2
					psa = PseudoAtom()
					psa.x, psa.y, psa.z, psa.atom_serial, psa.residue_sequence_number, psa.residue_name = x, y, z, i, i, "PSA"
					psa.chain_identifier = list(vChains.keys())[0]
					psa.charge = 0
					psa.reinitialize()
					v = Vertex4D((x, y, z, u))
					v.data = psa
					vList.append(v)
					i += 1

			# 2022.03.09
			# It would be best of the point set was determined to be in general position here...
			# Rather, than as the hull is built (which would ultimately be most efficient)...
			# However, that is an expensive calculation...
			# The easiest/fastest solution would be to jostle each point...
			# This does not formally guarantee GP, but may come close heuristically...

			# from math import floor
			# gp_vList, gpTestingDict2, gpTestingDict3, gpTestingDict4 = [], {}, {}, {}
			# for v in vList:
			# 	key1 = (floor(v.x), floor(v.y))
			# 	key2 = (floor(v.x), floor(v.y), floor(v.z))
			# 	key3 = (floor(v.x), floor(v.y), floor(v.z), floor(v.u))
			# 	if key1 not in gpTestingDict3:
			# 		gpTestingDict2[key1] = None
			# 		if key2 not in gpTestingDict3:
			# 			gpTestingDict3[key2] = None
			# 			if key3 not in gpTestingDict4:
			# 				# gpTestingDict2[key1] = None
			# 				# gpTestingDict3[key2] = None
			# 				gpTestingDict3[key3] = None
			# 				gp_vList.append(v)
			# print("B", len(vList))
			# vList = gp_vList

			# 2022.02.10 This loop heuristically seems to achieve general position in simple, elegant, and rapid way.
			# Notes:
			# If the point set is known to be in general position, jostling can be skipped for a deterministic result.
			# If the point set is jostled, the result is not deterministic, meaning the results of the calculation...
			# ...may vary slightly between sequential runs. This usually has no impact on the expected result.
			# Lastly, if it is known or expected that the point set may not be in general position, it must be jostled.
			if not skip_jostle:
				for v in vList:
					v.jostle(perturbLevel=0)

			# The primary triangulaton attributes.
			######################################
			self.triangulation, self.bottomHull = {}, {}
			self.s1Dict, self.s2Dict, self.s3Dict, self.s4Dict = {}, {}, {}, {}

			# Sort the vertices by their u-coordinate
			##########################################

			# Sort the vertices ensuring a unique 3D vertex list.
			#####################################################
			self.vDict, sortedVertices = {}, []
			for v in vList:
				key = (v.u, v.x, v.y, v.z)
				# if key in self.vDict:
				# 	continue
				self.vDict.update({key:v})
			keys = sorted(self.vDict)

			# Re-number the vertices.
			#########################
			new_vDict = {}
			for key in keys:

				# Vertex4D object.
				##################
				v = self.vDict[key]
				v.id = startID
				# Match Vertex4D objects.
				#########################
				v.reinitialize()
				v.data.v = v
				new_vDict[v.id] = v
				sortedVertices.append(v)
				startID += 1
			self.vDict = new_vDict
			sortedVertices = list(self.vDict.values())

			# Calculate zero based on coordinate scale
			global zero
			zero = geom_tol(sortedVertices)

			# Begin calculating the convex hull in 4D.
			##########################################
			self.hull4D = buildSimplex5(sortedVertices)

			# Create the initial hyper-planes.
			# Hyperplane dictionary lookups saves much time, as apposed to 
			# recalculating hyperplanes for every run of the loop.
			###############################################################
			self.hyperPlanes = {}
			for e in self.hull4D:
				f1 = self.hull4D[e][0]
				v1, v2, v3, v4 = f1[0], f1[1], f1[2], f1[3]
				c  = planeCoefficients4D(v1, v2, v3, v4)
				A, B, C, D, E = c[0], c[1], c[2], c[3], c[4]
				self.hyperPlanes.update({(v1,v2,v3,v4):(A,B,C,D,E)})

				f2 = self.hull4D[e][1]
				v1, v2, v3, v4 = f2[0], f2[1], f2[2], f2[3]
				c  = planeCoefficients4D(v1, v2, v3, v4)
				A, B, C, D, E = c[0], c[1], c[2], c[3], c[4]
				self.hyperPlanes.update({(v1,v2,v3,v4):(A,B,C,D,E)})
			
			print("\nconvexHull4D_1_22: convexHull4D.__init__(): Triangulating", len(sortedVertices) + 5, "vertices...")
			printFrequency = 500
			if len(sortedVertices) <= 1000:
				printFrequency = 100

			# Calculate the hull in 4D.
			###########################
			startEdge = sorted(self.hull4D)[0] # start edge for traversing horizon ridge.
			i, searchList = 0, [startEdge,]
			while i < len(sortedVertices):
				
				if (i + 6) % printFrequency == 0.:
					print("convexHull4D_1_22: convexHull4D.__init__(): Triangulating Vertex", i + 5, "of", len(sortedVertices) + 5)

				# Get the next vertex in the set.
				#################################
				vT = sortedVertices[i]

				startHorizonEdge = findHorizonEdge(vT, startEdge, self.hull4D, self.hyperPlanes)
				if startHorizonEdge == -1:
					i += 1
					continue

				newEdgeDict, tiers = -1, 1
				budget = JostleBudget("the horizon ridge could not be closed for vertex %d" % i)
				while newEdgeDict == -1:
					horizonEdges = getHorizonRidge(vT, startHorizonEdge, self.hull4D, self.hyperPlanes, tiers=tiers)
					if horizonEdges == -1:
						budget.jostle(vT)
						# print(vT.x, vT.y, vT.z, tiers)
						continue
					newEdgeDict = self.buildHorizonCone(vT, horizonEdges)
					tiers += 1

				# Set the next start edge to keep the iterative search for
				# the horizon ridge away from old hull edges.
				##########################################################
				keys = sorted(newEdgeDict)
				keys.reverse()
				# If startEdge is empty, vT was inside of the current convex hull, so continue...
				if not keys:
					i += 1
					continue
				startEdge = keys[0]

				# This should fully update the hull.
				# It should also automatically unlink old hull edges.
				#####################################################
				self.hull4D.update(newEdgeDict)

				# Increment to next vertex index.
				#################################
				i += 1

			# Set the final hull representation.
			####################################
			finalHull = getEmergingEdges(startEdge, self.hull4D)
			self.hull4D = finalHull
			
			# Orient the s2 edges of each s1Dict.
			#####################################
			self.orientS2s()

			# Identify the bottom facets of the hull.
			#########################################
			self.hasseDiagram()
			
			# Derive the triangulation from the hull. 
			#########################################
			self.makeFinalEdgeRepresentation()

	def buildHorizonCone(self, vT, horizonEdges):

		filteredHorizonEdges = {}
		for e in horizonEdges:

			t1, t2 = horizonEdges[e][0], horizonEdges[e][1]
			
			signs = getSigns(vT, t1, t2, self.hyperPlanes)
			if not signs:
				return -1
			v1, v2 = signs[0], signs[1]

			# print("signs......", v1, v2)

			# Determine which facet of the horizon edge to keep...
			if v1 > 0.:
				filteredHorizonEdges.update({e:[[],t2]})
			else:
				filteredHorizonEdges.update({e:[t1,[]]})

		# Update with corrected version of horizon edges.
		#################################################
		horizonEdges = filteredHorizonEdges

		# Make the new facets incident on the horizon edges.
		####################################################
		newEdgeDict, newHyperPlanes = {}, {}
		for e in horizonEdges:

			# Identify the facet that needs to be built...
			# Preserving this index maintains proper orientation...
			updateIndex, dataIndex = horizonEdges[e].index([]), 0
			if updateIndex == 0:
				dataIndex = 1
			else:
				dataIndex = 0

			# Preserving indices to replace with new vertices...
			le = list(e)
			# print("------------------------------------------------>", type(le), vT)
			f = horizonEdges[e][dataIndex]
			fIDs = [f[0].id,f[1].id,f[2].id,f[3].id]
			ppc, refV = [0,0,0], None
			for fID in fIDs:

				if fID in le:
					lei = le.index(fID)
					ppc[lei] = fIDs.index(fID)
				else:
					refVi = fIDs.index(fID)
					refV = f[refVi]

			# 2022.03.10
			# There has been a intermittent bug in some surface calculation...
			# The bug has revealed itself over the years on particularly challengings surfaces...
			# such as those calculated from larger point sets and lower circumsphere cutoff values...
			# This bug has likely persisted since the inception of the software many years ago...
			# It was causing major problems in the giga-docking algorithms that calculate surfaces of ligands and sampling voids...
			# It appears I have tracked the bug to here and fixed it...
			# I believe the problem was that the reference vertex when building new convexhull facets in 4D was inconsistently set...
			# This lead to intermittent incorrect assignments in which the reference vertex was contained within the edge of it's neighboring 3D facet...
			# The conditional statement that follows ensures this cannot happen...
			# This seems to have fixed the code and leads to beautiful and consistent surfaces, independent of input parameters...
			ppc0, ppc1, ppc2 = ppc[0], ppc[1], ppc[2]
			v1, v2, v3 = self.vDict[fIDs[ppc0]], self.vDict[fIDs[ppc1]], self.vDict[fIDs[ppc2]]
			vTest = vT # XXXX 3/9/22
			if vT.id in (v1.id, v2.id, v3.id): # XXXX 3/9/22
				vTest = refV

			# ppc0, ppc1, ppc2 = ppc[0], ppc[1], ppc[2]
			#v1, v2, v3, v4 = f[ppc0],f[ppc1],f[ppc2],vT
			v1, v2, v3, v4 = self.vDict[fIDs[ppc0]], self.vDict[fIDs[ppc1]], self.vDict[fIDs[ppc2]], vTest # Related to v vs v.id # XXXX 3/9/22

			# testSet = set([v1.id, v2.id, v3.id, v4.id])
			# if len(testSet) < 4:
			# 	print("test set...", testSet, [v1.id, v2.id, v3.id, v4.id])
			# 	raise SystemExit

			# This code builds the new horizon cone, but does not fully link it.
			####################################################################
			newF = [0,0,0,0]
			if abs(ppc0 - ppc1) == 1:
				# The swap.
				newF[ppc0] = v2
				newF[ppc1] = v1
				# End swap.
				newF[ppc2] = v3
				# print(1, newF)
				# print(1, ppc0, ppc1, ppc2)
				lastIndex = newF.index(0)
				# print("lastIndex", lastIndex)
				newF[lastIndex] = vTest # XXXX 3/9/22
				v1, v2, v3, v4 = newF[0], newF[1], newF[2], newF[3]
				# # XXXX 3/9/22
				# testAdd = set([v1.id, v2.id, v3.id, v4.id]) # difference between using vertex and vertex id...problem identified with id
				# if len(testAdd) < 4:
				# 	print("Problem.......1", [v1.id, v2.id, v3.id, v4.id], testAdd)
				# 	raise SystemExit
				# # print(1, v4, newF, lastIndex)
			else:
				newF[ppc0] = v1
				# The swap. 
				newF[ppc1] = v3
				newF[ppc2] = v2
				# End swap.
				# print(2, newF)
				# print(2, ppc0, ppc1, ppc2)
				lastIndex = newF.index(0)
				# print("lastIndex", lastIndex)
				newF[lastIndex] = vTest # XXXX 3/9/22 	
				v1, v2, v3, v4 = newF[0], newF[1], newF[2], newF[3]
				# testAdd = set([v1.id, v2.id, v3.id, v4.id])
				# if len(testAdd) < 4:
				# 	print("Problem.......2", [v1.id, v2.id, v3.id, v4.id], testAdd)
				# 	raise SystemExit
				# # print(2, v4, vT, newF)

			# New
			c = planeCoefficients4D(v1, v2, v3, v4)
			newHyperPlanes.update({(v1,v2,v3,v4):c})
			# End New

			# Pair the new cone tetraderon through their incident triangle edges.
			#####################################################################
			nF_simplex3 = simplex3_of_simplex4(v1,v2,v3,v4)
			for key in nF_simplex3.keys():
				if key not in newEdgeDict:
					newEdgeDict.update({key:nF_simplex3[key]})
				else:
					newEdgeDict[key] += nF_simplex3[key]

			# Update the horizon edge with the newly built facet...
			horizonEdges[e][updateIndex] = nF_simplex3[e][0]

			# Update the hull4D, eventually...
			##################################
			# a = horizonEdges[e]
			# print("******************", a)
			# testSet1 = set([a[0][0].id, a[0][1].id, a[0][2].id, a[0][3].id])
			# print(testSet1)
			# testSet2 = set([a[1][0].id, a[1][1].id, a[1][2].id, a[1][3].id])
			# print(testSet2)
			newEdgeDict.update({e:horizonEdges[e]})
	
		# Check for linkage mistakes.
		#############################
		linkageMistake = 0
		for e in newEdgeDict:
			if len(newEdgeDict[e]) == 1:
				linkageMistake = 1
				return -1

		self.hyperPlanes.update(newHyperPlanes)
		return newEdgeDict

	def orientS2s(self):
		
		# @triangulation really is an s1Dict of a triangulation
		for node in self.s1Dict:
			orientedS2s = []
			for s2 in self.s1Dict[node].s2s:
				if node == s2[0]:
					orientedS2s.append(s2)
				else:
					orientedS2s.append((s2[1],s2[0]))
			# Update the triangulation.
			self.s1Dict[node].s2s = orientedS2s

	def hasseDiagram(self):

		self.bottomHull = {}
		self.s1Dict, self.s2Dict, self.s3Dict, self.s4Dict = {}, {}, {}, {}
		
		"""This function also identifies the bottom facets of the 4d-hull.
			It is these facets whose coordinates comprise the Delaunay 
			triangulation when they are projected down to dimension-3."""
	
		# The proper reference for the 1-Simplex is critical for identifying the bottom of the hull in dimension 4.
		###########################################################################################################
		vR = Vertex4D((0., 0., 0., (-1000**2 + -1000**2 + -1000**2)))
		for e in self.hull4D:

			f1 = self.hull4D[e][0]
			s41 = Simplex4(f1[0],f1[1],f1[2],f1[3],
							self.s1Dict,self.s2Dict,self.s3Dict,self.s4Dict)
			c  = planeCoefficients4D(f1[0], f1[1], f1[2], f1[3])
			A, B, C, D, E = c[0], c[1], c[2], c[3], c[4]
			f1AboveBelow = A*vR.x + B*vR.y + C*vR.z + D*vR.u + E
			if f1AboveBelow > 0:
				s41.setBottom(1,self.s1Dict,self.s2Dict,self.s3Dict)
				if s41.id not in self.bottomHull:
					self.bottomHull.update({s41.id:s41})

			f2 = self.hull4D[e][1]
			s42 = Simplex4(f2[0],f2[1],f2[2],f2[3],
							self.s1Dict,self.s2Dict,self.s3Dict,self.s4Dict)
			c  = planeCoefficients4D(f2[0], f2[1], f2[2], f2[3])
			A, B, C, D, E = c[0], c[1], c[2], c[3], c[4]
			f2AboveBelow = A*vR.x + B*vR.y + C*vR.z + D*vR.u + E
			if f2AboveBelow > 0:
				s42.setBottom(1,self.s1Dict,self.s2Dict,self.s3Dict)
				if s42.id not in self.bottomHull:
					self.bottomHull.update({s42.id:s42})

	def makeFinalEdgeRepresentation(self):
		
		# This function returns the all of the edges (3-simplices) associated
		# with facets (4-simplices) that are visible (by hyperplane calculation
		# in dimension 4) from the bottom of the hull. This subset of edges
		# comprise the data structure that is the Delaunay triangulation in dimension 3.
		################################################################################
		self.triangulation = {}
		edges = {}
		for s4 in self.bottomHull.values():
			for s3 in s4.subKeys:
				if not self.s3Dict[s3].bottom:
					continue
				if s3 not in edges:
					edges.update({s3:[s4]})
				else:
					edges[s3].append(s4)
		
		self.triangulation = edges

	def stringReprSimplices2(self, cutOff=None, cutOffLow=None):

		# 2013.10.13
		# When visualizing the pHinder Megalo triangulation I discovered
		# that the original version of this function was not working
		# properly.  Although the self.triangulation object is correct, 
		# my approach for displaying the edges of the triangulation were
		# flawed.  The result was that some s2 edges of the triangulation were
		# not being appended to the string object intend to contain the coordinates 
		# of the triangulation s2 edges.  To fix this probem, I now use the s3 edges 
		# of the self.triangulation object to access 1) each s2 object, and 2)
		# the s1 object that comprise each s2 oject.  Because the self.triangulation
		# represents the subset of s3 object visible from the "bottom" of the hull, 
		# there is no longer a need to check the bottom status of the s2 object, 
		# as was done in the original function, and was likely the source of the issue.
		###############################################################################

		s2Dict = {} # Dictionary of the s2 object of the triangulation.
		string = ""
		format = "%15.6f%15.6f%15.6f "
		for s3Key in self.triangulation:
			for s2Key in self.s3Dict[s3Key].subKeys:
				s2 = self.s2Dict[s2Key]
				if s2.id not in s2Dict and (s2.id[1], s2.id[0]) not in s2Dict:
					if cutOff:
						if cutOffLow and (cutOffLow < cutOff):
							v1, v2 = s2.s2[0], s2.s2[1]
							d = distance(v1, v2, returnSquareDistance=1)
							if d < cutOff**2 and d > cutOffLow**2:
								string += (format % (v1.x,v1.y,v1.z) + 
											format % (v2.x,v2.y,v2.z) + "\n")  
								s2Dict.update({s2.id:s2})
						else:
							v1, v2 = s2.s2[0], s2.s2[1]
							if distance(v1, v2, returnSquareDistance=1) < cutOff**2:
								string += (format % (v1.x,v1.y,v1.z) + 
											format % (v2.x,v2.y,v2.z) + "\n")  
								s2Dict.update({s2.id:s2})
					else:
						string += (format % (v1.x,v1.y,v1.z) + 
									format % (v2.x,v2.y,v2.z) + "\n") 
						s2Dict.update({s2.id:s2})
		return string

	def getS2s(self, cutOff=100000000, cutOffLow=0):

		s2Dict = {} # Dictionary of the s2 object of the triangulation.
		for s3Key in self.triangulation:
			for s2Key in self.s3Dict[s3Key].subKeys:
				s2Instance = self.s2Dict[s2Key]
				v1, v2 = self.s1Dict[s2Key[0]].s1, self.s1Dict[s2Key[1]].s1
				if distance(v1,v2) >= cutOffLow and distance(v1,v2) <= cutOff:
					s2Dict.update({s2Instance.id:s2Instance})
		return s2Dict

	def s1gets2(self, s1Key, cutOff=None, cutOffLow=None):

		# Get the s2 objects incident on the s1 object represented by s1Key.
		####################################################################
		s2List = []
		for s2Key in self.s1Dict[s1Key].s2s: 
			simplex2 = self.s2Dict[s2Key]
			if simplex2.bottom:
				v1, v2 =  simplex2.s2[0], simplex2.s2[1]
				d = distance(v1, v2, returnSquareDistance=1)
				if cutOff:
					if cutOffLow and (cutOffLow < cutOff):
						if d < cutOff**2 and d > cutOffLow**2:
							s2List.append(s2Key)
					else:
						if d < cutOff**2:
							s2List.append(s2Key)
				else:
					s2List.append(s2Key)
		return s2List

	def s1gets4(self,s1Key,stringRepr=0,cutOff=None):

			# 1-Simplex
			simplex1 = self.s1Dict[s1Key]

			# 2-Simplices
			s2sList = []
			for s2Key in self.s1Dict[s1Key].s2s:
				s2sList.append(self.s2Dict[s2Key]) 

			# 3-Simplices
			s3sList = []
			for s2 in s2sList:
				for s3Key in s2.s3s:
					s3sList.append(self.s3Dict[s3Key])

			# 4-Simplices
			s4sDict, s4Keep = {}, {}
			for s3 in s3sList:
				for s4Key in s3.s4s:
					s4sDict.update({s4Key:self.s4Dict[s4Key]})
				for s4Key in s4sDict:
					skip = 0
					for v1 in s4sDict[s4Key].s4:
						for v2 in  s4sDict[s4Key].s4:
							if v1 != v2:
								d = distance(v1.data, v2.data, returnSquareDistance=1 )
								if d > cutOff**2:
									skip = 1
									break
						if skip == 1:
							break
					if not skip:
						s4Keep.update({s4Key:s4sDict[s4Key]})
					if skip:
						print("skipping")

			return s4Keep
	
	def s4toKeep(self,cutOff):
			
		s4sDict, s4Keep = {}, {}
		for s4Key in self.s4Dict:
			skip = 0
			for v1 in self.s4Dict[s4Key].s4:
				for v2 in  self.s4Dict[s4Key].s4:
					if v1 != v2:
						d = distance(v1.data, v2.data, returnSquareDistance=1)
						if d > cutOff**2:
							skip = 1
							break
				if skip == 1:
					break
			if not skip:
				s4Keep.update({s4Key:self.s4Dict[s4Key]})
			if skip:
				print("skipping")

		return s4Keep
