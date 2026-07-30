zero = 1e-12

from pHinder._vendor.compGeometry import geom_tol, Vertex, Edge, Triangle, gp2D, gp3D, planeCoefficients3D, distance
import pHinder._vendor.compGeometry as cg

class triangulation2D:

	# This class is used to extract and performed operations on 
	# a 2D triangulation calculated using a 3D convex hull.

	def __init__(self, bottom_edges, convex_hull):

		self.triangulation = bottom_edges
		self.convex_hull = convex_hull
		self.all_edges = {}
		self.edges_same_labels = {}
		self.edges_different_labels = {}

	def get_all_edges(self, d_cutoff=0):

		self.all_edges = {}
		for e in self.triangulation:
			e = self.triangulation[e]
			if d_cutoff:
				import math
				x1, y1 = e.origin_vertex.x, e.origin_vertex.y
				x2, y2 = e.destination_vertex.x, e.destination_vertex.y
				d = math.sqrt((x1-x2)**2 + (y1-y2)**2)
				if d < d_cutoff:
					self.all_edges[e.id] = e 
			else:
				self.all_edges[e.id] = e 

	def get_edges_by_labels(self, labels=[], d_cutoff=0):
		
		if not self.all_edges:
			self.get_all_edges(d_cutoff=d_cutoff)

		for e in self.all_edges:
			e = self.all_edges[e]
			if e.bottom:
				origin = e.origin_vertex
				destination = e.destination_vertex
				origin_label = e.origin_vertex.data.label
				destination_label = e.destination_vertex.data.label
				if origin_label in labels and destination_label in labels:
					if origin_label != destination_label:
						if d_cutoff:
							import math
							x1, y1 = e.origin_vertex.x, e.origin_vertex.y
							x2, y2 = e.destination_vertex.x, e.destination_vertex.y
							d = math.sqrt((x1-x2)**2 + (y1-y2)**2)
							if d < d_cutoff:
								e.origin_vertex.data.mixed = True
								e.destination_vertex.data.mixed = True
						else:
							e.origin_vertex.data.mixed = True
							e.destination_vertex.data.mixed = True

class convexHull3D:

	def __init__(self, list_of_vertices, minArea=zero, sortBy="x", noChangeID=0, for_2d_triangulation=False): 

		#A list of Cartesian coordinates in dimension 3.
		self.list_of_vertices = list_of_vertices
		self.non_mutated_vertices = list_of_vertices
		self.edges = {}
		self.horizon_edges = {}
		self.horizon_cone = {}
		self.i = 0
		self.loop_count = 0
		self.number_of_cumulative_search_steps = 0

		# Flag for if the convex hull is going to be used for a 2D triangulation
		self.for_2d_triangulation = for_2d_triangulation

		# Minimum facet area for collinearity testing.
		self.minArea = minArea
		cg.zero = self.minArea # ZZZZ

		self.encountered = {}

		self.euler_v = 0
		self.euler_e = 0
		self.euler_t = 0

		"""The class attribute self.start_edge is critical to the
		   efficiency of the calculation of the convex hull.  It
		   ensures that traversals of the edges of the current
		   convex hull (the hull changes throughout the calculation)
		   begin within and remain within the edges that comprise
		   the current convex hull.  Rather than deleting edges that
		   are not part of the convex hull these edges are simply
		   disconnected in the function self.build_cone.  As a result,
		   the class attribute self.edges contains a subset of old
		   edges that are not part of the current convex hull, and a
		   subset of edges, old and new, that are a part of the
		   convex hull.  The attribute self.start_edge is an edge
		   that is always in the subset of new edges that are part
		   of the current convex hull.  For those class methods that
		   loop over the edges of the current convex hull, e.g., the
		   rotation functions, the attribute self.start_edges ensures
		   that the loop is initialized with an edge that is in the
		   current convex hull.
		"""
		self.start_edge = None

		#Call the class methods that calculate the hull.
		self.sort_vertices(sortBy=sortBy,noChangeID=noChangeID)
		self.make_simplex_3d()
		self.calculate_the_convex_hull()
		self.facets = None
		self.hyperPlanes = None

		#PPPP: Find the right place to put this reminder.
		#!!!!Be careful transforming the data.
		#Translating it by 100 gives the same answer.
		#Translating it by 1000 givea a slightly different answer,
		#because the coordinate become large, the determinant becomes
		#larger, which introduces errors into the mathematics.
		#See O'rourke for an explanation. The main point here is that
		#if the coordinates are large (how large?) the code may not
		#behave properly.  Should I take the time to account for this
		#using integer math?

		# The option for saving the 2D triangulation calculated from the 3D hull
		self.triangulation = None
	
	def sort_vertices(self, sortBy="x",noChangeID=0):
		temp_vertices = {}
		i = 1
		for v in self.list_of_vertices:
			nV = Vertex((v.x,v.y,v.z))
			# Assign new vertex ID.
			if not noChangeID:
				nV.id = i
			# Retain original vertex ID.
			else:
				nV.id = v.id
			try:
				nV.data = v.data
			except:
				nV.data = None

			"""Map the vertices by x-coordinate.
			   The y- and z-coordinate are included in the key of the
			   mapping to ensure that the key is unique."""
			if sortBy == "x":
				temp_vertices.update({(nV.x,nV.y,nV.z,i):nV})
			elif sortBy == "y":
				temp_vertices.update({(nV.y,nV.x,nV.z,i):nV})
			elif sortBy == "z":
				temp_vertices.update({(nV.z,nV.x,nV.y,i):nV})
			else:
				# Default is sort by x coordinate.
				temp_vertices.update({(nV.x,nV.y,nV.z,i):nV})
			i += 1

		# Sort the set of vertices by x-coordinate.
		keys = sorted(temp_vertices)

		# Replace the unsorted vertex list with the sorted vertex list.
		self.list_of_vertices = []
		for k in keys:
			self.list_of_vertices.append(temp_vertices[k])

		# Calculate zero based on coordinate scale
		self.minArea = geom_tol(self.list_of_vertices)
		global zero 
		zero = self.minArea

	def make_simplex_3d(self):

		v1, v2, v3, v4 = self.list_of_vertices[0], self.list_of_vertices[1], self.list_of_vertices[2], self.list_of_vertices[3]
		# v5, v6, v7, v8 = self.list_of_vertices[4], self.list_of_vertices[5], self.list_of_vertices[6], self.list_of_vertices[7] 

		# v2.x, v2.y, v2.z = v1.x, v1.y, v1.z
		# v3.x, v3.y, v3.z = v1.x, v1.y, v1.z
		# v4.x, v4.y, v4.z = v1.x, v1.y, v1.z
		# v5.x, v5.y, v5.z = v1.x, v1.y, v1.z
		# v6.x, v6.y, v6.z = v1.x, v1.y, v1.z
		# v7.x, v7.y, v7.z = v1.x, v1.y, v1.z
		# v8.x, v8.y, v8.z = v1.x, v1.y, v1.z

		# self.list_of_vertices = [v1, v2, v3, v4, v5, v6, v7, v8]

		# Ensure v1 and v2 are not identical
		d = distance(v1, v2)
		while not d:
			if self.for_2d_triangulation:
				v2.jostle(parabaloid=True)
			else:
				v2.jostle()		
			d = distance(v1, v2)
			if d:
				v1.no_more_jostling = True
				v2.no_more_jostling = True

		# Ensure v1, v2, and v3 are in general position
		vGp = 0
		while not vGp:
			vGp = gp2D(v1, v2, v3)
			if not vGp:
				if self.for_2d_triangulation:
					v3.jostle(parabaloid=True)
				else:
					v3.jostle()
			if vGp:
				v3.no_more_jostling = True

		# Ensure v1, v2, v3, and v4 are in general position
		vGp = 0
		while not vGp:
			vGp = gp3D(v1,v2,v3,v4)
			if not vGp:
				if self.for_2d_triangulation:
					v4.jostle(parabaloid=True)
				else:
					v4.jostle()
			if vGp:
				v4.no_more_jostling = True

		# Delete used vertices.
		#######################
		for v in [v1,v2,v3,v4]:
			self.list_of_vertices.remove(v)

		# Now that general position has been established,
		# build the initial hull in dimension 4. 
		#################################################

		# General position vertices.
		############################

		"""There are four vertices (0-faces) in the intial simplex:
		   0, 1, 2, 3.
		   There are six edges (1-faces) in the initial simplex:
		   01, 02, 03, 12, 13, 23.
		   There are four facets (d-1-faces) in the initial simplex:
		   012, 013, 023, 123.
		   There is one volume (d-face) in the intial simplex:
		   0123."""

		#Make the facets of the simplex.
		t012 = Triangle(v1, v2 , v3, VertexR=v4, zero=self.minArea, parabaloid=self.for_2d_triangulation)
		t012.make_ltriangle(); t012.set_edge_topology()
		self.edges.update(t012.get_edges())

		t023 = Triangle(v1, v3, v4, VertexR=v2, zero=self.minArea, parabaloid=self.for_2d_triangulation)
		t023.make_ltriangle(); t023.set_edge_topology()
		self.edges.update(t023.get_edges())

		t031 = Triangle(v1, v4, v2, VertexR=v3, zero=self.minArea, parabaloid=self.for_2d_triangulation)
		t031.make_ltriangle(); t031.set_edge_topology()
		self.edges.update(t031.get_edges())

		t123 = Triangle(v2, v3, v4, VertexR=v1, zero=self.minArea, parabaloid=self.for_2d_triangulation)
		t123.make_ltriangle(); t123.set_edge_topology()
		self.edges.update(t123.get_edges())

		t130 = Triangle(v2, v4, v1, VertexR=v3, zero=self.minArea, parabaloid=self.for_2d_triangulation)
		t130.make_ltriangle(); t130.set_edge_topology()
		self.edges.update(t130.get_edges())

		t231 = Triangle(v3, v4, v2, VertexR=v1, zero=self.minArea, parabaloid=self.for_2d_triangulation)
		t231.make_ltriangle(zero=self.minArea); t231.set_edge_topology()
		self.edges.update(t231.get_edges())

		#Randomly choose an edge from the simplex.
		#This edge is the start point for iterating over the edges
		#of the current convex hull.
		self.start_edge = self.edges[list(self.edges)[0]]
		self.horizon_cone = self.edge_iterator(self.start_edge)

	def calculate_the_convex_hull(self, for_2d_triangulation=False):

		# True of the convex hull is being used for a 2D triangulation
		if not for_2d_triangulation:
			for_2d_triangulation = self.for_2d_triangulation
		
		i_skip_vertices = True
		list_of_vertices = self.list_of_vertices
		additional_round = 0
		while i_skip_vertices:
			i_skip_vertices = False
			i, skipped_vertices = 0, []
			while i < len(list_of_vertices):
				
				# Get the next vertex in the set.
	            #################################
				vT, ivGp = list_of_vertices[i], 0
				
				number_of_search_steps = 0
				#XXXX: Explain what determines a horizon ridge.  Negative
				#volume, the vertex is below the plane of triangle face;
				#therefore, the face is a member of the hull.  Positive volume,
				#the vertex is above the plane of the triangle face; therefore
				#the face is visible and is not a member of the new hull.
				#PPPP: Move through the qes by topology. Start the traversal of
				#of the edges of the hull with a horizon ridge of the previous
				#hull.  Further explanation. Negative: its outside? Vice Versa?
				"""Traverse the set of horizon edges in a counterclockwise
				   direction. Guarantee that the right face of the Edge
				   object is the visible face of the horizon edge. This
				   ensures that traversal of the horizon edges is
				   counterclockwise in its direction, and it ensures that
				   the horizon cone is comprised of the right face of each
				   edge object."""
				# if not self.number_of_cumulative_search_steps:
				# 	print("here A")
				# 	edges = self.edge_iterator(self.start_edge)
				if additional_round:
					edges = self.edge_iterator(self.start_edge)
				else:
					"""The set of vertices is sorted by their x-ccordinate.
					   As a result, at least one edge of the horizon cone
					   built during the last iteration of this while loop
					   will be visible.  Because of this property, the
					   visibility of every edge of the convex hull need not
					   be considered, and the calculation of the hull becomes
					   much more efficient.

					   It is important to note that the horizon cone is
					   comprised of both the subset of horizon edges and the
					   subset of the new horizon cone edges from the last
					   iteration of this while loop.

					   It is also important to note that this added
					   efficiency is based on a key assumption:  if one or
					   more facets are visible from the current vertex, the
					   previous horizon cone contains an edge that is
					   incident upon one of these visible faces.  In other
					   words, it is assumed that the previous horizon cone
					   contains a current horizon edge.  This assumption
					   holds only when the set of vertices has been sorted
					   so that the next vertex is outside of the current
					   hull.
					"""
					edges = {}
					for e in self.horizon_cone:
						e = self.edges[e]
						edges[e.lf.e1.id] = e.lf.e1
						edges[e.lf.e2.id] = self.edges[e.lf.e2.id]
						edges[e.lf.e3.id] = self.edges[e.lf.e3.id]
						edges[e.rf.e1.id] = e.rf.e1
						edges[e.rf.e2.id] = self.edges[e.rf.e2.id]
						edges[e.rf.e3.id] = self.edges[e.rf.e3.id]

				iResults, skip_vertex_because_not_in_gp = {}, False
				for e in edges.values():

					number_of_search_steps += 1

					lfAboveBelow = e.lf.beneath_beyond(vT)
					rfAboveBelow = e.rf.beneath_beyond(vT)

					if ((not lfAboveBelow) or (not rfAboveBelow)):
						skip_vertex_because_not_in_gp = True
						break

					# Save results for classification in case the vertex is in GP.
					##############################################################
					iResults.update({e:((lfAboveBelow,e.lf),(rfAboveBelow,e.rf))})

				if skip_vertex_because_not_in_gp:
					i += 1
					if for_2d_triangulation:
						vT.jostle(parabaloid=True)
					else:
						vT.jostle()
					skipped_vertices.append(vT)
					continue

				# vT is in GP so classify iResults.
				####################################
				horizon_edge_found = False
				for e in iResults:

					# Unpack the result.
					####################
					f1AboveBelow = iResults[e][0][0]
					f1 = iResults[e][0][1]
					f2AboveBelow = iResults[e][1][0]
					f2 = iResults[e][1][1]

					# Classify result.
					# This loop finds the first horizon edge, then efficiently 
					# identifies the remaining horizon edges using a function 
					# call to traverse_the_horizon_edges(). The function returns
					# all of the new horizon edges and then "breaks" the loop.
					# Thus, the behavior of this loop is to progress until the
					# first horizon edge is encountered, then it is terminated.
					############################################################
					l, r = f1AboveBelow, f2AboveBelow
					if l != r:
						if l == 1:
							self.traverse_the_horizon_edges(vT,self.edges[e.flipid])
						else:
							self.traverse_the_horizon_edges(vT,e)
						self.number_of_cumulative_search_steps += number_of_search_steps
						horizon_edge_found = True
						break

				# When using the convex hull for a 2D triangulation,
				# This iterative condition make sure each x,y point maps to the convex hull of the parabaloid embedded in 3D
				if for_2d_triangulation:
					if not horizon_edge_found and additional_round:
						i += 1
						vT.jostle(parabaloid=True)
						skipped_vertices.append(vT)
						continue				

				# Iterate to next vertex index.
				###############################
				#print("Vertex added: ", i)
				i += 1
				if not additional_round:
					print("Vertex added: ", i+4)
				else:
					print("Jostled Vertex added:", i)
				#print("breaking here....")
				#raise SystemExit 

			print(f"Number of vertices skipped because not in general position: {len(skipped_vertices)}")
			list_of_vertices = skipped_vertices
			if list_of_vertices:
				additional_round = 1
				i_skip_vertices = True

			# if list_of_vertices:
			# 	if len(list_of_vertices) <= 2:
			# 		v = list_of_vertices[0]
			# 		print(v.id, v.x, v.y, v.z)

		self.test_eulers_rule()

	def traverse_the_horizon_edges(self,vertex,start_edge):

		self.horizon_edges = [start_edge,]

		"""Interestingly, counterclockwise rotation over the edges
		   incident on the origin vertex of the current horizon edge
		   to find the next horizon edge is more efficient than
		   clockwise rotation.  In other words, when rotating
		   counterclockwise over the edges incident on the origin
		   vertex, then next horizon edge is almost always
		   encountered in fewer iterations than if the search for
		   the next horizon edge proceeds in the clockwise direction.
		   Ties are not uncommon: the next horizon edge is encountered
		   in the same number of iterations whether the search proceeds
		   in either a counterclockwise or clockwise direction.
		"""
		ccw_origin_edges = self.rotate_ccw_about_the_origin_vertex(start_edge)

		i = 1
		e = ccw_origin_edges[i]
		ePrev = e.id
		while e.id != start_edge.flipid:
			"""This is the general solution for traversing the horizon
			   edges in a counterclockwise direction.  Given a horizon
			   edge (start_edge), rotate over the set of edges that are
			   incident upon the origin vertex of start_edge.  The
			   orientation of these edges is not predictable, i.e., each
			   each edge is a directed line that points back to,
			   or away from the origin vertex.  In this algorithm, it is
			   necessary to select for the edge that is a directed line
			   that points back to the origin vertex.  The direction of
			   line is subsequently reversed.  This ensures that
			   as the algorithm moves from horizon edge to horizon edge
			   it proceeds in a counterclockwise direction.  This also
			   ensures that the right face of each horizon edge is the
			   visible face, which simplifies the process of building
			   the new horizon cone"""
			if e.destination_vertex != self.edges[ePrev].origin_vertex:
				e = self.edges[e.flipid]
			ePrev = e.id

			l = self.edges[e.id].lf.beneath_beyond(vertex)
			r = self.edges[e.id].rf.beneath_beyond(vertex)

			if l != r:
				self.horizon_edges.append(self.edges[e.id])
				ccw_origin_edges =  self.rotate_ccw_about_the_origin_vertex(self.edges[e.id])
				i = 1
				e = ccw_origin_edges[i]
			else:
				i += 1
				e = ccw_origin_edges[i]

		self.build_the_horizon_cone(vertex)
		self.link_the_horizon_cone(vertex)

	def build_the_horizon_cone(self,vertex):
		updated_horizon_edges = [] 
		for e in self.horizon_edges:
			t = Triangle(e.lf.v1, e.lf.v2, e.lf.v3, VertexR=vertex, zero=self.minArea, skip_orientation=1)
			if t.id == "not in general position":
				while t.id == "not in general position":
					if self.for_2d_triangulation:
						vertex.jostle(parabaloid=True)
					else:
						vertex.jostle()
					t = Triangle(e.lf.v1, e.lf.v2, e.lf.v3, VertexR=vertex, zero=self.minArea, skip_orientation=1)
			t.make_ltriangle(zero=self.minArea); t.set_edge_topology()
			updated_edges = t.get_edges()
			updated_edge = updated_edges[e.id]
			updated_sym_edge = updated_edges[e.flipid]
			updated_horizon_edges.append(updated_edge)
			self.edges.update({updated_edge.id:updated_edge})
			self.edges.update({updated_sym_edge.id:updated_sym_edge})
		self.horizon_edges = updated_horizon_edges

	def link_the_horizon_cone(self,vertex):
		"""Copy the first horizon edge in the list self.horizon_edges
		   and append it to the end of the list.  This ensures that
		   the last two faces of the horizon cone are properly linked.
		"""
		self.horizon_edges.append(self.horizon_edges[0])

		i = 0
		#A dictionary of the new cone edges.
		horizon_cone = {}
		while i < len(self.horizon_edges)-1:
			#Update the new horizon cone with the current horizon edges.
			horizon_cone.update({self.horizon_edges[i].id:
                    			 self.horizon_edges[i]})

			#Origin vertex of the current edge.
			o = self.horizon_edges[i].origin_vertex
			#Destination vertex of the current edge.
			d = self.horizon_edges[i].destination_vertex
			#The new vertex of the hull. 
			n = vertex
			#The origin vertex of the next edge in a ccw direction.
			o_ccw = self.horizon_edges[i+1].origin_vertex

			#Orientation is important.
			#From horizon edges identify new directed edge (n,o).
			#Make the right face of this new edge first; therefore,
			#change the direction of the directed edge to (o.n).
			t = Triangle(o, n, o_ccw, VertexR=d, zero=self.minArea, skip_orientation=1)
			if t.id == "not in general position":
				while t.id == "not in general position":
					if self.for_2d_triangulation:
						n.jostle(parabaloid=True)
					else:
						n.jostle()
					t = Triangle(o, n, o_ccw, VertexR=d, zero=self.minArea, skip_orientation=1)
			t.make_ltriangle(zero=self.minArea); t.set_edge_topology()
			self.edges.update(t.get_edges())
			#Update the new horizon cone with the newly built edges of
			#the horizon cone. 
			horizon_cone.update(t.get_edges())
			i += 1

		"""Replace the old horizon cone with the newly built horizon
		   cone.  This update is at the heart of the efficiency of
		   the calculation.  
		"""
		self.horizon_cone = horizon_cone
		"""Relatedly, update the class attribute self.start_edge with
		   the first edge of the horizon edge.  This ensures that the
		   search for an initial horizon edge is limited to the
		   edges of the current hull.  Moreover, this ensures that
		   the search for an intial horizon edge is limited to the
		   edges of the horizon cone.
		"""
		self.start_edge = self.horizon_edges[0]

	def test_eulers_rule(self):
		#Euler's rule for polyhedra (i.e., d=3): V - E + F = 2.
		edges = self.edge_iterator(self.start_edge)

		nVertices = {}
		nEdges = {}
		nTriangles = {}
		for e in edges.values():
			nVertices.update({e.origin_vertex.id:e.origin_vertex})
			nVertices.update({e.destination_vertex.id:e.destination_vertex})
			nEdges.update({e.id:e})
			nTriangles.update({e.lf.id:e.lf})
			nTriangles.update({e.rf.id:e.rf})

		v = self.euler_v =  len(nVertices)
		e = self.euler_e = len(nEdges)/2.
		t = self.euler_t = len(nTriangles)

		# A successful test evaluates to 2.
		euler = v - e + t

		result = ""
		result += "\n%-14.14s% 4s %59.59s\n" % ("TEST OUTCOME:",
                                        		"",
                                    		   "test_eulers_rule")
		result += "%-79.79s\n" % (79*"-")
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF VERTICES OF THE HULL:",v)
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF EDGES OF THE HULL:",e)
		result += "%-40.40s%-39.0f\n" % ("NUMBER OF TRIANGLES OF THE HULL:",t)
		result += "%-40.40s%-39.0f\n" % ("EULER\'S RULE:",v - e + t)
		result += 2*"\n"
		return result

	def get_triangulation(self):

		# 1) Extract coordinates
		vertices = self.get_the_vertices_of_the_hull()
		xs = [v.x for v in vertices]
		ys = [v.y for v in vertices]
		zs = [v.z for v in vertices]

		# 2) Compute centroid in XY
		centroid_x = sum(xs) / len(vertices)
		centroid_y = sum(ys) / len(vertices)

		# 3) Find bottom-most z and subtract 10000
		bottom_z = min(zs)
		centroid_z = bottom_z - 10000

		# 4) Bottom reference Vertex
		bottom_vertex = Vertex((centroid_x, centroid_y, centroid_z))

		# 5) Create and return the 2D triangulation object of the 3D convex hull
		bottom_edges, edges = {}, self.getEdges()
		for e in edges:
			l = e.lf.beneath_beyond(bottom_vertex)
			r = e.rf.beneath_beyond(bottom_vertex)
			if l == 1 and r == 1:
				e.bottom = True
				bottom_edges[tuple(sorted(e.id))] = e
				#bottom_edges[tuple(sorted(e.id))] = self.rotate_ccw_about_the_origin_vertex(e)
		self.triangulation = triangulation2D(bottom_edges, self)
		return self.triangulation

		# for e in t.edges:
		# 	neighbors = self.rotate_ccw_about_the_origin_vertex()
		# 	for neighbor in neighbors:
		# 		origin = neighbor.origin_vertex
		# 		destination = neighbor.destination_vertex
		# 		origin_label = neighbor.origin_vertex.data.label
		# 		destination_label = neighbor.destination_vertex.data.label
		# 		if origin_label != destination_label:
		# 			diff += 1
		# 			diff_edges.add(neighbor)
		# 		else:
		# 			same += 1
		# 			same_edges.add(neighbor)

	def rotate_ccw_over_the_right_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Rnext
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Dprev
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Rnext
			n += 1
		return rotate

	def rotate_ccw_over_the_left_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Lnext
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Oprev
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Lnext
			n += 1
		return rotate

	def rotate_cw_over_the_right_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Rprev
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Onext
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Rprev
			n += 1
		return rotate

	def rotate_cw_over_the_left_face(self,start_edge):
		rotate = [start_edge,]
		next_edge = start_edge.Lprev
		n = 0
		while next_edge != start_edge.id:
			if n%2:
				next_edge = self.edges[next_edge].flipid
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Dnext
			else:
				rotate.append(self.edges[next_edge])
				next_edge = self.edges[next_edge].Lprev
			n += 1
		return rotate

	def rotate_ccw_about_the_origin_vertex(self,start_edge):
		#Rotate ccw over the edges that share the same origin vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Onext
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Onext
		rotate.append(start_edge)
		return rotate

	def rotate_ccw_about_the_destination_vertex(self,start_edge):
		#Rotate ccw over the edges that share the same destination vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Dnext
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Dnext
		return rotate

	def rotate_cw_about_the_origin_vertex(self,start_edge):
		#Rotate cw over the edges that share the same origin vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Oprev
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Oprev
		rotate.append(start_edge)
		return rotate

	def rotate_cw_about_the_destination_vertex(self,start_edge):
		#Rotate cw over the edges that share the same destination vertex.
		rotate = [start_edge,]
		next_edge = start_edge.Dprev
		n = 0
		while next_edge != start_edge.id:
			rotate.append(self.edges[next_edge])
			next_edge = self.edges[next_edge].Dprev
		return rotate

	# def edge_iterator(self,start_edge): # XXXX resolve with getEdges
	# 	"""If the unique_ids of two vertices match, this code will
	# 	   enter an infinite loop."""
	# 	self.encountered = {}
	# 	self.encountered.update({start_edge.id:self.edges[start_edge.id]})
	# 	self.encountered.update({start_edge.flipid:self.edges[start_edge.flipid]})
	# 	self.get_the_edges_of_the_hull(start_edge)
	# 	edges = self.encountered
	# 	self.encountered = {}
	# 	#print("NUMBER OF ACCESSED EDGES: ", self.iterator_count)
	# 	return edges

	def edge_iterator(self, start_edge):
	    """
	    Collect all unique half-edges reachable from `start_edge`.
	    Returns a {edge_id: Edge} mapping. Non-recursive, no instance side-effects.
	    """
	    seen = set()
	    out = {}

	    stack = [start_edge.id, start_edge.flipid]
	    while stack:
	        eid = stack.pop()
	        if eid in seen:
	            continue
	        seen.add(eid)

	        e = self.edges[eid]
	        out[eid] = e

	        # Walk the Rnext ring of this face once, enqueueing unseen edges and twins
	        ring_start = eid
	        nid = e.Rnext
	        while nid not in (ring_start, self.edges[ring_start].flipid):
	            if nid not in seen:
	                stack.append(nid)
	            twin = self.edges[nid].flipid
	            if twin not in seen:
	                stack.append(twin)
	            nid = self.edges[nid].Rnext

	    return out

	# def edge_iterator(self, start_edge, max_steps=None):
	#     """
	#     Traverse all unique half-edges reachable from start_edge using only:
	#       - twin links (flipid)
	#       - face membership (lf/rf -> e1,e2,e3)
	#     Returns {edge_id: Edge}. Robust against malformed Rnext rings.
	#     """
	#     seen = set()
	#     out = {}
	#     stack = [start_edge.id, start_edge.flipid]

	#     # Reasonable cap: ~10x the edge count, or a fixed large number if unknown
	#     EDGE_COUNT = len(self.edges) or 1
	#     budget = max_steps if max_steps is not None else (10 * EDGE_COUNT + 1000)
	#     steps = 0

	#     while stack:
	#         eid = stack.pop()
	#         if eid in seen:
	#             continue
	#         if eid not in self.edges:
	#             continue  # defensive: skip stale pointers

	#         seen.add(eid)
	#         e = self.edges[eid]
	#         out[eid] = e

	#         # explore twin
	#         if e.flipid in self.edges and e.flipid not in seen:
	#             stack.append(e.flipid)

	#         # explore incident faces (each should contribute its 3 half-edges)
	#         for face in (getattr(e, "lf", None), getattr(e, "rf", None)):
	#             if not face:
	#                 continue
	#             for he in (face.e1, face.e2, face.e3):
	#                 hid = he.id
	#                 if hid in self.edges and hid not in seen:
	#                     stack.append(hid)
	#                 twin = self.edges.get(hid).flipid if hid in self.edges else None
	#                 if twin in self.edges and twin not in seen:
	#                     stack.append(twin)

	#         steps += 1
	#         if steps > budget:
	#             raise RuntimeError("edge_iterator: step budget exceeded; topology likely corrupted")

	#     return out

	# def edge_iterator(self, start_edge, max_steps=None):
	#     """
	#     Traverse all unique half-edges reachable from start_edge using only
	#     neighbors that are defined for stored e1 edges:
	#       flipid, Rnext/Rprev, Dnext/Dprev, Onext/Oprev, Lnext/Lprev.

	#     Returns: {edge_id: Edge}
	#     """
	#     seen = set()
	#     out = {}
	#     stack = [start_edge.id]  # start from the given edge id

	#     # Safe upper bound, avoids infinite traversal on corrupted topology
	#     edge_count = max(len(self.edges), 1)
	#     budget = max_steps if max_steps is not None else (10 * edge_count + 1000)
	#     steps = 0

	#     while stack:
	#         eid = stack.pop()
	#         if eid in seen:
	#             continue
	#         e = self.edges.get(eid)
	#         if e is None:
	#             continue  # stale pointer; skip defensively

	#         seen.add(eid)
	#         out[eid] = e

	#         # Only push neighbors that exist in self.edges
	#         neighbors = (
	#             e.flipid, e.Rnext, e.Rprev, e.Dnext, e.Dprev,
	#             e.Onext, e.Oprev, e.Lnext, e.Lprev
	#         )
	#         for nid in neighbors:
	#             if nid in self.edges and nid not in seen:
	#                 stack.append(nid)

	#         steps += 1
	#         if steps > budget:
	#             raise RuntimeError(
	#                 f"edge_iterator: step budget exceeded "
	#                 f"(visited {len(seen)} / ~{edge_count}). "
	#                 f"Topology likely corrupted near edge {eid}."
	#             )

	#     return out

	# def edge_iterator(self, start_edge=None, strict=False):
	#     """
	#     Return a mapping {edge_id: Edge} of well-formed half-edges by
	#     scanning self.edges. No Rnext/Onext/Lnext walking is used.

	#     A 'well-formed' edge passes:
	#       - twin symmetry: flip exists and flips back
	#       - faces set on both sides (lf and rf are not None)
	#       - non-degenerate id (origin != destination)
	#     If 'strict', also asserts that the edge id appears on its faces.
	#     """

	#     out = {}
	#     bad = 0

	#     for eid, e in self.edges.items():
	#         # 1) Non-degenerate (avoid self-loops from duplicate vertex ids)
	#         if e.origin_vertex.id == e.destination_vertex.id:
	#             bad += 1
	#             continue

	#         # 2) Twin symmetry
	#         twin = self.edges.get(e.flipid)
	#         if twin is None or twin.flipid != eid:
	#             bad += 1
	#             continue

	#         # 3) Faces exist
	#         if getattr(e, "lf", None) is None or getattr(e, "rf", None) is None:
	#             bad += 1
	#             continue

	#         if strict:
	#             # 4) (Optional) The edge id shows up on its own faces
	#             def face_has_edge(face, eid_check):
	#                 try:
	#                     ids = (face.e1.id, face.e2.id, face.e3.id)
	#                 except Exception:
	#                     return False
	#                 return eid_check in ids
	#             if not face_has_edge(e.lf, eid) and not face_has_edge(e.lf, e.flipid):
	#                 bad += 1
	#                 continue
	#             if not face_has_edge(e.rf, eid) and not face_has_edge(e.rf, e.flipid):
	#                 bad += 1
	#                 continue

	#         out[eid] = e

	#     # Optional: if you want visibility during debugging
	#     # print(f"edge_iterator: kept {len(out)} edges, filtered {bad}")

	#     return out

	# def edge_iterator(self, start_edge, max_factor=10):
	#     """
	#     Collect all reachable half-edges starting from start_edge.
	#     Uses only IDs that actually exist in self.edges.
	#     Bounded and cycle-safe; never infinite-loops even on corrupt rings.
	#     Returns: {edge_id: Edge}
	#     """
	#     edges_map = self.edges
	#     if not edges_map:
	#         return {}

	#     seen = set()
	#     out = {}
	#     stack = [start_edge.id, getattr(start_edge, "flipid", None)]
	#     stack = [eid for eid in stack if eid in edges_map]

	#     # Hard safety cap: proportional to current edge set size
	#     nE = max(len(edges_map), 1)
	#     budget = max_factor * nE + 1024
	#     steps = 0

	#     # Neighbor fields that MAY be valid on stored e1 half-edges
	#     neighbor_fields = ("flipid", "Rnext", "Rprev", "Onext", "Oprev",
	#                        "Dnext", "Dprev", "Lnext", "Lprev")

	#     while stack:
	#         eid = stack.pop()
	#         if eid in seen or eid not in edges_map:
	#             continue

	#         # Walk the local right-face ring safely (guards inside)
	#         for rid in self._walk_ring_safe(eid, next_attr="Rnext", cap=64):
	#             if rid not in seen and rid in edges_map:
	#                 seen.add(rid)
	#                 out[rid] = edges_map[rid]
	#                 # cross to twin for later exploration
	#                 twin = edges_map[rid].flipid
	#                 if twin in edges_map and twin not in seen:
	#                     stack.append(twin)

	#         # Also enqueue direct neighbors (only those that exist)
	#         e = edges_map[eid]
	#         for f in neighbor_fields:
	#             nid = getattr(e, f, None)
	#             if nid in edges_map and nid not in seen:
	#                 stack.append(nid)

	#         steps += 1
	#         if steps > budget:
	#             # Bail out with what we have; topology is likely corrupted
	#             break

	#     return out

	def _walk_ring_safe(self, start_eid, next_attr="Rnext", cap=64):
	    """
	    Yield a ring by following 'next_attr' pointers, but stop if:
	      - an edge repeats in this ring, or
	      - pointer is missing, or
	      - we exceed 'cap' steps.
	    This prevents infinite loops on malformed rings.
	    """
	    edges_map = self.edges
	    if start_eid not in edges_map:
	        return

	    seen_in_ring = set([start_eid])
	    # include twin to handle older stop conditions gracefully
	    twin = edges_map[start_eid].flipid
	    if twin is not None:
	        seen_in_ring.add(twin)

	    nxt = getattr(edges_map[start_eid], next_attr, None)
	    steps = 0

	    while nxt in edges_map and nxt not in seen_in_ring and steps < cap:
	        yield nxt
	        seen_in_ring.add(nxt)
	        nxt = getattr(edges_map[nxt], next_attr, None)
	        steps += 1
	    # silent stop on loop/missing pointer/exceeded cap



	# def edge_iterator(self, start_edge):
	#     """
	#     Traverse all unique half-edges reachable from start_edge
	#     using only face membership (lf/rf with e1/e2/e3) and flipid.
	#     Returns {edge_id: Edge}.
	#     """
	#     seen = set()
	#     out = {}
	#     stack = [start_edge.id, start_edge.flipid]

	#     while stack:
	#         eid = stack.pop()
	#         if eid in seen:
	#             continue
	#         seen.add(eid)

	#         e = self.edges[eid]
	#         out[eid] = e

	#         # 1) Cross to the adjacent face via the twin
	#         if e.flipid not in seen:
	#             stack.append(e.flipid)

	#         # 2) From each incident face, enqueue its three half-edges and their twins
	#         for face in (e.lf, e.rf):
	#             if not face:
	#                 continue
	#             for he in (face.e1, face.e2, face.e3):
	#                 hid = he.id
	#                 if hid not in seen:
	#                     stack.append(hid)
	#                 twin = self.edges[hid].flipid
	#                 if twin not in seen:
	#                     stack.append(twin)

	#     return out

	# 2012.10.09: Legacy function for collecting the edges of the hull.
	#             I will keep for now because of dependencies elsewhere. 
	#             I will eventually replace this call with getEdges.
	def get_the_edges_of_the_hull(self,start_edge): # XXXX
		e = start_edge.Rnext
		while e not in (start_edge.id, start_edge.flipid):
			e     = self.edges[e].Rnext
			esym  = self.edges[e].flipid
			if e not in self.encountered.keys():
				self.encountered.update({e:self.edges[e]})
				self.encountered.update({esym:self.edges[esym]})              
				self.get_the_edges_of_the_hull(self.edges[e])
				self.get_the_edges_of_the_hull(self.edges[esym])

	def get_the_vertices_of_the_hull(self):

		# Map all of the vertices used to calculate the convex hull.
		all_verts = {}
		for v in self.non_mutated_vertices:
			all_verts.update({v.id:v})

		# Map the vertices of the convex hull.
		temp_verts_of_hull = {}
		edges = self.edge_iterator(self.start_edge)
		for e in edges.values():
			ov = e.origin_vertex
			dv = e.destination_vertex
			temp_verts_of_hull.update({ov.id:ov, dv.id:dv})
			if ov.id in all_verts: del all_verts[ov.id]
			if dv.id in all_verts: del all_verts[dv.id]

		keys = sorted(temp_verts_of_hull)
		verts_of_hull = []
		for key in keys:
			verts_of_hull.append(temp_verts_of_hull[key])
		return verts_of_hull 

	def get_the_vertices_inside_of_the_hull(self):

		# hull_verts = self.get_the_vertices_of_the_hull()
		# inside_verts = {}
		# #for v in self.list_of_vertices:
		# for v in self.non_mutated_vertices:
		# 	if v not in hull_verts:
		# 		inside_verts.update({v.id:v})

		inside_verts = {}
		for v in self.list_of_vertices:
			result = self.test_if_this_vertex_is_in_the_hull(v)
			# print(result, v.data)
			if result:
				inside_verts.update({v.id:v})
				# print(v.data)
		# for v in self.non_mutated_vertices:
		# 	if v not in hull_verts:
		# 		inside_verts.update({v.id:v})

		keys = sorted(inside_verts)
		verts_inside_of_hull = []
		for key in keys:
			verts_inside_of_hull.append(inside_verts[key])
		return verts_inside_of_hull

	def empty_hull(self,vertex_list):
		for vertex in vertex_list:
			try:
				if self.test_if_this_vertex_is_in_the_hull(vertex):
					"""Return false if the hull is empty"""
					return 0
			except:
				print("\n1--The test vertex is not in general position",)
				print("relative to the hull. " )
				print("2--Attempt to jostle the point.")
				vertex.jostle()
				#Try to jostle the point.
				#from random import random
				#vertex.x += (random()/5.)
				#vertex.y += (random()/5.)
				#vertex.z += (random()/5.)
				try:
					print("3--Attempt to jostle was a success.")
					if self.test_if_this_vertex_is_in_the_hull(vertex):
						return 0
				except:
					print("3--The test vertex is STILL not in general position",)
					print("relative to the hull." )
					print("4--Returning false.\n")
					return 0
		"""Return true if the hull is empty"""
		return 1

	# 2012.10.09: Legacy function for collecting the facets of the hull.
	#             I will keep for now because of dependencies elsewhere. 
	#             I will eventually replace this call with getFacets.
	def get_the_facets_of_the_hull(self):
		facets = {}
		edges = self.edge_iterator(self.start_edge)
		for e in edges.values():
			leftFacetId = list(e.lf.id)
			leftFacetId.sort()
			rightFacetId = list(e.rf.id)
			rightFacetId.sort()
			facets.update({tuple(leftFacetId):e.lf})
			facets.update({tuple(rightFacetId):e.rf})
		return facets.values()

	# 2012.10.09: New, more concise class method call. 
	def getEdges(self):
		return self.edge_iterator(self.start_edge).values()

	# 2012.10.09: New, more concise class method call.
	def getFacets(self):
		if self.facets:
			return self.facets
		else:
			facets = {}
			edges = self.edge_iterator(self.start_edge)
			for e in edges.values():
				leftFacetId = list(e.lf.id)
				leftFacetId.sort()
				rightFacetId = list(e.rf.id)
				rightFacetId.sort()
				facets.update({tuple(leftFacetId):e.lf})
				facets.update({tuple(rightFacetId):e.rf})
			self.facets = facets.values()
			return self.facets

	def getHullAsString(self):
		facets, facetString = self.getFacets(), ""
		for facet in facets:
			fString = "%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n"
			facetString += fString % (facet.v1.x, facet.v1.y, facet.v1.z,facet.v2.x, facet.v2.y, facet.v2.z, facet.v3.x, facet.v3.y, facet.v3.z) 
		return facetString


	def getHyperPlanes(self):
		if self.hyperPlanes:
			return self.hyperPlanes
		else:
			#from compGeometry import planeCoefficients3D
			hyperPlanes, facets = {}, self.getFacets()
			for facet in facets:
				c = planeCoefficients3D(facet.v1,facet.v2,facet.v3)
				hyperPlanes.update({facet:c})
			self.hyperPlanes = hyperPlanes
			return self.hyperPlanes

	# def test_if_this_vertex_is_in_the_hull(self,v):
	# 	inside, hyperPlanes = 1, self.getHyperPlanes()
	# 	for facet in hyperPlanes:
	# 		hP = hyperPlanes[facet]
	# 		sign = hP[0]*v.x + hP[1]*v.y + hP[2]*v.z + hP[3]
	# 		# The sign convention is the opposite of what it seems 
	# 		# because of the manner in which the topology of the hull
	# 		# facets has been arranged. 
	# 		if sign > 0: # if sign < 0: # original, but likely wrong topology
	# 			inside = 0
	# 			break
	# 	return inside

	def test_if_this_vertex_is_in_the_hull(self, v, returnVolume=0):

		# Only need to test beneath_beyond for the left face because all edges (i.e., the e.flipid edges) are accessed in the loop.
		# In other words, in the quad edge, what is left face for e.id is right face for e.flipid, and vice versa. 
		###########################################################################################################################
		minVolume = 100000.
		edges = self.getEdges()
		for e in edges:
			l = e.lf.beneath_beyond(v)
			if l > 0: 
				return 0
			elif l < 0 and returnVolume:
				volume = e.lf.beneath_beyond(v, returnVolume=returnVolume)
				if volume < minVolume:
					minVolume = volume

		if returnVolume: 
			return minVolume
		else:
			return 1

	def modulateHull(self, multiplier=1, swell=1, shrink=0):

		from pHinder._vendor.pdbFile import PseudoAtom
		from pHinder._vendor.compGeometry import centroid3D, crossProductCoordinates

		edges = self.getEdges()
		swellFacetsString = ""
		vs = {}
		for edge in edges:

			e = edge
			psa = PseudoAtom()
			v1, v2, v3 = e.rf.v1, e.rf.v2, e.rf.v3
			c = centroid3D([v1, v2, v3], psa)

			psa1 = PseudoAtom()
			if swell:
				psaResult = crossProductCoordinates(v2, psa.v, v1, psa1, multiplier=multiplier, norm=1)
				psaResult.atom_serial_number = v2.data.atom_serial
				psaResult.v.id = v2.data.atom_serial
			else:
				psaResult = crossProductCoordinates(v1, psa.v, v2, psa1, multiplier=multiplier, norm=1)
				psaResult.atom_serial_number = v1.data.atom_serial
				psaResult.v.id = v1.data.atom_serial
			vs[psaResult.v.id] = psaResult.v

		# Calculates and returns a new convex hull instance...
		vs = list(vs.values())
		return convexHull3D(vs)

	def get_LO(self):
		edges = self.edge_iterator(self.start_edge)
		result_class = ConvexHull3D_LO(edges)
		result_class.of_hull = self.get_the_vertices_of_the_hull()
		result_class.in_hull = self.get_the_vertices_inside_of_the_hull()
		return result_class

class ConvexHull3D_LO:

    #THE SIMPLEST VERSION OF A CONVEX HULL INSTANCE THAT REDUCES OVERHEAD.
    
    def __init__(self,edges):
        self.edges   = edges
        self.of_hull = None
        self.in_hull = None

    def test_if_this_vertex_is_in_the_hull(self,vertex):
        for e in self.edges.values():
            l = e.lf.beneath_beyond(vertex)
            r = e.rf.beneath_beyond(vertex)
            if l+r != -2 and l+r != 0:
                return 0
        return 1



        
