global zero
zero = 1e-12

# Import dependencies here for enhanced performance.
####################################################
from random import choice
from math import sqrt, cos, sin, acos, degrees
from pHinder._vendor.determinants import det3x3, det4x4
from decimal import *
import numpy as np
from collections.abc import Sequence, Iterable
import hashlib
import os
import random
import struct
from typing import Any, Tuple, Union, Optional

# ---------------------------------------------------------------------------
# Where jostling gets its randomness
#
# Jostling exists to put points in general position: without it the orientation
# predicates are ambiguous at exact zero and the rotation around facets can
# cycle instead of terminating. That is not in question here and none of it
# changes below -- same epsilon, same uniform draw in [-eps, +eps], same call
# sites, same number of draws in the same order.
#
# The only thing that changes is where the numbers come from. Drawing them from
# the process-global unseeded RNG made every run of the same structure produce a
# slightly different answer. Drawing them from an RNG seeded on the vertex's own
# identity and coordinates and its jostle count makes the perturbation a
# deterministic function of the input: the same structure gets the same nudge,
# each retry still gets a different one, two vertices at the same point still
# come apart, and degeneracy is escaped exactly as before.
#
# TO REVERT, either of these restores the previous behaviour exactly:
#     * set the environment variable PHINDER_JOSTLE=random
#     * or set JOSTLE_DETERMINISTIC = False below
# Nothing else in the file depends on the choice.
# ---------------------------------------------------------------------------

JOSTLE_DETERMINISTIC = os.environ.get("PHINDER_JOSTLE", "deterministic").lower() != "random"


# The Vertex/Vertex4D default for unique_id.  A vertex carrying this has no
# identity of its own, and is seeded from position and count alone -- exactly as
# before identity was folded in -- so anything that does not name its vertices
# keeps the behaviour it had.
_NO_IDENTITY = "no id"


def _jostle_rng(x, y, z, count, identity=_NO_IDENTITY):
    """An RNG determined by which vertex this is, where it is, and how often it
    has been moved.

    Position and count alone are not enough. They are *identical* for two
    vertices that occupy the same point, so both would draw the same nudge, move
    together, and stay coincident however many times they were jostled -- the one
    degeneracy a perturbation exists to break, and the one case the unseeded
    global RNG handled for free by separating them on the first call.

    Mixing the vertex's own id in restores that separation without giving up
    determinism: the id is a property of the input (in the hulls it is the
    vertex's position in a sort of the input coordinates), not of the process, so
    the same structure still yields the same perturbations run after run.

    blake2b rather than hash() so the seed does not depend on the interpreter's
    hash seed, the platform, or the Python version.
    """
    payload = struct.pack("<dddq", x, y, z, count)
    if identity != _NO_IDENTITY:
        payload += repr(identity).encode("utf-8")
    seed = int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little")
    return random.Random(seed)



# Determine the appropriate value for zero from a point set
def geom_tol(*points: Any, factor: float = 10.0) -> float:
    """
    Estimate a near-zero tolerance from a set of 3D points.
    Accepts either:
      - geom_tol(v1, v2, v3, ...)
      - geom_tol([v1, v2, v3, ...])
    Each point can be an object with .x/.y/.z or a sequence/array of length 3.
    """
    def unpack_inputs(objs):
        # If single iterable (and not a point with x/y/z), unpack it
        if len(objs) == 1:
            single = objs[0]
            if (isinstance(single, Iterable)
                and not (hasattr(single, "x") and hasattr(single, "y") and hasattr(single, "z"))
                and not (isinstance(single, (str, bytes)))):
                return list(single)
        return list(objs)

    def point_to_array(p):
        if hasattr(p, "x") and hasattr(p, "y") and hasattr(p, "z"):
            return np.array([p.x, p.y, p.z], dtype=float)
        arr = np.asarray(p, dtype=float).ravel()
        if arr.size == 3:
            return arr.copy()
        raise TypeError(f"Cannot interpret {p!r} as a 3D point")

    pts = unpack_inputs(points)
    if not pts:
        raise ValueError("Need at least one point to compute tolerance")

    coords = np.vstack([point_to_array(p) for p in pts])  # shape (n,3)
    scale = max(np.max(np.abs(coords)), 1.0)  # floor to 1.0 to avoid too-small scales
    # print("------------------------->", factor * scale * np.finfo(float).eps)
    return factor * scale * np.finfo(float).eps


# Classes for basic geometry objects.
#####################################################################################################################

class Vertex:
    
    """DESCRIPTION: This class represents a point in dimension 3.
    """
    def __init__(self,coordinate,data=None,unique_id="no id"):
        """DESCRIPTION OF ARGUMENTS:

           coordinates: a tuple of Cartesian coordinates (x,y,z).

           data: an optional attribute of any type that contains
                 supplementary information about the vertex.

           unique_id: an optional attribute of any type that is unique
                      to only this vertex.
        """
        self.x = coordinate[0]
        self.y = coordinate[1]
        self.z = coordinate[2]
        self.data = data
        self.id = unique_id
        self.coordinate_tuple = (self.x,self.y,self.z)
        self.nJostles = 0
        self.no_more_jostling = False

    def reinitialize(self):
        self.coordinate_tuple = (self.x,self.y,self.z)
        self.data.x, self.data.y, self.data.z = self.x, self.y, self.z
        self.data.v = self

    def jostle(self, eps=1e-3, parabaloid=False):#, perturbLevel=1):
        """
        Perturb x and y by a random amount in [–eps, +eps],
        and—for a paraboloid lift—recompute z if desired.
        """
        if not self.no_more_jostling:
            if JOSTLE_DETERMINISTIC:
                self._jostle_count = getattr(self, "_jostle_count", 0) + 1
                r = _jostle_rng(self.x, self.y, self.z, self._jostle_count,
                                self.id).random
            else:
                r = random.random  # local ref is slightly faster than full lookup
            self.x += (r()*2 - 1) * eps
            self.y += (r()*2 - 1) * eps

            if parabaloid:
                # If you’re using z = x**2 + y**2, update it here:
                self.z = self.x**2 + self.y**2
            else:
                self.z += (r()*2 - 1) * eps

    # def jostle_old(self, perturbLevel=1):

    #     # 2022.03.08
    #     # These functions should now not be necessary very often...
    #     # Because general position is ensured in the convexHull4D_1_22 constructor...
    #     # They may be necessary in surface refinement...i.e., highResolutionSurface=1
        
    #     self.nJostles += 1

    #     # Use the random module to generate random signs and dividends.
    #     ###############################################################
    #     signs = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    #     if perturbLevel == 0:
    #         seq = [100., 500., 1000.]
    #     elif perturbLevel == 1:
    #         seq = [10., 20., 30., 40., 50., 60., 70., 80., 90.]
    #     elif perturbLevel == 2:
    #         seq = [2., 3., 4., 5., 6., 7., 8., 9., 10.]

    #     # Jostle X.
    #     ###########
    #     factor = choice(seq)
    #     sign = choice(signs)
    #     if sign:
    #         self.x += 1./factor
    #     else:
    #         self.x -= 1./factor

    #     # Jostle Y.
    #     ###########
    #     factor = choice(seq)
    #     sign = choice(signs)
    #     if sign:
    #         self.y += 1./factor
    #     else:
    #         self.y -= 1./factor

    #     # # Jostle Z.
    #     # ###########
    #     # factor = choice(seq)
    #     # sign = choice(signs)
    #     # if sign:
    #     #     self.z += 1./factor
    #     # else:
    #     #     self.z -= 1./factor

    #     # Reinitialize the vertex.
    #     ##########################
    #     self.reinitialize()

#	def __repr__(self):
#        
#		return `self.id`

class Vertex4D:
    
    """DESCRIPTION: This class represents a point in dimension 4.
    """
    
    def __init__(self,coordinate,data=None,unique_id="no id",setC=0):
        
        """DESCRIPTION OF ARGUMENTS:

           coordinates: a tuple of coordinates (x,y,z,u).

           data: an optional attribute of any type that contains
                 supplementary information about the vertex.

           unique_id: an optional attribute of any type that is unique
                      to only this vertex.
        """
        
        # Original unjostled coordinates.
        self.xo = coordinate[0]
        self.yo = coordinate[1]
        self.zo = coordinate[2]
        self.uo = coordinate[3]

        # Coordinates that represent the object and may have been jostled.
        self.x = coordinate[0]
        self.y = coordinate[1]
        self.z = coordinate[2]
        self.u = coordinate[3]


        self.setC = setC
        #if setC:
        #	self.cx = c_double(self.x)
        #	self.cy = c_double(self.y)
        #	self.cz = c_double(self.z)
        #	self.cu = c_double(self.u)
        #	self.cxp = pointer(self.cx)
        #	self.cyp = pointer(self.cy)
        #	self.czp = pointer(self.cz)
        #	self.cup = pointer(self.cu)
        #	point4D = c_double * 4
        #	self.cPoint = point4D(self.x,self.y,self.z,self.u)
        #	self.cPointer = pointer(self.cPoint)
        self.data = data
        self.id = unique_id
        self.coordinate_tuple = (self.x,self.y,self.z,self.u)
        self.nJostles = 0
        self.core = 0
        self.margin = 0
        self.exposed = 0 

    def reinitialize(self):
        self.u = self.x**2 + self.y**2 + self.z**2
        self.coordinate_tuple = (self.x,self.y,self.z,self.u)
        self.data.x, self.data.y, self.data.z = self.x, self.y, self.z
        self.data.v = self
        #if self.setC:
        #	point4D = c_double * 4
        #	self.cPoint = point4D(self.x,self.y,self.z,self.u)
        #	self.cPointer = pointer(self.cPoint)

    def jostle(self, perturbLevel=1):

        # 2022.03.08
        # These functions should now not be necessary very often...
        # Because general position is ensured in the convexHull4D_1_22 constructor...
        # They may be necessary in surface refinement...i.e., highResolutionSurface=1

        # print("Vertex4D ", self.id, " jostling...")
        
        self.nJostles += 1

        # Same switch as Vertex.jostle above: deterministic by default, drawn
        # from this vertex's position and its jostle count, so the sequence of
        # perturbations is a function of the input rather than of the process.
        # Set PHINDER_JOSTLE=random (or JOSTLE_DETERMINISTIC = False) to go back
        # to the global unseeded choice().
        _choice = (_jostle_rng(self.x, self.y, self.z, self.nJostles,
                               self.id).choice
                   if JOSTLE_DETERMINISTIC else choice)

        # Use the random module to generate random signs and dividends.
        ###############################################################
        signs = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
        if perturbLevel == 0:
            seq = [100., 500., 1000.]
        elif perturbLevel == 1:
            seq = [10., 20., 30., 40., 50., 60., 70., 80., 90.]
        elif perturbLevel == 2:
            seq = [2., 3., 4., 5., 6., 7., 8., 9., 10.]

        # Jostle X.
        ###########
        factor = _choice(seq)
        sign = _choice(signs)
        if sign:
            self.x += 1./factor
        else:
            self.x -= 1./factor

        # Jostle Y.
        ###########
        factor = _choice(seq)
        sign = _choice(signs)
        if sign:
            self.y += 1./factor
        else:
            self.y -= 1./factor

        # Jostle Z.
        ###########
        factor = _choice(seq)
        sign = _choice(signs)
        if sign:
            self.z += 1./factor
        else:
            self.z -= 1./factor

        # Jostle U.
        ###########
        self.u = self.x**2 + self.y**2 + self.z**2

        # # Measure the degree of the jostle's perturbation. # XXXX 3/9/22
        # ##################################################
        # d = sqrt((self.x-self.xo)**2 + (self.y-self.yo)**2 + (self.z-self.zo)**2)

        # # Reject the jostle. # XXXX 3/9/22
        # ####################
        # if d > 0.1:
        #     self.x, self.y, self.z = self.xo, self.yo, self.zo
        # Accept the jostle.
        # ####################
        # else:
        #     # Reinitialize the vertex.
        #     ##########################
        self.reinitialize()

#    def __repr__(self):
#        
#        return `self.id`

class Edge:
    
    """DESCRIPTION: This class represents a directed edge comprised of
                    two Vertex objects.
    """
    def __init__(self, Vertex1, Vertex2):
        
        """DESCRIPTION OF ARGUMENTS:

           Vertex1 and Vertex2: the two vertex objects that comprise
                                the directed edge Vertex1-->--Vertex2. 
        """

        # The origin and destination vertices of the directed edge.
        ##########################################################
        self.origin_vertex = Vertex1
        self.destination_vertex = Vertex2

        # Calculate the midpoint of the edge.
        #####################################
        from pHinder._vendor.pdbFile import PseudoAtom
        psa = PseudoAtom()
        self.midPoint = centroid3D([self.origin_vertex, self.destination_vertex], psa)
        
        # Optional attribute to link the Edge object to its closet protein atom.
        ########################################################################
        self.closestAtom = None

        # Unique class attributes that identify the directed edge.
        ##########################################################
        self.id = (self.origin_vertex.id, self.destination_vertex.id)
        self.flipid = (self.destination_vertex.id, self.origin_vertex.id)

        # For triangulation calculations, set to True if it is a triangulation edge.
        ############################################################################
        self.bottom = False

         # For triangulation calculations, set to True if it is a surface edge.
        ############################################################################
        self.visible = False       
        self.surface = False

        # # For network clustering from triangulation calculations. 
        # # Set to True if the triangulation edge was trimmed from the network.
        # ############################################################################
        # self.trimmed = False

        """The following group of attributes describe the topology of
           the directed edge O-->--D, where O and D correspond to the
           attributes origin_vertex and destination_vertex, respectively.
           Directed edge O-->--D is a directed edge of a convex polytope
           in dimension 3.  A quad-edge datastructure describes the
           topology of the directed edges of a convex polytope
           regardless of the number of edges of each face of the
           polytope, i.e., the faces of the polytope can be triangles,
           pentagons, hexagons, etc.  A quad-edge data structure is
           comprised of the two faces (lf and rf) that share the
           directed edge O-->--D, and eight directed edges:

           COUNTERCLOCKWISE TOPOLOGY ATTRIBUTES
           Rnext ()
           Dnext ()
           Onext ()
           Lnext ()

           CLOCKWISE TOPOLOGY ATTRIBUTES
           Rprev ()
           Dprev ()
           Oprev ()
           Lprev ()

           Here is a shematic of the quad-edge datastructure:

                   COUNTERCLOCKWISE                 CLOCKWISE

                   Lnext       Dnext            Dprev       Rprev
                   a---<---D---<---c            a--->---D--->---c
                           |                            |
                      lf   ^   rf                  lf   ^   rf
                           |                            |
                   b---<---O---<---d            b--->---O--->---d
                   Onext       Rnext            Lprev       Oprev

           The characters ^, < and > indicate the directionality of each
           edge in the quad-edge.  In the special case that the faces of
           the convex polytopes are triangles, which is the case we are
           concerned with when calculating the convex hull of a set of
           points in dimension 3, vertices a and b are the same vertex,
           and vertices c and d are the same vertex.  REFINE: By organizing the
           directed edges of a convex polytope into a set of quad-edge
           datastructures, it is relative easy and efficient to rotate
           over the edges of a face or... 
        """
        self.lf, self.rf = None, None  
        self.Lnext, self.Dnext = None, None
        self.Onext, self.Rnext = None, None
        self.Dprev, self.Rprev = None, None
        self.Lprev, self.Oprev = None, None

    def length(self):
        return sqrt((self.origin_vertex.x-self.destination_vertex.x)**2+ \
                    (self.origin_vertex.y-self.destination_vertex.y)**2+ \
                    (self.origin_vertex.z-self.destination_vertex.z)**2)

    def __repr__(self): 
        return str(self.id)

    def topological__repr__(self):
        result = ""
        result += "\n%-30.30s%49.49s\n" % ("TEST EDGE TOPOLOGY:",
                                           "Edge.topological__repr__()")
        result += "%-79.79s\n" % (79*"-")
        result += "%-36s%-8s%35s\n" % ("LEFT FACE","EDGE","RIGHT FACE")
        result += "%-35s%-9s%35s\n" % (str(self.lf),self.id,str(self.rf))
        result += "%-79.79s\n" % (79*"-")
        result += "%-33s%-11s%35s\n" % ("","CW CIRCUIT","")
        result += "%-21s%-19s%20s%19s\n" % ("DPREV","LPREV","RPREV","OPREV")
        result += "%-20s%-20s%20s%19s\n" % (str(self.Dprev),str(self.Lprev),
                                            str(self.Rprev),str(self.Oprev))
        result += "%-33s%-11s%35s\n" % ("","CCW CIRCUIT","")
        result += "%-21s%-19s%20s%19s\n" % ("ONEXT","LNEXT","RNEXT","DNEXT")
        result += "%-20s%-20s%20s%19s\n" % (str(self.Onext),str(self.Lnext),
                                            str(self.Rnext),str(self.Dnext))
        result += 2*"\n"
        return result

class Triangle:

    """DESCRIPTION: This class represents a triangle with an oriented
                    topology.  The vertices of the triangle are oriented
                    in a counterclockwise direction when viewed from the
                    side of the triangle opposite a reference vertex.
    """
    def __init__(self, Vertex1, Vertex2, Vertex3, VertexR=None, skip_orientation=0, zero=zero, parabaloid=False):


        """DESCRIPTION OF ARGUMENTS:

           Vertex1, Vertex2, and Vertex3: the three vertex objects that comprise the triangle.

           VertexR:

           skip_orientation: 
        """
        self.id = "no id"
        self.v1 = Vertex1
        self.v2 = Vertex2
        self.v3 = Vertex3
        self.vR = VertexR
        self.zero = zero
        self.ltriangle = None
        self.gp = 1

        # Facet id.
        ###########
        self.id = tuple(sorted([self.v1.id, self.v2.id, self.v3.id]))

        # # First, test gp in 2D...
        # if not gp2D(self.v1, self.v2, self.v3):
        #     self.gp = 0

        # i = 0
        # while not self.gp:
        #     print("jostling in Triangle __init__ due to gp2D")
        #     if parabaloid:
        #         self.v3.jostle(parabaloid=True)
        #         print("-----------------------------")
        #         print(self.v1.x, self.v1.y, self.v1.z)
        #         print(self.v2.x, self.v2.y, self.v2.z)
        #         print(self.v3.x, self.v3.y, self.v3.z)
        #     else:
        #         self.v3.jostle()
        #     if gp2D(self.v1, self.v2, self.v3):
        #         self.gp = 1
        #         self.v3.no_more_jostling = True
        #     i += 1
        #     if i > 10:
        #         self.gp = 0
        #         print(self.v1.id, self.v2.id, self.v3.id)
        #         print(self.v1.coordinate_tuple)
        #         print(self.v2.coordinate_tuple)
        #         print(self.v3.coordinate_tuple)
        #         print("Warning, could not achieve general position for triangle...")
        #         break

        # Set reference vertex...
        if not VertexR: 
            if not skip_orientation:
                self.vR = Vertex((0.,0.,0))

        # Test for general position of coordinates in 3D.
        if self.vR:
            if not gp3D(self.v1, self.v2, self.v3, self.vR, zero=self.zero):
                test = gp3D(self.v1, self.v2, self.v3, self.vR, returnVolume=1, zero=self.zero)
                self.gp = 0

            # If necessary, jostle vR coordinates to achieve general position.
            while not self.gp:
                test = gp3D(self.v1, self.v2, self.v3, self.vR, returnVolume=1, zero=self.zero)
                if parabaloid:
                    self.vR.jostle(parabaloid=True)
                else:
                     self.vR.jostle()
                if gp3D(self.v1, self.v2, self.v3, self.vR, zero=self.zero):
                    self.gp = 1 
                    self.vR.no_more_jostling = True

        # Orient the vertices of the triangle relative to the reference vertex vR. 
        # The function orient_vertices() also test for general position.
        ##########################################################################

        if not skip_orientation and self.gp:
            if self.vR:
                self.orient_vertices()

        # Edge attributes. Directed line, CW circuit. 
        #############################################
        self.e1 = Edge(self.v1,self.v2)
        self.e2 = Edge(self.v2,self.v3)
        self.e3 = Edge(self.v3,self.v1)

        # For triangulation calculations, set to True if it is a triangulation facet.
        ############################################################################
        self.bottom = False

         # For triangulation calculations, set to True if it is a surface triangle.
        ############################################################################       
        self.visible = False 
        self.surface = False

    def orient_vertices(self):

        """ The three vertices of the triangle object --
            Vertex1, Vertex2, and Vertex3 -- trace a counterclockwise
            circuit when viewed from the side away from VertexR if the
            sign of the volume of the tetrahedron comprised of these
            four vertices is positive.  By convention, i.e., the
            right-hand rule, this means that VertexR is located on the
            negative side (i.e., below) of the hyperplane (simply a
            plane in dimension 3) formed by Vertex1, Vertex2, and
            Vertex3.  By convention, VertexR should be located on the
            positive side of the hyperplane for calculations.  
        """
        """ Calculate the volume of the tetrahedron assuming that the
            directed edges (Vertex1,Vertex2), (Vertex3,Vertex1), and
            (Vertex3,Vertex2) trace a counterclockwise circuit as
            described above. 
        """

        """Check if the four vertices of the triangle are in general 
           position."""
        sign_of_volume = gp3D(self.v1, self.v2, self.v3, self.vR, returnVolume=1, zero=self.zero)

        if sign_of_volume > 0:
            """ Reversing the orientation of the directed line
                (Vertex1,Vertex2) to (Vertex2,Vertex1) will ensure that
                the volume of the tetrahedron, relative to VertexR, is
                negative.
            """
            swap = self.v1
            self.v1, self.v2 = self.v2, swap

            # Directed line, CW circuit.
            ############################
            self.e1 = Edge(self.v1,self.v2)
            self.e2 = Edge(self.v2,self.v3)
            self.e3 = Edge(self.v3,self.v1)

            # Also missing...?
            ##################
            self.make_ltriangle(zero=self.zero)
            self.set_edge_topology()

    def make_ltriangle(self, zero=zero):

        if self.vR:
            # print(f"Making left triangle ({self.v2.id, self.v1.id, self.vR.id}) {self.v3.id}")
            # print(f"Right {gp3D(self.v1, self.v2, self.v3, self.vR, returnVolume=1)}")
            # print(f"Left  {gp3D(self.v2, self.v1, self.vR, self.v3, returnVolume=1)}")
            self.ltriangle = Triangle(self.v2, self.v1, self.vR, VertexR=self.v3, zero=zero, skip_orientation=1)
        self.ltriangle.ltriangle = self

    def set_edge_topology(self,kill=0):

        # Facet topologies.
        ###################
        self.e1.rf = self

        # Edge topologies.
        ##################
        self.e1.Rprev = self.e2.id
        self.e1.Dnext = self.e2.flipid
        self.e1.Rnext = self.e3.id
        self.e1.Oprev = self.e3.flipid

        if self.ltriangle:

            lt = self.ltriangle

            #Set the remaining topology information for self.
            #################################################
            self.e1.lf = lt

            # Set the topology information for self.ltriangle.
            ##################################################
            if not kill:
                self.ltriangle.set_edge_topology(kill=1)

            # Edge topologies.
            ##################
            self.e1.Lprev = lt.e2.flipid
            self.e1.Onext = lt.e2.id
            self.e1.Dprev = lt.e3.id
            self.e1.Lnext = lt.e3.flipid

    def update_triangle(self, zero=zero):

        # Why not just call self.__init__()? # 3/8/22

        # Update Edge ojbects.
        ######################
        self.e1 = Edge(self.v1,self.v2)
        self.e2 = Edge(self.v2,self.v3)
        self.e3 = Edge(self.v3,self.v1)

        # Update triangle id.
        #####################
        self.id = [self.v1.id, self.v2.id, self.v3.id]
        self.id.sort()
        self.id = tuple(self.id)

        # Make left triangle...XXXX added 2020.05.15
        self.make_ltriangle(zero=self.zero)

        # Update edge topology.
        #######################
        self.set_edge_topology()

    def get_vertices(self):

        return [self.v1, self.v2, self.v3]
        #return {self.v1.id:self.v1, self.v2.id:self.v2, self.v3.id:self.v3}

    def get_edges(self):

        if self.ltriangle:
            return {self.e1.id:self.e1,
                    self.ltriangle.e1.id:self.ltriangle.e1}
        else:
            return {self.e1.id:self.e1}

    def beneath_beyond(self, vertex, returnVolume=0, returnVolumeTuple=0):

        sign_of_volume = gp3D(self.v1, self.v2, self.v3, vertex, returnVolume=1, zero=self.zero)
        absolute_volume = abs(sign_of_volume)

        # Default: Return the absolute volume.
        #####################################
        if returnVolume:
            return absolute_volume #sign_of_volume # XXXX 2017.09.15
        
        # Fixed on 2026.01.18 (see commented-out code below for difference)
        # Return volume tuple.
        ######################
        if returnVolumeTuple:
            if absolute_volume == self.zero:
                return (0, absolute_volume)
            elif sign_of_volume > self.zero:
                return (1, absolute_volume)
            else:
                return (-1, absolute_volume)
        # Return orientation to facet plane.
        ####################################
        if absolute_volume == self.zero:
            # Vertex is not in general position.
            return 0 
        elif sign_of_volume > self.zero: 
            # Vertex outside of facet plane.
            return 1
        elif sign_of_volume < self.zero: 
            # Vertex inside of facet plane.
            return -1

        # # Return volume tuple.
        # ######################
        # if returnVolumeTuple:
        # 	if absolute_volume <= self.zero:
        # 		return (0, absolute_volume)
        # 	elif sign_of_volume > 0:
        # 		return (1, absolute_volume)
        # 	else:
        # 		return (-1, absolute_volume)
        # # Return orientation to facet plane.
        # ####################################
        # if absolute_volume <= self.zero:
        #     # Vertex is not in general position.
        #     return 0 
        # elif sign_of_volume > 0: 
        #     # Vertex outside of facet plane.
        #     return 1
        # elif sign_of_volume < 0: 
        #     # Vertex inside of facet plane.
        #     return -1
        
    def circuit_length(self):
        return (self.e1.length() + self.e2.length() + self.e3.length())
        
    def area(self):
        """The area of the triangle is one-half the area of the 
           parallelogram made by the cross product (or 3x3 determinant)
           of the directed lines v1 --> v2 and v1 --> v3.
        """
        area = cross_product(self.v1,self.v2,self.v3)
        return (1./2.)*abs(area)

    def print_vertex_information(self):
        print(self.v1.x, self.v1.y, self.v1.z)
        print(self.v2.x, self.v2.y, self.v2.z)
        print(self.v3.x, self.v3.y, self.v3.z)

    def __repr__(self): 
        return str(self.id)

    def topological__repr__(self):
        string =  'R triangle edges (CCW circuit): '+ str(self.id)+' : '+ \
                str(self.e1)+' --> '+str(self.e2)+' --> '+str(self.e3)+'\n'
        string += 'R triangle edges (CW  circuit): '+ str(self.id)+' : '+ \
                str(self.e1)+' --> '+str(self.e3)+' --> '+str(self.e2)+'\n'
        string += self.e1.topological__repr__() + '\n'
        lt = self.ltriangle
        if lt:
            string += 'L triangle edges (CCW circuit): '+str(lt)+' : '+str(lt.e1)+ \
                      ' --> '+str(lt.e2)+' --> '+str(lt.e3)+'\n'
            string += 'L triangle edges (CW  circuit): '+str(lt)+' : '+str(lt.e1)+ \
                      ' --> '+str(lt.e3)+' --> '+str(lt.e2)+'\n'
            string += lt.e1.topological__repr__() + '\n'
        else:
            string += '*L triangle has not be intitialized.*\n'
        return string

class Tetrahedron:

    def __init__(self, v1, v2, v3, v4, zero=zero):

        self.v1 = v1
        self.v2 = v2
        self.v3 = v3
        self.v4 = v4
        self.v = (self.v1, self.v2, self.v3, self.v4)
        self.id = tuple(sorted([self.v1.id, self.v2.id, self.v3.id, self.v4.id]))
        self.zero = zero

        self.edges = {}

        # Make the facets of the tetrahedron....
        self.t123 = Triangle(v2,v3,v4,VertexR=v1, zero=self.zero)
        self.t123.make_ltriangle(zero=self.zero); self.t123.set_edge_topology()
        self.edges.update(self.t123.get_edges())
        print(self.t123.get_edges())

        self.t013 = Triangle(v4,v2,v1,VertexR=v3, zero=self.zero)
        self.t013.make_ltriangle(zero=self.zero); self.t013.set_edge_topology()
        self.edges.update(self.t013.get_edges())
        print(self.t013.get_edges())

        self.t012 = Triangle(v1,v2,v3,VertexR=v4, zero=self.zero)
        self.t012.make_ltriangle(zero=self.zero); self.t012.set_edge_topology()
        self.edges.update(self.t012.get_edges())
        print(self.t012.get_edges())

        self.t320 = Triangle(v4,v3,v1,VertexR=v2, zero=self.zero)
        self.t320.make_ltriangle(zero=self.zero); self.t320.set_edge_topology()
        self.edges.update(self.t320.get_edges())
        print(self.t320.get_edges())

        self.t023 = Triangle(v1,v3,v4,VertexR=v2, zero=self.zero)
        self.t023.make_ltriangle(zero=self.zero); self.t023.set_edge_topology()
        self.edges.update(self.t023.get_edges())
        print(self.t023.get_edges())

        self.t230 = Triangle(v1,v4,v3,VertexR=v2, zero=self.zero)
        self.t230.make_ltriangle(zero=self.zero); self.t230.set_edge_topology()
        self.edges.update(self.t230.get_edges())
        print(self.t230.get_edges())


        self.facets = (self.t123, self.t013, self.t012, self.t320, self.t023, self.t230)

    def __repr__(self):
        return str(self.id)

# My library of common computational geometry needs.
#####################################################################################################################

def distance2D(Vertex1,Vertex2):

	#from math import sqrt
	return sqrt((Vertex1.x-Vertex2.x)**2 + (Vertex1.y-Vertex2.y)**2 )

def distance(atom1, atom2, returnSquareDistance=0):

	if returnSquareDistance:
		return (atom1.x-atom2.x)**2 + (atom1.y-atom2.y)**2 + (atom1.z-atom2.z)**2
	else:
		return sqrt((atom1.x-atom2.x)**2 + (atom1.y-atom2.y)**2 + (atom1.z-atom2.z)**2)

def to_vec(p: Any) -> np.ndarray:
    """Convert an object with .x/.y/.z or a sequence/array to a 3D numpy vector."""
    if hasattr(p, "x") and hasattr(p, "y") and hasattr(p, "z"):
        return np.array([p.x, p.y, p.z], dtype=float)
    arr = np.asarray(p, dtype=float)
    if arr.ndim == 1 and arr.size == 3:
        return arr.copy()
    raise TypeError("Point must have .x/.y/.z or be a sequence of length 3")

def cross_product(
    ref: Any,
    a: Any,
    b: Any,
    *,
    return_normal: bool = False,
    normal_length: Optional[float] = None
    ) -> Union[float, Tuple[float, np.ndarray]]:
    """
    Compute the area of the parallelogram defined by vectors (ref - a) and (ref - b).
    Optionally return the normal vector, optionally scaled to `normal_length`.

    Parameters
    ----------
    ref, a, b :
        Points (objects with .x/.y/.z or 3-element sequences)
    return_normal : bool
        If True, also return the normal vector.
    normal_length : float or None
        If provided, scale the unit normal to this length. If None and return_normal is True,
        returns the unit normal.

    Returns
    -------
    area : float
        Magnitude of cross product (parallelogram area).
    (area, normal) : tuple
        If return_normal is True, returns (area, normal_vector). The normal vector is
        unit-length unless normal_length is given, in which case it has that magnitude.
    """
    ref_v = to_vec(ref)
    u = ref_v - to_vec(a)
    v = ref_v - to_vec(b)
    cross_vec = np.cross(u, v)
    area = np.linalg.norm(cross_vec)

    if not return_normal:
        return area

    if area == 0:
        # Degenerate: define normal as zero vector
        normal = np.zeros(3, dtype=float)
    else:
        unit_normal = cross_vec / area
        if normal_length is not None:
            normal = unit_normal * normal_length
        else:
            normal = unit_normal

    return area, normal

def crossProductCoordinates(atomObject1, atomObject2, atomObject3, psa, multiplier=1, norm=0, returnUnitVector=0):

    # Think about folding this into cross_product

    # This cross product function returns the coordinates of the cross product vector defined by vectors u and v. 
    #############################################################################################################

    # The vertex chosen to serve as the reference, or common vertex, incident to both vectors u and v.
    ##################################################################################################
    x0, y0, z0 = atomObject1.x, atomObject1.y, atomObject1.z

    # Vector u pseudo-translated to the origin via the common vertex.
    #################################################################
    u1, u2, u3 = x0-atomObject2.x, y0-atomObject2.y, z0-atomObject2.z

    # Vector v pseudo-translated to the origin via the common vertex.
    #################################################################
    v1, v2, v3 = x0-atomObject3.x, y0-atomObject3.y, z0-atomObject3.z

    # The calculated determinants for u x v.
    # a x b = (a2b3 - a3b2)i + (a3b1 - a1b3)j + (a1b2 - a2b1)k
    ###########################################################
    x = u2*v3 - u3*v2
    y = u3*v1 - u1*v3
    z = u1*v2 - u2*v1

    # The calculated magnitude of the cross product vector || u x v ||.
    ###################################################################
    from math import sqrt
    crossMag = sqrt((x*x)+(y*y)+(z*z))

    # The calculated magnitude of vector | u |. 
    ###########################################
    uMag = sqrt((u1*u1)+(u2*u2)+(u3*u3))

    # The calculated magnitude of vector | v |. 
    ###########################################
    vMag = sqrt((v1*v1)+(v2*v2)+(v3*v3))

    # Build the cross product vector.
    #################################
    # Return the cross product vector as a unit vector.
    ###################################################
    if returnUnitVector:
        psa.x = x/crossMag
        psa.y = y/crossMag
        psa.z = z/crossMag
        psa.v = Vertex4D((psa.x, psa.y, psa.z, (psa.x**2 + psa.y**2 + psa.z**2)))
        psa.v.data = psa
    # Return the normalized cross product vector translated to (x0, y0, z0).
    ########################################################################
    elif norm:
        psa.x = (multiplier*(x/crossMag)) + x0
        psa.y = (multiplier*(y/crossMag)) + y0
        psa.z = (multiplier*(z/crossMag)) + z0
        psa.v = Vertex4D((psa.x, psa.y, psa.z, (psa.x**2 + psa.y**2 + psa.z**2)))
        psa.v.data = psa
    # Return the cross product vector translated to (x0, y0, z0).
    #############################################################
    else:
        psa.x = (multiplier*x) + x0
        psa.y = (multiplier*y) + y0
        psa.z = (multiplier*z) + z0
        psa.v = Vertex4D((psa.x, psa.y, psa.z, (psa.x**2 + psa.y**2 + psa.z**2)))
        psa.v.data = psa

    return psa

# def to_vec(p: Any) -> np.ndarray:
#     """Convert an object with .x/.y/.z or a sequence/array to a 3D numpy vector."""
#     if hasattr(p, "x") and hasattr(p, "y") and hasattr(p, "z"):
#         return np.array([p.x, p.y, p.z], dtype=float)
#     arr = np.asarray(p, dtype=float)
#     if arr.ndim == 1 and arr.size == 3:
#         return arr.copy()
#     raise TypeError("Point must have .x/.y/.z or be a sequence of length 3")

def dot_product(
    ref: Any,
    a: Any,
    b: Any,
    *,
    return_angle: bool = False,
    return_unit_vectors: bool = False,
    projection: Optional[str] = None,  # "u_on_v" or "v_on_u"
    angle_in_degrees: bool = False
    ) -> Union[
    float,
    Tuple[float, float],
    Tuple[float, float, np.ndarray, np.ndarray],
    Tuple[float, Optional[float], Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]
    ]:
    """
    Compute the dot product of (ref - a) and (ref - b), with optional extras.

    Parameters
    ----------
    ref, a, b :
        Points (objects with .x/.y/.z or 3-element sequences)
    return_angle :
        If True, also return the angle between the vectors.
    return_unit_vectors :
        If True, also return the unit vectors of u and v.
    projection :
        If "u_on_v", returns the vector projection of u onto v.
        If "v_on_u", returns the vector projection of v onto u.
    angle_in_degrees :
        If True and return_angle is requested, returns angle in degrees instead of radians.

    Returns
    -------
    Depending on flags:
      - Always returns dot (float).
      - If return_angle: also angle.
      - If return_unit_vectors: also unit_u and unit_v.
      - If projection: also projected vector.
    """
    u = to_vec(ref) - to_vec(a)
    v = to_vec(ref) - to_vec(b)
    norm_u = np.linalg.norm(u)
    norm_v = np.linalg.norm(v)
    dot = float(np.dot(u, v))

    angle = None
    if return_angle:
        if norm_u == 0 or norm_v == 0:
            angle = None  # undefined if one vector is zero
        else:
            # numerical safety for arccos
            cos_theta = dot / (norm_u * norm_v)
            cos_theta = max(-1.0, min(1.0, cos_theta))
            angle = np.arccos(cos_theta)
            if angle_in_degrees:
                angle = np.degrees(angle)

    unit_u = unit_v = None
    if return_unit_vectors:
        unit_u = None if norm_u == 0 else u / norm_u
        unit_v = None if norm_v == 0 else v / norm_v

    proj_vec = None
    if projection is not None:
        if projection == "u_on_v":
            if norm_v == 0:
                proj_vec = None
            else:
                proj_vec = (np.dot(u, v) / (norm_v ** 2)) * v
        elif projection == "v_on_u":
            if norm_u == 0:
                proj_vec = None
            else:
                proj_vec = (np.dot(v, u) / (norm_u ** 2)) * u
        else:
            raise ValueError("projection must be 'u_on_v', 'v_on_u', or None")

    # Build return tuple based on requested pieces
    result = (dot,)
    if return_angle:
        result += (angle,)
    if return_unit_vectors:
        result += (unit_u, unit_v)
    if projection is not None:
        result += (proj_vec,)

    # Unpack single element
    if len(result) == 1:
        return result[0]
    return result

def findMinSidechainDistance(a, b):
    
    # !!!! Re-evaluate need for special treatment of Arg and Lys.
    
    # The parameters a and b are atom distances.
    # This function finds the shortest distance between the atoms of two sidechains.
    # The calculation does not include the backbone atoms of the sidechain.
    #################################################################################
    
    # A pseudo-atom has no residue: convexHull4D seeds the triangulation with a
    # template tetrahedron of PseudoAtom vertices, and with a small enough query
    # set -- the narrow residue selections, acidicSet and the like -- those
    # survive pruning and turn up here as a neighbour. It has no side chain, so
    # it stands for itself, exactly as the HOH/ARG/LYS case below does.
    ################################
    # Define atoms for side chain a.
    ################################
    aAtoms = None
    if a.residue is not None and a.residue.name not in ["HOH", "ARG", "LYS"]:
        aAtoms = a.residue.get_sidechain_atoms()
        # If this is the case, no atoms are associated with the side chain.
        ###################################################################
        if not aAtoms:
            aAtoms = [a,]
    else:
        aAtoms = [a,]

    # Define atoms for side chain b.
    ################################
    bAtoms = None
    if b.residue is not None and b.residue.name not in ["HOH", "ARG", "LYS"]:
        bAtoms = b.residue.get_sidechain_atoms()
        # If this is the case, no atoms are associated with the side chain.
        ###################################################################
        if not bAtoms:
            bAtoms = [b,]
    else:
        bAtoms = [b,]

    # Calculate the minimum distance between the side chain atoms for a and b.
    ##########################################################################
    minDistance = 10000000.
    for aAtom in aAtoms:
        for bAtom in bAtoms:
            
            # Calculate the distance between aAtom and bAtom.
            #################################################
            d = distance(aAtom, bAtom)
            
            # Make fuzzy sidechains due to reduced representations. This "softens" the network calculation.
            # A pseudo-atom has no residue and gets no softening -- it does not
            # stand for a reduced side chain, it is a triangulation seed.
            ###############################################################################################
            if aAtom.residue is not None and aAtom.residue.name in ["HOH", "ARG", "LYS", "CYS"]:
                d -= 0.5
            if bAtom.residue is not None and bAtom.residue.name in ["HOH", "ARG", "LYS", "CYS"]:
                d -= 0.5
            if d < minDistance:
                minDistance = d
            
    return minDistance

def collinear(Vertex1, Vertex2, Vertex3, minArea=zero):

	# No, not collinear.
	####################
	if (1./2.)*abs(cross_product(Vertex1, Vertex2, Vertex3)) >= minArea:
		#print((1./2.)*abs(crossProduct(Vertex1, Vertex2, Vertex3)))
		return 0
	# Yes, collinear and a problem.
	###############################
	else:
		#print((1./2.)*abs(crossProduct(Vertex1, Vertex2, Vertex3)))
		return 1

def centroid3D(vertex3DList,psa):

    n, xSum, ySum, zSum = 0., 0., 0., 0.
    for vertex in vertex3DList:
        n += 1.0
        xSum += vertex.x
        ySum += vertex.y
        zSum += vertex.z

    # Make the Vertex object representing the centroid.
    ###################################################
    result = Vertex((xSum/n,ySum/n,zSum/n))
    result.id = int((xSum/n) * (ySum/n) * (zSum/n))
    psa.atom_serial = int((xSum/n) * (ySum/n) * (zSum/n))
    psa.v = result # Overwrite the original Vertex() object with new Vertex() object.
    psa.x, psa.y, psa.z = xSum/n, ySum/n, zSum/n
    psa.reinitialize()
    result.data = psa

    return result

def centroid4D(vertex4DList, psa):

    n, xSum, ySum, zSum, uSum = 0., 0., 0., 0., 0.
    for vertex in vertex4DList:
        n += 1.0
        xSum += vertex.x
        ySum += vertex.y
        zSum += vertex.z
        uSum += vertex.u

    ###################################################
    result = Vertex4D((xSum/n,ySum/n,zSum/n,uSum/n))
    result.id = abs(int((xSum/n) * (ySum/n) * (zSum/n) * (uSum/n)))
    psa.atom_serial = abs(int((xSum/n) * (ySum/n) * (zSum/n) * (uSum/n)))
    psa.v = result # Overwrite the original Vertex() object with new Vertex4D() object.
    psa.x, psa.y, psa.z = xSum/n, ySum/n, zSum/n
    psa.reinitialize()
    result.data = psa

    return result

def gp2D(v1, v2, v3, zero=zero):
    """Are the three vertices non-collinear?  0 if collinear, 1 if not.

    This is a 3D test, and has to be.  It was written as a 2x2 determinant on x
    and y alone:

        bx, by = x2 - x1, y2 - y1
        det = 2.0 * (bx * cy - by * cx)

    Both terms of that determinant carry a factor of bx or by, so when v1 and v2
    share an x and a y -- a column of points on a vertical line -- det is
    identically zero for every possible v3.  The callers all respond to a zero
    by jostling v3 and testing again, and no movement of v3 can change bx or by,
    so the loop never terminates.  Three points that are plainly non-collinear
    in space, such as (0,0,2), (0,0,2.2) and (35,13,23), were reported collinear.

    Twice the triangle area from the true cross product is the test the original
    code used, and it is zero only when the three points really are collinear.
    """
    area = 0.5 * abs(float(cross_product(v1, v2, v3)))
    # Decimal, as in gp3D and gp4D, for consistent behaviour near zero.
    if Decimal(area) <= zero:
        return 0
    else:
        return 1

def planeCoefficients3D(v1, v2, v3):
    p1 = np.array([v1.x, v1.y, v1.z], dtype=float)
    p2 = np.array([v2.x, v2.y, v2.z], dtype=float)
    p3 = np.array([v3.x, v3.y, v3.z], dtype=float)
    normal = np.cross(p2 - p1, p3 - p1)  # (A, B, C)
    A, B, C = normal
    D = -normal.dot(p1)
    return A, B, C, D

def distanceToPlane(v1, v2, v3, v0):

    c = planeCoefficients3D(v1,v2,v3)
    numerator = c[0]*v0.x + c[1]*v0.y + c[2]*v0.z + c[3]
    denominator = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2])
    return (numerator/denominator)

def gp3D(v1, v2, v3, vT, returnVolume=0, zero=zero):

    # Create the hyperplane for the three vertex4D instances.
    c = planeCoefficients3D(v1,v2,v3)
    # Calculate the volume of the 4-simplex.
    v1 = (c[0]*vT.x + c[1]*vT.y + c[2]*vT.z + c[3])/6.0

    if returnVolume:
        # Return the signed volume of the hypervolume.
        return v1
    elif Decimal(abs(v1)) <= zero: 
        # Not in general position.
        return 0
    else:
        # Convert 3D hyper-volumes to binary score (1 or -1).
        return v1/abs(v1)

# def distance_to_hyperplane4D(v1, v2, v3, v4, v0):
#     A, B, C, D, E = planeCoefficients4D(v1, v2, v3, v4)
#     numerator = A * v0.x + B * v0.y + C * v0.z + D * v0.u + E
#     denom = np.sqrt(A * A + B * B + C * C + D * D)
#     return numerator / denom

def planeCoefficients4D(v1,v2,v3,v4):

    # Think about trying another version of this that avoids det4x4. 

    v1x, v1y, v1z, v1u = v1.x, v1.y, v1.z, v1.u
    v2x, v2y, v2z, v2u = v2.x, v2.y, v2.z, v2.u
    v3x, v3y, v3z, v3u = v3.x, v3.y, v3.z, v3.u
    v4x, v4y, v4z, v4u = v4.x, v4.y, v4.z, v4.u
    A = det4x4([1., 1., 1., 1., v1y, v2y, v3y, v4y, v1z, v2z, v3z, v4z, v1u, v2u, v3u, v4u])
    B = det4x4([v1x, v2x, v3x, v4x, 1., 1., 1., 1., v1z, v2z, v3z, v4z, v1u, v2u, v3u, v4u])
    C = det4x4([v1x, v2x, v3x, v4x, v1y, v2y, v3y, v4y, 1., 1., 1., 1., v1u, v2u, v3u, v4u])
    D = det4x4([v1x, v2x, v3x, v4x, v1y, v2y, v3y, v4y, v1z, v2z, v3z, v4z, 1., 1., 1., 1.])
    E = -det4x4([v1x, v2x, v3x, v4x, v1y, v2y, v3y, v4y, v1z, v2z, v3z, v4z, v1u, v2u, v3u, v4u])
    return (A,B,C,D,E) 

def gp4D(v1,v2,v3,v4,vT, returnVolume=0, zero=zero):

    # Create the hyperplane for the four vertex4D instances.
    c = planeCoefficients4D(v1,v2,v3,v4)
    # Calculate the volume of the 5-simplex.
    v1 = (c[0]*vT.x + c[1]*vT.y + c[2]*vT.z + c[3]*vT.u + c[4])/24.

    if returnVolume:
        # Return the signed volume of the hypervolume.
        return v1
    elif Decimal(abs(v1)) <= zero: # 2019.06.25 Decimal...an attempt to deal with zero inconsistencies on the cluster
        # Not in general position.
        return 0
    else:
        # Convert 4D hyper-volumes to binary score.
        return v1/abs(v1)

# Version to use with faceted surfaces.
#######################################
def circumCircle(facet, psa1, psa2, psa3):
    
    # Have to pass PseudoAtom instances (psa1, psa2, psa3) instead of importing
    # them because other modules (e.g. pdb.py) import this compGeometry module.
    ############################################################################

    # Calculate the centroid of the three atoms.
    ############################################
    v1, v2, v3 = facet.v1, facet.v2, facet.v3
    mp12 = centroid([v1,v2])
    mp13 = centroid([v1,v3])
    mp23 = centroid([v2,v3])

    # fN = crossProductCoordinates(v1,v2,v3,returnUnitVector=1)
    fN = cross_product(v1, v2, v3, return_normal=True)

    # Bisector normal 12.
    #####################
    mp12N = psa1
    mp12N.x = mp12.x + fN.x
    mp12N.y = mp12.y + fN.y
    mp12N.z = mp12.z + fN.z

    # Bisector normal 13.
    #####################
    mp13N = psa2
    mp13N.x = mp13.x + fN.x
    mp13N.y = mp13.y + fN.y
    mp13N.z = mp13.z + fN.z

    # fN12 = crossProductCoordinates(v1,mp12,mp12N,returnUnitVector=1)
    fN12 = cross_product(v1, mp12, mp12N, return_normal=True)
    a, b, c = fN12.x, fN12.y, fN12.z
    # fN13 = crossProductCoordinates(v1,mp13,mp13N,returnUnitVector=1)
    fN13 = cross_product(v1, mp13, mp13N, return_normal=True)
    d, e, f = fN13.x, fN13.y, fN13.z

    x1, y1, z1 = mp12.x, mp12.y, mp12.z
    x2, y2, z2 = mp13.x, mp13.y, mp13.z
    t = ((e*x1)/d - (e*x2)/d - y1 + y2)/(b - (e*a)/d)
    s = ((b*x2)/a - (b*x1)/a - y2 + y1)/(e - (b*d)/a)

    circumCenter = psa3
    circumCenter.x = x1 + a*t
    circumCenter.y = y1 + b*t
    circumCenter.z = z1 + c*t

    r = sqrt((circumCenter.x-v1.x)**2 + \
         (circumCenter.y-v1.y)**2 + \
         (circumCenter.z-v1.z)**2)

    return (r,circumCenter)

# # Version to use with a set of atoms.
# #####################################
# def circumCircleFromAtoms(atom1, atom2, atom3, psa1, psa2, psa3, psa4, psa5, psa6, psa7, psa8, psa9):
    
#     # Have to pass PseudoAtom instances (psa1, psa2, etc.) instead of importing
#     # them because other modules (e.g. pdb.py) import this compGeometry module.
#     ############################################################################

#     # Calculate the centroid of the three atoms.
#     ############################################
#     a1, a2, a3 = atom1, atom2, atom3
#     mp12 = centroid([a1,a2], psa1)
#     mp13 = centroid([a1,a3], psa2)
#     mp23 = centroid([a2,a3], psa3)

#     fN = crossProductCoordinates(a1, a2, a3, psa4, returnUnitVector=1)
    
#     # Bisector normal 12.
#     #####################
#     mp12N = psa5
#     mp12N.x = mp12.x + fN.x
#     mp12N.y = mp12.y + fN.y
#     mp12N.z = mp12.z + fN.z
    
#     # Bisector normal 13.
#     #####################
#     mp13N = psa6
#     mp13N.x = mp13.x + fN.x
#     mp13N.y = mp13.y + fN.y
#     mp13N.z = mp13.z + fN.z
    
#     fN12 = crossProductCoordinates(a1, mp12, mp12N, psa7, returnUnitVector=1)
#     a, b, c = fN12.x, fN12.y, fN12.z
#     fN13 = crossProductCoordinates(a1, mp13, mp13N, psa8, returnUnitVector=1)
#     d, e, f = fN13.x, fN13.y, fN13.z
    
#     x1, y1, z1 = mp12.x, mp12.y, mp12.z
#     x2, y2, z2 = mp13.x, mp13.y, mp13.z
#     t = ((e*x1)/d - (e*x2)/d - y1 + y2)/(b - (e*a)/d)
#     s = ((b*x2)/a - (b*x1)/a - y2 + y1)/(e - (b*d)/a)
    
#     circumCenter = psa9
#     circumCenter.x = x1 + a*t
#     circumCenter.y = y1 + b*t
#     circumCenter.z = z1 + c*t

#     r = sqrt((circumCenter.x-a1.x)**2 + \
#              (circumCenter.y-a1.y)**2 + \
#              (circumCenter.z-a1.z)**2)
             
#     return (r,circumCenter)

def _to_vec3(v):
    if hasattr(v, "x") and hasattr(v, "y") and hasattr(v, "z"):
        return np.array([v.x, v.y, v.z], dtype=float)
    arr = np.asarray(v, dtype=float)
    if arr.shape == (3,):
        return arr.copy()
    raise TypeError("Point must have .x/.y/.z or be a length-3 sequence")

def circumSphere(
    v1, v2, v3, v4, psa,
    *,
    returnSquareDistance: bool = False,
    jitter_attempts: int = 5,
    cond_thresh: float = 1e12,
    jitter_scale: float = 1e-8
    ):
    """
    Compute the circumsphere of four 3D points. Updates psa.x/y/z to center.
    Returns (radius or sqrt-dist, psa). If points are degenerate, attempts jostling.
    """
    for attempt in range(jitter_attempts):
        a = _to_vec3(v1)
        b = _to_vec3(v2)
        c = _to_vec3(v3)
        d = _to_vec3(v4)

        # Build linear system M @ O = rhs
        M = np.vstack([b - a, c - a, d - a])  # shape (3,3)
        rhs = np.array([
            np.dot(b, b) - np.dot(a, a),
            np.dot(c, c) - np.dot(a, a),
            np.dot(d, d) - np.dot(a, a),
        ], dtype=float) * 0.5

        # Check conditioning
        cond = np.linalg.cond(M)
        if cond < cond_thresh:
            try:
                center = np.linalg.solve(M, rhs)
            except np.linalg.LinAlgError:
                center = None
        else:
            center = None

        if center is not None:
            diff = center - a
            sqr_rad = np.dot(diff, diff)
            r = sqr_rad if returnSquareDistance else np.sqrt(sqr_rad)
            psa.x, psa.y, psa.z = center.tolist()
            return r, psa

        # if we get here, degeneracy / ill-conditioning: jostle the vertices slightly
        for v in (v1, v2, v3, v4):
            if hasattr(v, "jostle"):
                v.jostle()
            else:
                # fallback: add tiny random perturbation if no jostle method
                if hasattr(v, "x") and hasattr(v, "y") and hasattr(v, "z"):
                    if JOSTLE_DETERMINISTIC:
                        rng = _jostle_rng(v.x, v.y, v.z, attempt + 1)
                        perturb = [(rng.random() * 2 - 1) * jitter_scale
                                   for _ in range(3)]
                    else:
                        perturb = (np.random.uniform(-1, 1, size=3) * jitter_scale)
                    v.x += perturb[0]
                    v.y += perturb[1]
                    v.z += perturb[2]

    raise ValueError("Failed to compute non-degenerate circumsphere after jostling")

# def circumSphere(v1,v2,v3,v4,psa,returnSquareDistance=0):

#     # Ensure general position.
#     ##########################
#     #while testGp4D(v1,v2,v3,v4):
#     while not gp3D(v1,v2,v3,v4):
#         v1.jostle()
#         v2.jostle()
#         v3.jostle()
#         v4.jostle()
#         #print("jostling from circumSphere...")

#     x1, y1, z1 = v1.x, v1.y, v1.z
#     x2, y2, z2 = v2.x, v2.y, v2.z
#     x3, y3, z3 = v3.x, v3.y, v3.z
#     x4, y4, z4 = v4.x, v4.y, v4.z
#     m11 = det4x4([x1,y1,z1,1.,x2,y2,z2,1.,x3,y3,z3,1.,x4,y4,z4,1.])
#     m12 = det4x4([(x1*x1 + y1*y1 + z1*z1),y1,z1,1.,(x2*x2 + y2*y2 + z2*z2),y2,z2,1.,
#                (x3*x3 + y3*y3 + z3*z3),y3,z3,1.,(x4*x4 + y4*y4 + z4*z4),y4,z4,1.])

#     m13 = det4x4([(x1*x1 + y1*y1 + z1*z1),x1,z1,1.,(x2*x2 + y2*y2 + z2*z2),x2,z2,1.,
#                (x3*x3 + y3*y3 + z3*z3),x3,z3,1.,(x4*x4 + y4*y4 + z4*z4),x4,z4,1.])

#     m14 = det4x4([(x1*x1 + y1*y1 + z1*z1),x1,y1,1.,(x2*x2 + y2*y2 + z2*z2),x2,y2,1.,
#               (x3*x3 + y3*y3 + z3*z3),x3,y3,1.,(x4*x4 + y4*y4 + z4*z4),x4,y4,1.])

#     m15 = det4x4([(x1*x1 + y1*y1 + z1*z1),x1,y1,z1, (x2*x2 + y2*y2 + z2*z2),x2,y2,z2,
#               (x3*x3 + y3*y3 + z3*z3),x3,y3,z3,(x4*x4 + y4*y4 + z4*z4),x4,y4,z4])

#     x0 =  0.5*(m12/m11)
#     y0 = -0.5*(m13/m11)
#     z0 =  0.5*(m14/m11)
#     origin = psa
#     origin.x, origin.y, origin.z = x0, y0, z0

#     r = 0.
#     if returnSquareDistance:
#         r = x0**2 + y0**2 + z0**2 - (m15/m11)
#     else:
#         r = sqrt(x0**2 + y0**2 + z0**2 - (m15/m11))

#     return (r,origin)

# # Some functions to follow may not working properly...
# ####################################################################################################
# def translate(refAtom, atomListCentroid, atomList):

#     x0, y0, z0 = refAtom.x, refAtom.y, refAtom.z
#     x1, y1, z1 = atomListCentroid.x, atomListCentroid.y, atomListCentroid.z
#     transX, transY, transZ = (x0 - x1), (y0 - y1), (z0 - z1)
#     for atom in atomList:
#         atom.x = transX + atom.x
#         atom.y = transY + atom.y
#         atom.z = transZ + atom.z

# class Rotations:

#     # Adapted from:
#     # https://www.meccanismocomplesso.org/en/3d-rotations-and-euler-angles-in-python/

#     def __init__(self):

#         import numpy as np
#         import math as m
#         # Default values...
#         self.phi = m.pi/4
#         self.theta = m.pi/4
#         self.psi = m.pi/4
#         self.poses = {}
     
#     def Rx(self, rads):
#         import numpy as np
#         import math as m
#         return np.matrix([[ 1, 0           , 0           ],
#                             [ 0, m.cos(rads),-m.sin(rads)],
#                             [ 0, m.sin(rads), m.cos(rads)]])
     
#     def Ry(self, rads):
#         import numpy as np
#         import math as m
#         return np.matrix([[ m.cos(rads), 0, m.sin(rads)],
#                             [ 0           , 1, 0           ],
#                             [-m.sin(rads), 0, m.cos(rads)]])
     
#     def Rz(self, rads):
#         import numpy as np
#         import math as m
#         return np.matrix([[ m.cos(rads), -m.sin(rads), 0 ],
#                             [ m.sin(rads), m.cos(rads) , 0 ],
#                             [ 0           , 0            , 1 ]])

#     def rotationMatrix(self, psi, theta, phi):
#         import numpy as np
#         import math as m
#         return self.Rx(psi) * self.Ry(theta) * self.Rz(phi)

#     def rotationAtomSet(self, atomSet, psi, theta, phi):
#         import numpy as np
#         import math as m
#         from pdbFile import Atom
#         atomSet, pdbString = atomSet[:], "" # Don't alter atomSet in main program...
#         for atom in atomSet:
#             v1 = np.array([[atom.x], [atom.y], [atom.z]])
#             R = self.rotationMatrix(psi, theta, phi)
#             v2 = R * v1
#             atom.x, atom.y, atom.z = float(v2[0]), float(v2[1]), float(v2[2])
#             atom.reinitialize()
#             pdbString += str(atom)
#         key = (psi, theta, phi)
#         self.poses[key] = pdbString

#     def rotationAtomSetPoses(self, atomSet, psi_start, theta_start, phi_start, cut, stop=None):
#         import numpy as np
#         import math as m
#         self.poses = {} # Reset poses...
#         stop = m.pi
#         psi_step, theta_step, phi_step = (stop-psi_start)/cut, (stop-theta_start)/cut, (stop-phi_start)/cut
#         while psi_start <= stop:
#             i = theta_start
#             while i <= stop:
#                 i += theta_step
#                 j = phi_start
#                 k = 0
#                 while j <= stop:
#                     self.rotationAtomSet(atomSet, psi_start, i, j)
#                     j += phi_step
#             psi_start += psi_step

# def dotProduct(atomObject1, atomObject2, atomObject3, angle=0, skipGp=0):

# 	# This function accommodates both vertex and atom objects.
# 	# This accommodates test of general position.
# 	# All geometry-related functions should posses this generality.
# 	###############################################################
# 	try:
# 		c1 = atomObject1.v
# 	except:
# 		c1 = atomObject1
# 	try:
# 		c2 = atomObject2.v
# 	except:
# 		c2 = atomObject2
# 	try:
# 		c3 = atomObject3.v
# 	except:
# 		c3 = atomObject3

# 	# Ensure general position.
# 	##########################
# 	if not skipGp:
# 		while not gp2D(c1, c2, c3):
# 			print("Warning: had to jostle to general position in dotProduct().")
# 			c1.jostle()
# 			c2.jostle()
# 			c3.jostle()

# 	# The vertex chosen to serve as the reference, or common vertex, incident to both vectors u and v.
# 	##################################################################################################
# 	x0, y0, z0 = c1.x, c1.y, c1.z

# 	# Vector u pseudo-translated to the origin via the common vertex.
# 	#################################################################
# 	u1, u2, u3 = x0-c2.x, y0-c2.y, z0-c2.z

# 	# Vector v pseudo-translated to the origin via the common vertex.
# 	#################################################################
# 	v1, v2, v3 = x0-c3.x, y0-c3.y, z0-c3.z

# 	"""The order of the atomObjects matters in this function.
# 	   The first argument (atomObject0) is the reference point.
# 	   Given this reference point, there are two directed lines (i.e., vectors):
	   
# 	                u = atomObject0 ---> atomObject1
# 	                v = atomObject0 ---> atomObject2
	                
# 	    The equation for the dot product of two vectors is,
	    
# 	                u.v = |u||v|cos(theta)
	                
# 	    where, 
	                
# 	                |u| is the magnitude of the vector u
# 	                |v| is the magnitude of the vector v
# 	                theta is the angle between the two vectors.
	    
# 	    The result of a dot product is a scalar that is equal to the
# 	    area of the parallelogram made by the two vectors. 
	    
# 	    The dot product is cummutative: u.v = v.u
# 	    Therefore, the order of the arguments atomObject1 and atomObject2
# 	    does not matter. 
	    
# 	    The optional (angle) argument allows the user to return the angle
# 	    between the two vectors, instead of the scalar value of the dot product.
# 	    A value of 0, returns the scalar value of the dot product. A value of 1,
# 	    returns the angle between the two vectors in radians.
# 	 """

# 	# The dot product.
# 	##################
# 	uv = u1*v1 + u2*v2 + u3*v3

# 	# Return the scalar dot product.
# 	################################
# 	if not angle:
# 		return uv
# 	# Return the angle.
# 	###################
# 	else:
# 		uMag = sqrt((u1*u1) + (u2*u2) + (u3*u3))
# 		vMag = sqrt((v1*v1) + (v2*v2) + (v3*v3))
# 		from math import degrees
# 		return degrees(acos(uv/(uMag*vMag)))

# def crossProduct(atomObject1, atomObject2, atomObject3):

# 	# This cross product function returns the area of the parellogram defined by vectors u and v. 
# 	##############################################################################################

# 	# The vertex chosen to serve as the reference, or common vertex, incident to both vectors u and v.
# 	##################################################################################################
# 	x0, y0, z0 = atomObject1.x, atomObject1.y, atomObject1.z

# 	# Vector u pseudo-translated to the origin via the common vertex.
# 	#################################################################
# 	u1, u2, u3 = x0-atomObject2.x, y0-atomObject2.y, z0-atomObject2.z

# 	# Vector v pseudo-translated to the origin via the common vertex.
# 	#################################################################
# 	v1, v2, v3 = x0-atomObject3.x, y0-atomObject3.y, z0-atomObject3.z

# 	# The calculated determinants for u x v.
# 	# a x b = (a2b3 - a3b2)i + (a3b1 - a1b3)j + (a1b2 - a2b1)k
# 	###########################################################
# 	x = u2*v3 - u3*v2
# 	y = u3*v1 - u1*v3
# 	z = u1*v2 - u2*v1

# 	# Return the parrallelogram area: the magintude of the cross product vector. 
# 	############################################################################
# 	from math import sqrt
# 	return sqrt((x*x)+(y*y)+(z*z))


