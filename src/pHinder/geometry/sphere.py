""" Spherical polar coordinates: a system of curvilinear coordinates.

    x = r cos(theta) sin(phi)
    y = r sin(theta) sin(phi)
    z = r cos(phi)

    Theta is defined as the azimuthal angle in the xy-plane from the
    x-axis with 0 <= theta < 2pi.
    
    Phi is defined as the polar angle from the z-axis with 0 <= phi < pi.
    
"""
from pHinder._vendor.compGeometry import Vertex
from math import pi

class Sphere:

    def __init__(self):

        #The basic attributes of a sphere
        self.radius = 1.
        self.origin_x = 0.
        self.origin_y = 0.
        self.origin_z = 0.

        #Class attributes for generating and storing
        #coordinates located on the surface of the sphere.
        self.theta_min = 0
        self.theta_max = 2*pi
        self.theta_increment = pi/2.
        self.phi_min = 0
        self.phi_max = pi
        self.phi_increment = pi/2.
        self.polar_coordinates     = []
        self.cartesian_coordinates = []

        #AN ATTRIBUTE FOR ASSOCIATING ADDITIONAL DATA WITH THE CLASS.
        self.xtra_info = None

    def generate_surface_points(self, theta_increment=pi/2., phi_increment=pi/2., general_position=1):

        """ Arguments
            theta_increment: the angle increment in the xy-plane in radians.
            phi_increment: the angle increment in the xz-plane in radians.
            general_position: if true, generate the coordinates of the sphere
                              such that no four coordinates are collinear. 
                              This ensures that the point set is in general
                              position for convex hull calculations. As a 
                              consequence of demanding that the vertices of the
                              sphere be in general position, the total number 
                              of sphere vertices will vary (for a given set of
                              calculation parameters) from run-to-run.
        """

        self.theta_increment = theta_increment
        self.phi_increment = phi_increment

        from math import pi, cos, sin
        from pHinder._vendor.pdbFile import PseudoAtom
        
        phi, jostle_cuts = self.phi_min, [10., 20., 30., 40., 50., 60., 70., 80., 90.]
        while phi <= self.phi_max:
            theta = self.theta_min
            while theta < self.theta_max:
                if theta == 0 and (phi == 0 or phi == pi):
                    x = self.origin_x + self.radius*cos(0)*sin(phi)
                    y = self.origin_y + self.radius*sin(0)*sin(phi)
                    z = self.origin_z + self.radius*cos(phi)
                    if general_position:
                        v = Vertex((x,y,z))
                        psa = PseudoAtom()
                        v.data = psa
                        v.jostle()
                        self.cartesian_coordinates.append((v.x,v.y,v.z))
                    else:
                        self.cartesian_coordinates.append((x,y,z))
                    break
                else:
                    x = self.origin_x + self.radius*cos(theta)*sin(phi)
                    y = self.origin_y + self.radius*sin(theta)*sin(phi)
                    z = self.origin_z + self.radius*cos(phi)
                    if general_position:
                        v = Vertex((x,y,z))
                        psa = PseudoAtom()
                        v.data = psa
                        v.jostle()
                        self.cartesian_coordinates.append((v.x,v.y,v.z))
                    else:
                        self.cartesian_coordinates.append((x,y,z))
                    theta += theta_increment
            phi += phi_increment

#    def generateRandomSurfacePoints(self,sampling,nPoints):
#      
#        # Calculate theta incremental steps
#        thetaStep = (self.theta_max - self.theta_min)/float(sampling)
#        theta, thetaValues = 0, []
#        while theta < self.theta_max:
#            thetaValues.append(theta)
#            theta += thetaStep
#
#        # Calculate phi incremental steps
#        phiStep = (self.phi_max - self.phi_min)/float(sampling)
#        phi, phiValues = 0, []
#        while phi < self.phi_max:
#            phiValues.append(phi)
#            phi += phiStep
#
#        # Generate random points on the surface of a sphere
#        sPoints = 0
#        while sPoints <= nPoints:
#            from random import choice
#            tC, pC = choice(thetaValues), choice(phiValues)
#            from math import cos, sin
#            x = self.origin_x + self.radius*cos(tC)*sin(pC)
#            y = self.origin_y + self.radius*sin(tC)*sin(pC)
#            z = self.origin_z + self.radius*cos(pC)
#            from random import random
#            x += random()/10.
#            y += random()/10.
#            z += random()/10.
#            self.cartesian_coordinates.append((x,y,z))
#            sPoints += 1

    def get_origin_psa(self):
        
        from pHinder._vendor.pdbFile import PseudoAtom
        psa = PseudoAtom()
        psa.x = self.origin_x
        psa.y = self.origin_y
        psa.z = self.origin_z
        psa.v.x = psa.x
        psa.v.y = psa.y
        psa.v.z = psa.z
        psa.v.data = psa
        psa.v.reinitialize()
    
        return psa
        
    def get_pseudoatoms(self, atom_serial=1, return_vertices=False):

        from pHinder._vendor.pdbFile import PseudoAtom
        toReturn = []
        for coordinate in self.cartesian_coordinates:
            psa = PseudoAtom()
            psa.x = coordinate[0]
            psa.y = coordinate[1]
            psa.z = coordinate[2]
            psa.atom_serial = atom_serial + 10000
            psa.residue_sequence_number = atom_serial
            psa.chain_identifier = "Z"
            psa.reinitialize()
            toReturn.append(psa)
            atom_serial += 1
        
        if not return_vertices:
            return toReturn
        else:
            vertices = []
            for psa in toReturn:
                vertices.append(psa.v)
            return(vertices)

    def intersection_with_line(self,v1,v2):

        #Using the convention Paul Bourke

        #P1
        x1, y1, z1 = v1.x, v1.y, v1.z
        #P2
        x2, y2, z2 = v2.x, v2.y, v2.z
        #P3
        x3, y3, z3 = self.origin_x, self.origin_y, self.origin_z
        #a
        a = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2
        #b
        b = 2*(((x2-x1)*(x1-x3)) + ((y2-y1)*(y1-y3)) + ((z2-z1)*(z1-z3)))
        #c
        c = x3**2 + y3**2 + z3**2 + x1**2 + y1**2 + z1**2 - \
            2*(x3*x1 + y3*y1 + z3*z1) - self.radius**2

        #Result
        result = b*b - (4*a*c)
        #If the result is negative, the line does not intersect the sphere.
        if result < 0:
            return 0
        #If the result is zero, the line is tangent to the sphere.
        #if the result is postive, the line intersects the sphere.
        else:
            return 1

    def intersection_with_line_if_between(self,p1,p2):

        """ Like intersection, but tests whether the sphere is
            between v1 and v2
        """
        #Get ray topology.
        from pHinder._vendor.determinants import det3x3
        matrix_p1_p2 = (((p1.x,p1.y,p1.z),
                         (p2.x,p2.y,p2.z),
                         (1.,1.,1.)))
        det_p1_p2 = det3x3(matrix_p1_p2)
        #If det_ab is negative, p1 is left of p2

        p3 = self.get_origin_psa()
            
        #Get ray -- sphere topology.
        matrix_p1_p3 = (((p1.x,p1.y,p1.z),
                         (p3.x,p3.y,p3.z),
                         (1.,1.,1.)))
        det_p1_p3 = det3x3(matrix_p1_p3)
        matrix_p2_p3 = (((p2.x,p2.y,p2.z),
                         (p3.x,p3.y,p3.z),
                         (1.,1.,1.)))
        det_p2_p3 = det3x3(matrix_p2_p3)

        go_p1 = 0
        go_p2 = 0
        #p1 is left of p2
        if det_p1_p2 < 0:
            #p1 is also left of p3
            if det_p1_p3 < 0:
                go_p1 = 1
            #p2 is to the right of p3
            if det_p2_p3 > 0:
                go_p2 = 1
        #p1 is right of p2
        else:
            #p1 is also right of p3
            if det_p1_p3 > 0:
                go_p1 = 1
            #p2 is to the left of p3
            if det_p2_p3 < 0:
                go_p2 = 1

        #!!!!!!This appears to be working well
        if not go_p1 or not go_p2:
            #Should register the same as no intersection
            return 0        
    
        #Using the convention Paul Bourke

        #P1
        x1, y1, z1 = p1.x, p1.y, p1.z
        #P2
        x2, y2, z2 = p2.x, p2.y, p2.z
        #P3
        x3, y3, z3 = self.origin_x, self.origin_y, self.origin_z
        #a
        a = (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2
        #b
        b = 2*(((x2-x1)*(x1-x3)) + ((y2-y1)*(y1-y3)) + ((z2-z1)*(z1-z3)))
        #c
        c = x3**2 + y3**2 + z3**2 + x1**2 + y1**2 + z1**2 - \
            2*(x3*x1 + y3*y1 + z3*z1) - self.radius**2

        #Result
        result = b*b - (4*a*c)
        #If the result is negative, the line does not intersect the sphere.
        if result < 0:
            return 0
        #If the result is zero, the line is tangent to the sphere.
        #if the result is postive, the line intersects the sphere.
        else:
            return 1

#Script code for testing this class
#Usage: python sphere.py [output_directory]   (default: ./spheres)
if __name__ == "__main__" :
    import sys
    from math import pi
    from pathlib import Path

    sphere = Sphere()
    sphere.generate_surface_points(pi/8,pi/8)
    psas = sphere.get_pseudoatoms()

    outpath = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("spheres")
    outpath.mkdir(parents=True, exist_ok=True)

    with open(outpath / "test.pdb", 'w') as output:
        for psa in psas:
            output.write(psa.__repr__())
    print("wrote", outpath / "test.pdb")
