def det2x2(det):
	return det[0][0]*det[1][1] - det[0][1]*det[1][0]

def det3x3(d):
	"""The sign chart for a 3x3 determinant is

		| + - + |                                 | 00 01 02 |
		| - + - |                                 | 10 11 12 |
		| + - + |, which corresponds to positions |.20 21 22 |.

		Here, the 3x3 determinant is solved by expanding it
		into its 3 minor 2x2 determinants in 20, 21, 22, i.e.,

		minor20 = 20*det((01,02),(11,12)),
		minor21 = 21*det((00,02),(10,12)),
		minor22 = 22*det((00,01),(10,11)).

		By the sign chart above,the signs of minor20 and minor22
		are unchanged.  The sign of minor21 should be inverted
		by multiplying it by -1. 

		"""

	"""The code runs faster if the 2x2 determinant is calculated
		explicitly here in this function, rather than repeatly calling
		the function det2x2().
	"""
	#minor20 = det[0][1]*det[1][2] - det[0][2]*det[1][1]
	#minor21 = det[0][0]*det[1][2] - det[0][2]*det[1][0]
	#minor22 = det[0][0]*det[1][1] - det[1][0]*det[0][1]
	#d = (det[0][0], det[0][1], det[0][2],
	#     det[1][0], det[1][1], det[1][2],
	#     det[2][0], det[2][1], det[2][2])
	# [0][0] -- 0 
	# [0][1] -- 1 
	# [0][2] -- 2 
	# [1][0] -- 3
	# [1][1] -- 4
	# [1][2] -- 5
	# [2][0] -- 6
	# [2][1] -- 7
	# [2][2] -- 8
	#minor20 = d[1]*d[5] - d[2]*d[4]
	#minor21 = d[0]*d[5] - d[2]*d[3]
	#minor22 = d[0]*d[4] - d[3]*d[1]
	return d[6]*(d[1]*d[5] - d[2]*d[4]) - d[7]*(d[0]*d[5] - d[2]*d[3]) + d[8]*(d[0]*d[4] - d[3]*d[1])


	#PPPP: Reduce the number of lookups.
	#minor20 = det[0][1]*det[1][2] - det[0][2]*det[1][1]
	#minor21 = det[0][0]*det[1][2] - det[0][2]*det[1][0]
	#minor22 = det[0][0]*det[1][1] - det[1][0]*det[0][1]
	#r2 = det[2][0]*minor20 - det[2][1]*minor21 + det[2][2]*minor22
	#print(r1, r2)
	#return det[2][0]*minor20 - det[2][1]*minor21 + det[2][2]*minor22

#def cDet(p1,p2,p3,p4,n):
#
#    # c pointers to doubles: p1, p2, p3, p4
#    aP = type(p1) * 4
#    d = aP(p1,p2,p3,p4)
#    libDet.Determinant.restype = c_double
#    return libDet.Determinant(d,c_int(n))

def det4x4Old1(det):

	"""The sign chart for a 4x4 determinant is

		| + - + - |                                 | 00 01 02 03 |
		| - + - + |                                 | 10 11 12 13 |
		| + - + - |                                 | 20 21 22 23 |
		| - + - + |, which corresponds to positions | 30 31 32 33 |.

		Here, the 4x4 determinant is solved by first deconstructing it
		into its four minor 3x3 determinants in 30, 31, 32, and 33, i.e.,

						| 01 02 03 |
		minor30 = -30*det| 11 12 13 |
						| 21 22 23 |

						| 00 02 03 |
		minor31 =  31*det| 10 12 13 |
						| 20 22 23 |

						| 00 01 03 |
		minor32 = -32*det| 10 11 13 |
						| 20 21 23 |

						| 00 01 02 |
		minor33 =  33*det| 10 11 12 |
						| 20 21 22 |


		By the sign chart above,the signs of minor31 and minor33
		are positive, and the signs of minor30 and minor32 are negative.
		Once the 4x4 determinant has been deconstructed into the four 3x3
		determinants, these 3x3 determinants are deconstructed into
		three minor 2x2 determinants.  For example, for minor30 these
		determinants are

							| 02 03 |      
		minor3020 = 21*det| 12 13 |

							| 01 03 |      
		minor3021 = 22*det| 11 13 |

							| 01 02 |      
		minor3022 = 23*det| 11 12 |.

		To greatly reduce the number of lookups, the coordinates of the
		determinant are explicitly assigned to a variable.  This greatly
		improves the efficiency of this function.
	"""

	x00, y01, z02, a03 = det[0][0], det[0][1], det[0][2], det[0][3]
	x10, y11, z12, a13 = det[1][0], det[1][1], det[1][2], det[1][3]
	x20, y21, z22, a23 = det[2][0], det[2][1], det[2][2], det[2][3]
	x30, y31, z32, a33 = det[3][0], det[3][1], det[3][2], det[3][3]

	#minor30
	m3020 = z02*a13 - a03*z12
	m3021 = y01*a13 - a03*y11
	m3022 = y01*z12 - y11*z02
	minor30 = y21*m3020 - z22*m3021 + a23*m3022

	#minor31
	m3120 = z02*a13 - a03*z12
	m3121 = x00*a13 - a03*x10
	m3122 = x00*z12 - x10*z02
	minor31 = x20*m3120 - z22*m3121 + a23*m3122

	#minor32
	m3220 = y01*a13 - a03*y11
	m3221 = x00*a13 - a03*x10
	m3222 = x00*y11 - x10*y01
	minor32 = x20*m3220 - y21*m3221 + a23*m3222

	#minor33
	m3320 = y01*z12 - z02*y11
	m3321 = x00*z12 - z02*x10
	m3322 = x00*y11 - x10*y01
	minor33 = x20*m3320 - y21*m3321 + z22*m3322

	return -x30*minor30 + y31*minor31 - z32*minor32 + a33*minor33


def det4x4Old2(det):

	"""The sign chart for a 4x4 determinant is

		| + - + - |                                 | 00 01 02 03 |
		| - + - + |                                 | 10 11 12 13 |
		| + - + - |                                 | 20 21 22 23 |
		| - + - + |, which corresponds to positions | 30 31 32 33 |.

		Here, the 4x4 determinant is solved by first deconstructing it
		into its four minor 3x3 determinants in 30, 31, 32, and 33, i.e.,

						| 01 02 03 |
		minor30 = -30*det| 11 12 13 |
						| 21 22 23 |

						| 00 02 03 |
		minor31 =  31*det| 10 12 13 |
						| 20 22 23 |

						| 00 01 03 |
		minor32 = -32*det| 10 11 13 |
						| 20 21 23 |

						| 00 01 02 |
		minor33 =  33*det| 10 11 12 |
						| 20 21 22 |


		By the sign chart above,the signs of minor31 and minor33
		are positive, and the signs of minor30 and minor32 are negative.
		Once the 4x4 determinant has been deconstructed into the four 3x3
		determinants, these 3x3 determinants are deconstructed into
		three minor 2x2 determinants.  For example, for minor30 these
		determinants are

							| 02 03 |      
		minor3020 = 21*det| 12 13 |

							| 01 03 |      
		minor3021 = 22*det| 11 13 |

							| 01 02 |      
		minor3022 = 23*det| 11 12 |.

		To greatly reduce the number of lookups, the coordinates of the
		determinant are explicitly assigned to a variable.  This greatly
		improves the efficiency of this function.
	"""

	p0 = det[0]
	x00, y01, z02, a03 = p0[0], p0[1], p0[2], p0[3]
	p1 = det[1]
	x10, y11, z12, a13 = p1[0], p1[1], p1[2], p1[3]
	p2 = det[2]
	x20, y21, z22, a23 = p2[0], p2[1], p2[2], p2[3]
	p3 = det[3]
	x30, y31, z32, a33 = p3[0], p3[1], p3[2], p3[3]

	#minor30
	m3020 = z02*a13 - a03*z12
	m3021 = y01*a13 - a03*y11
	m3022 = y01*z12 - y11*z02
	minor30 = y21*m3020 - z22*m3021 + a23*m3022

	#minor31
	m3120 = z02*a13 - a03*z12
	m3121 = x00*a13 - a03*x10
	m3122 = x00*z12 - x10*z02
	minor31 = x20*m3120 - z22*m3121 + a23*m3122

	#minor32
	m3220 = y01*a13 - a03*y11
	m3221 = x00*a13 - a03*x10
	m3222 = x00*y11 - x10*y01
	minor32 = x20*m3220 - y21*m3221 + a23*m3222

	#minor33
	m3320 = y01*z12 - z02*y11
	m3321 = x00*z12 - z02*x10
	m3322 = x00*y11 - x10*y01
	minor33 = x20*m3320 - y21*m3321 + z22*m3322

	return -x30*minor30 + y31*minor31 - z32*minor32 + a33*minor33

def det4x4(d):

	"""The sign chart for a 4x4 determinant is

		| + - + - |                                 | 00 01 02 03 |
		| - + - + |                                 | 10 11 12 13 |
		| + - + - |                                 | 20 21 22 23 |
		| - + - + |, which corresponds to positions | 30 31 32 33 |.

		Here, the 4x4 determinant is solved by first deconstructing it
		into its four minor 3x3 determinants in 30, 31, 32, and 33, i.e.,

						| 01 02 03 |
		minor30 = -30*det| 11 12 13 |
						| 21 22 23 |

						| 00 02 03 |
		minor31 =  31*det| 10 12 13 |
						| 20 22 23 |

						| 00 01 03 |
		minor32 = -32*det| 10 11 13 |
						| 20 21 23 |

						| 00 01 02 |
		minor33 =  33*det| 10 11 12 |
						| 20 21 22 |


		By the sign chart above,the signs of minor31 and minor33
		are positive, and the signs of minor30 and minor32 are negative.
		Once the 4x4 determinant has been deconstructed into the four 3x3
		determinants, these 3x3 determinants are deconstructed into
		three minor 2x2 determinants.  For example, for minor30 these
		determinants are

							| 02 03 |      
		minor3020 = 21*det| 12 13 |

							| 01 03 |      
		minor3021 = 22*det| 11 13 |

							| 01 02 |      
		minor3022 = 23*det| 11 12 |.

		To greatly reduce the number of lookups, the coordinates of the
		determinant are explicitly assigned to a variable.  This greatly
		improves the efficiency of this function.
	"""

		#p0 = det[0]
	#x00, y01, z02, a03 = p0[0], p0[1], p0[2], p0[3]
		#p1 = det[1]
	#x10, y11, z12, a13 = p1[0], p1[1], p1[2], p1[3]
		#p2 = det[2]
	#x20, y21, z22, a23 = p2[0], p2[1], p2[2], p2[3]
		#p3 = det[3]
	#x30, y31, z32, a33 = p3[0], p3[1], p3[2], p3[3]

		#d = det[0] + det[1] + det[2] + det[3]
		#print(d)
	x00, y01, z02, a03 = d[0], d[1], d[2], d[3]
		#p1 = det[1]
	x10, y11, z12, a13 = d[4], d[5], d[6], d[7]
		#p2 = det[2]
	x20, y21, z22, a23 = d[8], d[9], d[10], d[11]
		#p3 = det[3]
	x30, y31, z32, a33 = d[12], d[13], d[14], d[15]


	return -x30*(y21*(z02*a13 - a03*z12) - \
				z22*(y01*a13 - a03*y11) + \
				a23*(y01*z12 - y11*z02))     + \
			y31*(x20*(z02*a13 - a03*z12) - \
				z22*(x00*a13 - a03*x10) + \
				a23*(x00*z12 - x10*z02))     + \
			-z32*(x20*(y01*a13 - a03*y11) - \
				y21*(x00*a13 - a03*x10) + \
				a23*(x00*y11 - x10*y01))     + \
			a33*(x20*(y01*z12 - z02*y11) - \
				y21*(x00*z12 - z02*x10) + \
				z22*(x00*y11 - x10*y01))


#def testDet(v1,v2,v3,v4):
#
#    # c pointers to doubles: p1, p2, p3, p4
#    aO = type(v1.cx) * 16
#    a = aO(v1.cx,v1.cy,v1.cz,v1.cu,
#           v2.cx,v2.cy,v2.cz,v2.cu,
#           v3.cx,v3.cy,v3.cz,v3.cu,
#           v4.cx,v4.cy,v4.cz,v4.cu)
#    aP = pointer(a)
#    libDet.Determinant4x4.restype = c_double
#    return libDet.Determinant4x4(aP)
