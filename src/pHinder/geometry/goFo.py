# Dependencies.
###############
from pHinder._vendor.compGeometry import distance

def goFo(nStart, been, s2s, s1Dict, skipS2=None, psa=1, level=1, byEdge=0):
    
    # Recursive function that moves through a triangulation to find minimum network paths.
    # It starts at node (nStart), and is forced to move forward (cannot re-visit nStart).
    ######################################################################################
    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue_name == "PSA" and not psa:
        return 0

    s2s = sorted(s2s)
    
    # Record progress to force forward moves.
    #########################################
    if byEdge:
        if nStart == s2s[0][0]:
            been[nStart] = (level, nStart, s2s[0][1])
            #been[s2s[0][1]] = (level, s2s[0][1], nStart)
        else:
            been[nStart] = (level, nStart, s2s[0][0])
            #been[s2s[0][0]] = (level, s2s[0][0], nStart)
    else:
        been[nStart] = level

    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
        
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue

        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]
        
        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0

        v = None
        if v1 == nStart:
            v = v2
        else:
            v = v1

        if v not in been:
            goFo(v, been, s1Dict[v].s2s, s1Dict, psa=psa, level=level+1, byEdge=byEdge)

def goFoLevel(nStart, been, s2s, s1Dict, skipS2=None, psa=1, level=1, stopLevel=4, byEdge=0):
    
    # Recursive function that moves through a triangulation to find minimum network paths.
    # It starts at node (nStart), and is forced to move forward (cannot re-visit nStart).
    ######################################################################################
    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue_name == "PSA" and not psa:
        return 0

    if level == stopLevel:
        return 0

    s2s = sorted(s2s)
    
    # Record progress to force forward moves.
    #########################################
    if byEdge:
        if nStart == s2s[0][0]:
            been[nStart] = (level, nStart, s2s[0][1])
            #been[s2s[0][1]] = (level, s2s[0][1], nStart)
        else:
            been[nStart] = (level, nStart, s2s[0][0])
            #been[s2s[0][0]] = (level, s2s[0][0], nStart)
    else:
        been[nStart] = level

    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
        
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue

        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]
        
        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0

        v = None
        if v1 == nStart:
            v = v2
        else:
            v = v1

        if v not in been:
            goFoLevel(v, been, s1Dict[v].s2s, s1Dict, psa=psa, level=level+1, stopLevel=stopLevel, byEdge=byEdge)

def goFoCutOff(nStart, been, s2s, s1Dict, skipS2=None, psa=1, cutOff=1000.0):

    # Recursive function that moves through a triangulation to find minimum network paths.
    # It starts at node (nStart), and is forced to move forward (cannot re-visit nStart).
    ######################################################################################
    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue_name == "PSA" and not psa:
        return 0
    
    # Record progress to force forward moves.
    #########################################
    been.update({nStart:None})
    
    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
        
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue
    
        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]

        # Restrict goFo based on internode distance.
        ############################################
        networkEdgeLength = distance(s1Dict[v1].s1, s1Dict[v2].s1)

        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0
        
        # goFo from v1.
        ###############
        if v1 not in been and networkEdgeLength <= cutOff:
            goFoCutOff(v1, been, s1Dict[v1].s2s, s1Dict, psa=psa, cutOff=cutOff)

        # goFo from v2.
        ###############
        if v2 not in been and networkEdgeLength <= cutOff:
            goFoCutOff(v2, been, s1Dict[v2].s2s, s1Dict, psa=psa, cutOff=cutOff)

def goFoClassification(nStart, been, s2s, s1Dict, classification_nodes, skipS2=None, psa=1, level=1, byEdge=0):
    
    # Recursive function that moves through a triangulation to find minimum network paths.
    # It starts at node (nStart), and is forced to move forward (cannot re-visit nStart).
    ######################################################################################
    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue_name == "PSA" and not psa:
        return 0

    if nStart not in classification_nodes:
        return 0

    s2s = sorted(s2s)
    
    # Record progress to force forward moves.
    #########################################
    if byEdge:
        if nStart == s2s[0][0]:
            been[nStart] = (level, nStart, s2s[0][1])
            #been[s2s[0][1]] = (level, s2s[0][1], nStart)
        else:
            been[nStart] = (level, nStart, s2s[0][0])
            #been[s2s[0][0]] = (level, s2s[0][0], nStart)
    else:
        been[nStart] = level

    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
        
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue

        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]
        
        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0

        v = None
        if v1 == nStart:
            v = v2
        else:
            v = v1

        if v not in been:
            goFoClassification(v, been, s1Dict[v].s2s, s1Dict, classification_nodes, skipS2=skipS2, psa=psa, level=level+1, byEdge=byEdge)

def goFoParity(nStart, been, s2s, s1Dict, disulfideBonds, skipS2=None, psa=1, includeCys=1):
    
    # Same algorithm as goFo, but with additional conditional statements.
    #####################################################################

    # Get the node's residue key.
    #############################
    residueKey = s1Dict[nStart].s1.data.residue.key
    rearrangedKey = (residueKey[1], residueKey[0], residueKey[2])

    # Collect information on disulfide bonds: so that they can be skipped.
    ######################################################################
    disulfideBondResKeys = []
    for values in disulfideBonds.values():
        disulfideBondResKeys.extend(values)

    # Skip Cys residues involved in disulfide bonds.
    ################################################
    if rearrangedKey in disulfideBondResKeys:
        return 0
    
    # Relevant parity side chains.
    ##############################
    plusList = ("HIS", "LYS", "ARG")
    if includeCys:
        minusList = ("ASP", "GLU", "CYS")
    else:
        minusList = ("ASP", "GLU")

    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue.name == "PSA":
        return 0
    
    # Set node charge.
    ##################
    charge = "minus"
    if s1Dict[nStart].s1.data.residue.name in plusList:
        charge = "plus"
    
    # Record progress to force forward moves.
    #########################################
    been.update({nStart:None})

    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
    
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue

        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]
        
        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0
        
        for v in s2:

            if v != nStart:
                
                # Only goFo from neighboring networks nodes having opposite formal charge.
                ##########################################################################
                if charge == "plus":
                    # Node nStart is plus, so s2 node must be minus.
                    ################################################
                    if v not in been and s1Dict[v].s1.data.residue.name in minusList:
                        goFoParity(v, been, s1Dict[v].s2s, s1Dict, disulfideBonds, psa=psa)
                else:
                    # Node nStart is minus, so s2 node must be plus.
                    ################################################
                    if v not in been and s1Dict[v].s1.data.residue.name in plusList:
                        goFoParity(v, been, s1Dict[v].s2s, s1Dict, disulfideBonds, psa=psa)


def goFoAntiParity(nStart, been, s2s, s1Dict, charge, disulfideBonds, skipS2=None, psa=1, includeCys=1):
    

    # Same algorithm as goFo, but with additional conditional statements.
    # In this version of the parity function charge must be defined in the function call.
    #####################################################################################

    # Get the node's residue key.
    #############################
    residueKey = s1Dict[nStart].s1.data.residue.key

    # Collect information on disulfide bonds: so that they can be skipped.
    ######################################################################
    disulfideBondResKeys = list(disulfideBonds)

    # Skip Cys residues involved in disulfide bonds.
    ################################################
    if residueKey in disulfideBondResKeys:
        return 0
    
    # Relevant anti-parity side chains.
    ###################################
    plusList = ("HIS", "LYS", "ARG")
    if includeCys:
        minusList = ("ASP", "GLU", "CYS")
    else:
        minusList = ("ASP", "GLU")

    # Skip Cys bases on the includeCys argument.
    ############################################
    if s1Dict[nStart].s1.data.residue.name == "CYS" and not includeCys:
        return 0
        
    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue.name == "PSA" and not psa:
        return 0

    # Record progress to force forward moves.
    #########################################
    been.update({nStart:None})

    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
        
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue

        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]
        
        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0

        # Only goFo from neighboring networks nodes having the same formal charge.
        ##########################################################################
        if charge == "plus":
            if v1 not in been and s1Dict[v1].s1.data.residue.name in plusList:
                goFoAntiParity(v1, been, s1Dict[v1].s2s, s1Dict, charge, disulfideBonds, psa=psa)
            if v2 not in been and s1Dict[v2].s1.data.residue.name in plusList: #error, was s1Dict[v1].s1.data.residue.name
                goFoAntiParity(v2, been, s1Dict[v2].s2s, s1Dict, charge, disulfideBonds, psa=psa)
        elif charge == "minus":
            if v1 not in been and s1Dict[v1].s1.data.residue.name in minusList:
                goFoAntiParity(v1, been, s1Dict[v1].s2s, s1Dict, charge, disulfideBonds, psa=psa)
            if v2 not in been and s1Dict[v2].s1.data.residue.name in minusList: #error, was s1Dict[v1].s1.data.residue.name
                goFoAntiParity(v2, been, s1Dict[v2].s2s, s1Dict,charge, disulfideBonds, psa=psa)

def goFoChelators(nStart, been, s2s, s1Dict, psa=1, skipS2=None):
    
    # Same algorithm as goFo, but with additional conditional statements.
    # In this version of the parity function charge must be defined in the function call.
    #####################################################################################

    # Get the node's residue key.
    #############################
    residueKey = s1Dict[nStart].s1.data.residue.key
    rearrangedKey = (residueKey[1], residueKey[0], residueKey[2])
    
    # Relevant chelator side chains.
    #################################
    resList = ("HIS", "CYS")
        
    # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
    ##############################################################
    if s1Dict[nStart].s1.data.residue.name == "PSA":
        return 0

    # Record progress to force forward moves.
    #########################################
    been.update({nStart:None})

    # goFo along each edge incident on the network node "nStart".
    #############################################################
    for s2 in s2s:
        
        # Option to skip certain network edges (s2 Edge objects): seldom used.
        ######################################################################
        if skipS2:
            if s2 == skipS2:
                continue

        # goFo recursively from the network node "nStart".
        ##################################################
        # The two Vertex objects (either 3D or 4D) of the Edge object.
        ##############################################################
        v1, v2 = s2[0], s2[1]
        
        # Conditionally approve inclusion of PsuedoAtom "PSA" objects.
        ##############################################################
        if s1Dict[v1].s1.data.residue_name == "PSA" and not psa:
            return 0
        if s1Dict[v2].s1.data.residue_name == "PSA" and not psa:
            return 0

        # Only goFo from neighboring networks nodes having the same formal charge.
        ##########################################################################
        if v1 not in been and s1Dict[v1].s1.data.residue.name in resList:
            goFoChelators(v1, been, s1Dict[v1].s2s, s1Dict, psa=psa)
        if v2 not in been and s1Dict[v2].s1.data.residue.name in resList: 
            goFoChelators(v2, been, s1Dict[v2].s2s, s1Dict, psa=psa)

def pruneTriangulationGoFo(triangulation, maxNetworkEdgeLength, minNetworkEdgeLength=0, psa=1):
    
    # Dependencies.
    ###############
    from pHinder._vendor.compGeometry import findMinSidechainDistance
    
    # Prune the triangulation using the parameter maxNetworkEdgeLength.
    ###################################################################
    networkNodes = []
    for node in triangulation:
        
        # Save the node in the dictionary of complete triangulation nodes.
        ##################################################################
        nodeSaveResNum = triangulation[node].s1.data.residue.num
        nodeSaveResChn = triangulation[node].s1.data.residue.chn
        nodeSaveResName = triangulation[node].s1.data.residue.name
        
        # Identify network connections.
        ###############################
        s2Keep = []
        for s2 in triangulation[node].s2s:
            v1 = triangulation[s2[0]].s1
            v2 = triangulation[s2[1]].s1
            minRD = findMinSidechainDistance(v1.data, v2.data)
            if not minNetworkEdgeLength:
                if minRD < maxNetworkEdgeLength:
                    s2Keep.append(s2)
            else:
                if minRD < maxNetworkEdgeLength and minRD > minNetworkEdgeLength:
                    s2Keep.append(s2)

        # Update the s2s list of the node.
        ##################################
        if len(s2Keep) > 0:
            # Over-write the node s2s with only network s2s.
            ################################################
            triangulation[node].s2s = s2Keep
            networkNodes.append(node)
        else:
            # Erase s2s.
            ############
            triangulation[node].s2s = []
            # Add singular nodes.
            #####################
            networkNodes.append(node)

    # Use goFo to identify number of separate networks.
    ###################################################
    networks = {}
    for node in networkNodes:
        
        net = {}
        s2s = triangulation[node].s2s
        t = triangulation
        goFo(node, net, s2s, t, psa=psa)
        if len(net) > 0:
            net = sorted(net)
            net = tuple(net)
            networks.update({net:t})

    return (triangulation, networks)
