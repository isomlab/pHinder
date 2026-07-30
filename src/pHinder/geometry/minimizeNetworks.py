from copy import deepcopy
from pHinder._vendor.compGeometry import distance
from pHinder.geometry.goFo import goFo

def consolidated_networks(sub_networks, t):

    # Consolidate the sub-networks in the current network...
    #################################################################################################
    # 1) Match sub-networks to link...
    keys, consolidated_sub_networks, absorbed_sub_networks = list(sub_networks), {}, {}
    for k1 in keys:
        if k1 not in absorbed_sub_networks:
            for k2 in keys:
                if k1 != k2:
                    match = 0
                    for node in k1:
                        if node in k2:
                            match = 1
                            break
                    if match:
                        absorbed_sub_networks[k2] = None
                        if k1 not in consolidated_sub_networks:
                            consolidated_sub_networks[k1] = [k1, k2]
                        else:
                            consolidated_sub_networks[k1].append(k2)
    
    # 2) Consolidate sub-networks to link...
    consolidated_sub_networks_update = {}
    for consolidated_sub_network in consolidated_sub_networks:
        unique_nodes = {}
        for sub_network in consolidated_sub_networks[consolidated_sub_network]:
            for node in sub_network:
                unique_nodes[node] = t[node]
        consolidated_sub_networks_update[tuple(sorted(unique_nodes))] = unique_nodes
    consolidated_sub_networks = consolidated_sub_networks_update

    # 3) Internally link the sub-networks via s2s...
    for consolidated_sub_network in consolidated_sub_networks:
        for node_1 in consolidated_sub_network:
            for node_2 in consolidated_sub_network:
                if node_1 != node_2:
                    for s2 in t[node_2].s2s:
                        if node_1 in s2:
                            if s2 not in t[node_1].s2s and (s2[1], s2[0]) not in t[node_1].s2s:
                                t[node_1].s2s.append(s2)

    return consolidated_sub_networks

def join_networks(network_edges, sub_networks, mutable_triangulation, non_mutated_triangulation, max_edge_length=1e6):

    t = mutable_triangulation
    
    # 1) Identify the shortest edge shared by the sub-networks...
    i, keys, sub_network_matches = 0, sorted(sub_networks), {}
    while i < len(sub_networks) - 1:
        sub_network_1 = sub_networks[keys[i]]
        sub_network_1_neighbor_nodes = {}
        for node in sub_network_1:
            for s2 in t[node].s2s:
                for node_s2 in s2:
                    if node_s2 not in sub_network_1:
                        sub_network_1_neighbor_nodes[node_s2] = None
        j = i + 1
        while j < len(sub_networks):
            sub_network_2 = sub_networks[keys[j]]
            for node in sub_network_2:
                # The key here is to use the original, non mutated triangulation...
                for s2 in non_mutated_triangulation[node].s2s: 
                    for node_s2 in s2:
                        if node_s2 not in sub_network_2:
                            if node_s2 in sub_network_1:
                                d = distance(t[node].s1, t[node_s2].s1)
                                if d < max_edge_length:
                                    sub_network_matches[(d, (node, node_s2), tuple(sorted(sub_network_1)), tuple(sorted(sub_network_2)))] = None
            j += 1
        i += 1

    # 2) Externally link the consolidated sub-networks via their shortest shared sub-network edge...
    # 2a) Set the shortest linking edge...
    to_link, already_linked = {}, {}
    keys = sorted(sub_network_matches)
    for key in keys:
        sub_network_1, sub_network_2 = key[2], key[3]
        if sub_network_1 in already_linked or sub_network_2 in already_linked:
            continue
        already_linked[sub_network_1] = None
        already_linked[sub_network_2] = None
        to_link[key] = None
    # 2b) Create the network linkage...
    keys = sorted(to_link)
    for key in keys:
        node_1, node_2 = key[1][0], key[1][1]
        if (node_1, node_2) not in t[node_1].s2s and (node_2, node_1) not in t[node_1].s2s:
            t[node_1].s2s.append((node_1, node_2))
        if (node_1, node_2) not in t[node_2].s2s and (node_2, node_1) not in t[node_2].s2s:
            t[node_2].s2s.append((node_1, node_2))

def minimization_goFo(network, t):

    sub_networks, touched_nodes = {}, {}
    for node in network:

        # The minimum, most efficient goFo unit...
        if node in touched_nodes:
            continue

        net = {}
        s2s = t[node].s2s
        goFo(node, net, s2s, t, psa=1)
        if len(net) > 0:
            sub_networks[tuple(sorted(net))] = net
        touched_nodes.update(net) 

    return sub_networks

def minimizeNetworks(networks, hull4DCopy, max_edge_length=1e6):

    # Minimum network trace with no loops!
    print("\nStart -- optimizing networks")
    print(40*"-")

    # Find the triangle loops within the network path.
    ##################################################
    n_network = 1
    for network in networks:
        
        print("Optimizing network", n_network, "of", len(networks))
        n_network += 1

        nodes_with_edges, nodes_with_edge_distances = {}, {}
        for node in network:
            nodes_with_edges[(len(hull4DCopy[node].s2s), node)] = hull4DCopy[node].s2s[:]
            node_edge_distance_dict = {}
            for s2 in hull4DCopy[node].s2s[:]:
                v1 = hull4DCopy[s2[0]].s1
                v2 = hull4DCopy[s2[1]].s1
                d = distance(v1, v2)
                if d <= max_edge_length:
                    node_edge_distance_dict[d] = s2
            nodes_with_edge_distances[node] = node_edge_distance_dict

        # Mutable log of edges in the current network...
        network_edges = {}
        for node in network:
            for s2 in hull4DCopy[node].s2s:
                network_edges[s2] = None

        t = deepcopy(hull4DCopy)
        keys = sorted(nodes_with_edges)
        for k1 in keys:
            node = k1[1]
            node_distances = sorted(nodes_with_edge_distances[node])
            if node_distances:
                t[node].s2s = [nodes_with_edge_distances[node][node_distances[0]],]
            else:
                t[node].s2s = []
        for k1 in keys:
            for k2 in keys:
                if k1[1] != k2[1]:
                    if t[k1[1]].s2s == t[k2[1]].s2s:
                        t[k2[1]].s2s = []
        for k1 in keys:
            for s2 in t[k1[1]].s2s:
                del network_edges[s2]

        # Iteratively assemble the minimized sub-networks of the parent network...
        #################################################################################################
        # print("Optimizing sub-network")
        # print(40*"-")
        sub_networks, delta, i = minimization_goFo(network, t), 1, 0
        while delta:
            consolidated_sub_networks = consolidated_networks(sub_networks, t)
            new_sub_networks = minimization_goFo(network, t)
            join_networks(network_edges, new_sub_networks, t, hull4DCopy, max_edge_length=max_edge_length) 
            if new_sub_networks == sub_networks:
                delta = 0
            sub_networks = new_sub_networks
            i += 1
            # print("network minimization interation...", i)
        print()

        # Update the triangulation with the minimized network edges...
        #################################################################################################
        for node in network:
            hull4DCopy[node].s2s = t[node].s2s

    print("Done -- optimizing networks")

def getFinalMinimizedNetworks(networks, minNetworkSize=0):

    updated_networks = {}
    for network in networks:
        if minNetworkSize:
            if len(network) <= minNetworkSize:
                continue
        # print(40*"-")
        # print(network)
        updated_networks[network] = {}
        for node in network:
            updated_networks[network][node] = networks[network][node] 
            # print(node, networks[network][node])
    return updated_networks

