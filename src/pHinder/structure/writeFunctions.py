import gzip
import pHinder._vendor.pdbFile as pdbFile

# Have to write PDB formatted files instead of mmCIF formatted files,
# because mmCIF doesn't allow me to output CONECT data for the network edges.

def writeNetworkToFile(networks, triangulation, networkName, outputDirectory, pdbCode,
                       chain=None, include_residues=1, zipit=0):

    # Rank order networks by size...
    ###################################################################################
    network_keys = list(networks)
    ranked_network_keys = []

    for network_key in network_keys:
        ranked_network_keys.append((len(network_key), network_key))

    ranked_network_keys = sorted(ranked_network_keys)
    ranked_network_keys.reverse()

    # Write PDB and CONECT information order by decreasing network sizes...
    ###################################################################################
    residue_networks, network_connect = {}, {}

    for ranked_network_key in ranked_network_keys:
        network = networks[ranked_network_key[1]]
        residue_network_nodes, connect = {}, {}

        for node in sorted(network):
            # Need to include atom_serial in case residue key is the same from different PDB contributions...
            node_residue_key = tuple(
                [triangulation[node].s1.data.atom_serial] +
                list(triangulation[node].s1.data.residue.key)
            )
            node_residue = triangulation[node].s1.data.residue
            residue_network_nodes[node_residue_key] = node_residue
            connect[node_residue_key] = {node_residue: None}

            for s2 in triangulation[node].s2s:
                for neighbor_node in s2:
                    # Need to include atom_serial in case residue key is the same from different PDB contributions...
                    neighbor_residue_key = tuple(
                        [triangulation[neighbor_node].s1.data.atom_serial] +
                        list(triangulation[neighbor_node].s1.data.residue.key)
                    )
                    neighbor_residue = triangulation[neighbor_node].s1.data.residue
                    residue_network_nodes[neighbor_residue_key] = neighbor_residue
                    connect[node_residue_key][neighbor_residue] = None

        network_connect[ranked_network_key[1]] = connect
        residue_networks[ranked_network_key[1]] = residue_network_nodes

    # Write full amino acid side chain PDB version with network and CONECT information...
    ####################################################################################
    if include_residues:
        # Open output for residue-level network
        if zipit:
            if chain:
                outFile = (outputDirectory + pdbCode + "-" + networkName +
                           "ResidueNetwork-" + chain + ".pdb.gz")
                outputTopoPDB = gzip.open(outFile, "wt")
            else:
                outFile = (outputDirectory + pdbCode + "-" + networkName +
                           "ResidueNetwork.pdb.gz")
                outputTopoPDB = gzip.open(outFile, "wt")
        else:
            if chain:
                outFile = (outputDirectory + pdbCode + "-" + networkName +
                           "ResidueNetwork-" + chain + ".pdb")
                outputTopoPDB = open(outFile, "w")
            else:
                outFile = (outputDirectory + pdbCode + "-" + networkName +
                           "ResidueNetwork.pdb")
                outputTopoPDB = open(outFile, "w")

        network_number, atoms_string = 1, ""
        for residue_network in residue_networks:
            keys = sorted(residue_networks[residue_network])
            for key in keys:
                residue = residue_networks[residue_network][key]
                for atom_key in residue.atoms:
                    atom = residue.atoms[atom_key]

                    # Use the atom as-is; it already came from pdbFile.PDBfile
                    atom.segment_identifier = network_number
                    atoms_string += str(atom)

                network_number += 1

        outputTopoPDB.write(atoms_string)

    # Write CONECT information
    conectString, network_number = "", 1

    for network in network_connect:
        connect = network_connect[network]
        conectString += ("\nREMARK BEGIN NETWORK " + str(network_number) +
                         " : " + str(len(connect)) + " nodes\n")

        for network_node in sorted(connect):
            conectString += "%-6s" % "CONECT"
            network_nodes = connect[network_node]

            for node in network_nodes:
                tsc_atom = node.get_terminal_sidechain_atom()

                if tsc_atom:
                    # The -1 corrects for how the tscAtom is instantiated...
                    atom_serial = tsc_atom.atom_serial - 1
                else:
                    atom_key = list(node.atoms)[0]
                    tsc_atom = node.atoms[atom_key]
                    # tscAtom is PSA...
                    atom_serial = tsc_atom.atom_serial

                conectString += "%5i" % atom_serial

            conectString += "\n"

        conectString += ("REMARK END NETWORK " + str(network_number) +
                         " : " + str(len(connect)) + " nodes\n")
        network_number += 1

    if include_residues:
        outputTopoPDB.write(conectString)
        outputTopoPDB.close()

    # Write terminal side chain atom PDB version with network and CONECT information...
    ###################################################################################
    if zipit:
        if chain:
            outFile = (outputDirectory + pdbCode + "-" + networkName +
                       "AtomNetwork-" + chain + ".pdb.gz")
            outputTopoPDB = gzip.open(outFile, "wt")
        else:
            outFile = (outputDirectory + pdbCode + "-" + networkName +
                       "AtomNetwork.pdb.gz")
            outputTopoPDB = gzip.open(outFile, "wt")
    else:
        if chain:
            outFile = (outputDirectory + pdbCode + "-" + networkName +
                       "AtomNetwork-" + chain + ".pdb")
            outputTopoPDB = open(outFile, "w")
        else:
            outFile = (outputDirectory + pdbCode + "-" + networkName +
                       "AtomNetwork.pdb")
            outputTopoPDB = open(outFile, "w")

    network_number, atoms_string = 1, ""

    for residue_network in residue_networks:
        keys = sorted(residue_networks[residue_network])
        for key in keys:
            residue = residue_networks[residue_network][key]
            tsc_atom = residue.get_terminal_sidechain_atom()

            out_atom = tsc_atom
            out_atom.segment_identifier = network_number
            # Mirror the -1 correction applied on the CONECT side (~L109):
            # get_terminal_sidechain_atom() instantiates the PseudoAtom with
            # atom_serial = parent_atom.atom_serial + 1, so we strip the +1
            # here to keep the written atom_serial aligned with the serials
            # referenced by CONECT records.
            out_atom.atom_serial -= 1
            atoms_string += str(out_atom)

            network_number += 1

    outputTopoPDB.write(atoms_string)
    outputTopoPDB.write(conectString)
    outputTopoPDB.close()

