# Import dependencies here for enhanced performance.

from os import sep
import gzip
from pHinder._vendor.compGeometry import Vertex
from decimal import *

singleLetter = {
    "ALA": "A", "ASP": "D", "ASN": "N", "ARG": "R", "CYS": "C",
    "GLY": "G", "GLU": "E", "GLN": "Q", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PRO": "P", "PHE": "F",
    "SER": "S", "THR": "T", "TYR": "Y", "TRP": "W", "VAL": "V"
}

NONPOLAR = ["ALA", "ILE", "LEU", "MET", "PHE", "PRO", "VAL", "GLY"]
POLAR = ["ASN", "GLN", "SER", "THR", "TYR", "TRP"]
IONIZABLE = ["ASP", "GLU", "HIS", "CYS", "LYS", "ARG", "MLY", "MLZ", "M3L"]
ALL_SIDECHAINS = NONPOLAR + POLAR + IONIZABLE
ACTIVE = ["GLY", "ARG", "ASP", "GLU", "HIS", "LYS", "CYS", "SER", "TYR", "THR", "PHE", "TRP"]


def assess_cif_fields(cif_options):
    cif_options_all = {
        "_atom_site.group_PDB": False,
        "_atom_site.id": False,
        "_atom_site.type_symbol": False,
        "_atom_site.label_atom_id": False,
        "_atom_site.label_alt_id": False,
        "_atom_site.label_comp_id": False,
        "_atom_site.label_asym_id": False,
        "_atom_site.label_entity_id": False,
        "_atom_site.label_seq_id": False,
        "_atom_site.pdbx_PDB_ins_code": False,
        "_atom_site.Cartn_x": False,
        "_atom_site.Cartn_y": False,
        "_atom_site.Cartn_z": False,
        "_atom_site.occupancy": False,
        "_atom_site.B_iso_or_equiv": False,
        "_atom_site.Cartn_x_esd": False,
        "_atom_site.Cartn_y_esd": False,
        "_atom_site.Cartn_z_esd": False,
        "_atom_site.occupancy_esd": False,
        "_atom_site.B_iso_or_equiv_esd": False,
        "_atom_site.pdbx_formal_charge": False,
        "_atom_site.auth_seq_id": False,
        "_atom_site.auth_comp_id": False,
        "_atom_site.auth_asym_id": False,
        "_atom_site.auth_atom_id": False,
        "_atom_site.pdbx_PDB_model_num": False
    }
    for cif_option in cif_options.split("\n"):
        if cif_option in cif_options_all:
            cif_options_all[cif_option] = True
    return cif_options_all


class Atom_LO:
    def __init__(self):
        self.atom_serial = None
        self.residueKey = None
        self.repr = ""

    def __repr__(self):
        return self.repr


class PseudoAtom:
    def __init__(self):
        self.format = "pdb"
        self.pdbfileline = ""
        self.record_name = "ATOM  "
        self.atom_serial = 1
        self.atom_name = " O "
        self.atom_alternate_location = ""
        self.residue_name = "PSA"
        self.chain_identifier = "A"
        self.residue_sequence_number = 1
        self.residue_insertion_code = ""
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.occupancy = 0.0
        self.temperature_factor = 0.0
        self.segment_identifier = ""
        self.symbol = ""
        self.charge = 0
        self.atom_key = (self.residue_sequence_number, self.residue_name,
                         self.chain_identifier, self.atom_name)
        self.residue_key = (self.residue_sequence_number, self.residue_name,
                            self.chain_identifier)
        self.residue = None
        self.exposed, self.extended = 0, 0
        self.margin = 0
        self.core = 0
        self.pdbCode = ""
        self.Atom_LO = Atom_LO()
        self.Atom_LO.atom_serial = self.atom_serial
        self.Atom_LO.residueKey = self.residue_key
        self.Atom_LO.repr = self.__repr__()
        self.v = Vertex((self.x, self.y, self.z), data=self,
                        unique_id=self.atom_serial)

    def reinitialize(self):
        if self.residue:
            self.residue.num = self.residue_sequence_number
            self.residue.name = self.residue_name
            self.residue.chn = self.chain_identifier
        self.residue_key = (self.residue_sequence_number, self.residue_name,
                            self.chain_identifier)
        self.atom_key = (self.residue_sequence_number, self.residue_name,
                         self.chain_identifier, self.atom_name)
        self.Atom_LO = Atom_LO()
        self.Atom_LO.atom_serial = self.atom_serial
        self.Atom_LO.residueKey = self.residue_key
        self.Atom_LO.repr = self.__repr__()
        self.v = Vertex((self.x, self.y, self.z), data=self,
                        unique_id=self.atom_serial)

    def copy(self):
        return Atom(str(self))

    def __repr__(self, forceHex=0):
        chain_identifier = " " if self.chain_identifier == "NULL" else self.chain_identifier
        stringAtomSerial = "%5i" % self.atom_serial
        stringResidueNumber = "%3i" % self.residue_sequence_number
        pdbformat = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
        if len(chain_identifier) > 1:
            pdbformat = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
        pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2i\n"
        return pdbformat % (
            self.record_name,
            stringAtomSerial,
            self.atom_name,
            self.atom_alternate_location,
            self.residue_name,
            chain_identifier,
            stringResidueNumber,
            self.residue_insertion_code,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.temperature_factor,
            self.segment_identifier,
            self.symbol,
            self.charge
        )

class Atom:
    def __init__(self, pdbfileline="", cif_line="", cif_options=None,
                 legalize=0, twoCharacterChain=0):
        self.format = "pdb"  # "pdb" or "cif"
        self.pdbfileline = pdbfileline
        self.cif_line = cif_line
        self.cif_options = cif_options or {}
        self.pdb_code = " "

        self.record_name = " "
        self.atom_serial = 0
        self.atom_name = " "
        self.atom_alternate_location = " "
        self.residue_name = ""
        self.chain_identifier = " "
        self.entity_id = " "
        self.residue_sequence_number = 0
        self.residue_insertion_code = " "
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.occupancy = 0.0
        self.temperature_factor = 0.0
        self.sigma_x = 0.0
        self.sigma_y = 0.0
        self.sigma_z = 0.0
        self.sigma_occupancy = 0.0
        self.sigma_temperature_factor = 0.0
        self.segment_identifier = " "
        self.symbol = " "
        self.charge = 0
        self.atom_key = ""
        self.residue_key = ""
        self.residue = None
        self.exposed, self.extended = 0, 0
        self.margin = 0
        self.core = 0
        self.hex_atom_serial = 0
        self.hex_residue_sequence_number = 0
        self.cif_dictionary = {}

        if pdbfileline:
            self._init_from_pdb_line(pdbfileline)
            self.format = "pdb"
        elif cif_line and cif_options:
            self._init_from_cif_line(cif_line, cif_options)
            self.format = "cif"

        self.Atom_LO = Atom_LO()
        self.Atom_LO.atom_serial = self.atom_serial
        self.Atom_LO.residueKey = self.residue_key
        self.Atom_LO.repr = self.__repr__()
        self.v = Vertex((self.x, self.y, self.z), data=self,
                        unique_id=self.atom_serial)

    def _init_from_pdb_line(self, line):
        self.record_name = line[0:6]
        try:
            self.atom_serial = int(line[6:11])
        except Exception:
            self.atom_serial = int(line[6:11], 16)
            self.hex_atom_serial = 1

        self.atom_name = line[12:16]
        self.atom_alternate_location = line[16:17]
        self.residue_name = line[17:20]
        self.chain_identifier = line[21:22]
        try:
            self.residue_sequence_number = int(line[22:26])
        except Exception:
            self.residue_sequence_number = int(line[22:26], 16)
            self.hex_residue_sequence_number = 1

        self.residue_insertion_code = line[26:27]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])

        if line[54:60]:
            self.occupancy = float(line[54:60])
        if line[60:66]:
            self.temperature_factor = float(line[60:66])

        self.segment_identifier = line[72:76]
        self.symbol = line[76:78]
        try:
            self.charge = int(line[78:79])
        except Exception:
            self.charge = 0

        self.residue_key = (self.residue_sequence_number,
                            self.residue_name,
                            self.chain_identifier)
        self.atom_key = (self.residue_sequence_number,
                         self.residue_name,
                         self.chain_identifier,
                         self.atom_name)

    def _init_from_cif_line(self, cif_line, cif_options):
        atom_details = cif_line.strip().split()
        cif_option_keys = list(cif_options)
        i = 0

        while i < len(atom_details) and i < len(cif_option_keys):
            atom_detail = atom_details[i]
            cif_option = cif_option_keys[i]

            if cif_option == "_atom_site.group_PDB":
                # e.g. "ATOM" / "HETATM" → PDB-style 6-char field
                try:
                    self.record_name = "%-6s" % atom_detail
                    self.cif_dictionary["_atom_site.group_PDB"] = self.record_name
                except Exception:
                    pass

            if cif_option == "_atom_site.id":
                try:
                    self.atom_serial = int(atom_detail)
                    self.cif_dictionary["_atom_site.id"] = self.atom_serial
                except Exception:
                    pass

            if cif_option == "_atom_site.type_symbol":
                self.symbol = atom_detail
                self.cif_dictionary["_atom_site.type_symbol"] = self.symbol

            if cif_option == "_atom_site.label_atom_id":
                # Atom name; treat "." / "?" as empty
                self.atom_name = atom_detail
                if self.atom_name in [".", "?"]:
                    self.atom_name = ""
                self.cif_dictionary["_atom_site.label_atom_id"] = self.atom_name

            if cif_option == "_atom_site.label_alt_id":
                # AltLoc; "." / "?" → blank
                self.atom_alternate_location = atom_detail
                if self.atom_alternate_location in [".", "?"]:
                    self.atom_alternate_location = " "
                self.cif_dictionary["_atom_site.label_alt_id"] = self.atom_alternate_location

            if cif_option == "_atom_site.label_comp_id":
                # Residue name; normalize "." / "?"
                self.residue_name = atom_detail
                if self.residue_name in [".", "?"]:
                    self.residue_name = ""
                self.cif_dictionary["_atom_site.label_comp_id"] = self.residue_name

            if cif_option == "_atom_site.label_asym_id":
                self.chain_identifier = atom_detail
                self.cif_dictionary["_atom_site.label_asym_id"] = self.chain_identifier

            # auth_asym_id is what a viewer shows and what a person means by
            # "chain R". It sits after label_asym_id in the atom_site loop, so
            # assigning it here overwrites the label value for this row, and
            # falls back to the label when a file carries no author ids.
            if cif_option == "_atom_site.auth_asym_id":
                self.chain_identifier = atom_detail
                self.cif_dictionary["_atom_site.auth_asym_id"] = self.chain_identifier

            if cif_option == "_atom_site.label_entity_id":
                self.entity_id = atom_detail
                self.cif_dictionary["_atom_site.label_entity_id"] = self.entity_id

            if cif_option == "_atom_site.label_seq_id":
                try:
                    self.residue_sequence_number = int(atom_detail)
                    self.cif_dictionary["_atom_site.label_seq_id"] = self.residue_sequence_number
                except Exception:
                    pass

            # Author residue numbering, for the same reason and with the same
            # ordering. Chain and residue number have to come from one system
            # or a residue ends up identified by a pair that matches nothing.
            if cif_option == "_atom_site.auth_seq_id":
                try:
                    self.residue_sequence_number = int(atom_detail)
                    self.cif_dictionary["_atom_site.auth_seq_id"] = self.residue_sequence_number
                except Exception:
                    pass

            if cif_option == "_atom_site.pdbx_PDB_ins_code":
                # Insertion code; "." / "?" → blank
                self.residue_insertion_code = atom_detail
                if self.residue_insertion_code in [".", "?"]:
                    self.residue_insertion_code = " "
                self.cif_dictionary["_atom_site.pdbx_PDB_ins_code"] = self.residue_insertion_code

            if cif_option == "_atom_site.Cartn_x":
                try:
                    token = str(atom_detail).strip().split()[0]
                    if token != ".":
                        self.x = float(token)
                        self.cif_dictionary["_atom_site.Cartn_x"] = self.x
                except Exception:
                    pass

            if cif_option == "_atom_site.Cartn_y":
                try:
                    token = str(atom_detail).strip().split()[0]
                    if token != ".":
                        self.y = float(token)
                        self.cif_dictionary["_atom_site.Cartn_y"] = self.y
                except Exception:
                    pass

            if cif_option == "_atom_site.Cartn_z":
                try:
                    token = str(atom_detail).strip().split()[0]
                    if token != ".":
                        self.z = float(token)
                        self.cif_dictionary["_atom_site.Cartn_z"] = self.z
                except Exception:
                    pass

            if cif_option == "_atom_site.occupancy":
                try:
                    self.occupancy = float(atom_detail)
                    self.cif_dictionary["_atom_site.occupancy"] = self.occupancy
                except Exception:
                    pass

            if cif_option == "_atom_site.B_iso_or_equiv":
                try:
                    self.temperature_factor = float(atom_detail)
                    self.cif_dictionary["_atom_site.B_iso_or_equiv"] = self.temperature_factor
                except Exception:
                    pass

            # (sigma fields and other keys can be copied straight from CIFAtom as needed)

            i += 1

        self.residue_key = (
            self.residue_sequence_number,
            self.residue_name,
            self.chain_identifier,
        )
        self.atom_key = (
            self.residue_sequence_number,
            self.residue_name,
            self.chain_identifier,
            self.atom_name,
        )


    def reinitialize(self):
        if self.residue:
            self.residue.num = self.residue_sequence_number
            self.residue.name = self.residue_name
            self.residue.chn = self.chain_identifier
        self.residue_key = (self.residue_sequence_number, self.residue_name,
                            self.chain_identifier)
        self.atom_key = (self.residue_sequence_number, self.residue_name,
                         self.chain_identifier, self.atom_name)
        self.Atom_LO = Atom_LO()
        self.Atom_LO.atom_serial = self.atom_serial
        self.Atom_LO.residueKey = self.residue_key
        self.Atom_LO.repr = self.__repr__()
        self.v = Vertex((self.x, self.y, self.z), data=self,
                        unique_id=self.atom_serial)

    def copy(self):
        return Atom(str(self))

    def __repr__(self, forceHex=0):
        if self.chain_identifier == "NULL":
            chain_identifier = " "
        else:
            chain_identifier = self.chain_identifier
        if self.hex_atom_serial and forceHex:
            stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
        else:
            stringAtomSerial = "%5i" % self.atom_serial
        if self.hex_residue_sequence_number and forceHex:
            stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
        else:
            stringResidueNumber = "%3i" % self.residue_sequence_number
        pdbformat = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
        if len(chain_identifier) > 1:
            pdbformat = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
        pdbformat += "   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s%2s\n"
        return pdbformat % (
            self.record_name,
            stringAtomSerial,
            self.atom_name,
            self.atom_alternate_location,
            self.residue_name,
            chain_identifier,
            stringResidueNumber,
            self.residue_insertion_code,
            self.x,
            self.y,
            self.z,
            self.occupancy,
            self.temperature_factor,
            self.segment_identifier,
            self.symbol,
            self.charge
        )

# class Atom:
#     def __init__(self, pdbfileline="", legalize=0, twoCharacterChain=0):
#         self.format = "pdb"
#         self.pdbfileline = pdbfileline
#         self.record_name = ""
#         self.atom_serial = 0
#         self.atom_name = ""
#         self.atom_alternate_location = ""
#         self.residue_name = ""
#         self.chain_identifier = ""
#         self.residue_sequence_number = 0
#         self.residue_insertion_code = ""
#         self.x = 0.0
#         self.y = 0.0
#         self.z = 0.0
#         self.occupancy = 0.0
#         self.temperature_factor = 0.0
#         self.segment_identifier = ""
#         self.symbol = ""
#         self.charge = 0
#         self.atom_key = ""
#         self.residue_key = ""
#         self.residue = None
#         self.exposed, self.extended = 0, 0
#         self.margin = 0
#         self.core = 0
#         self.hex_atom_serial = 0
#         self.hex_residue_sequence_number = 0
#         if pdbfileline:
#             line = pdbfileline
#             self.record_name = line[0:6]
#             try:
#                 self.atom_serial = int(line[6:11])
#             except:
#                 self.atom_serial = int(line[6:11], 16)
#                 self.hex_atom_serial = 1
#             self.atom_name = line[12:16]
#             self.atom_alternate_location = line[16:17]
#             self.residue_name = line[17:20]
#             self.chain_identifier = line[21:22]
#             try:
#                 self.residue_sequence_number = int(line[22:26])
#             except:
#                 self.residue_sequence_number = int(line[22:26], 16)
#                 self.hex_residue_sequence_number = 1
#             self.residue_insertion_code = line[26:27]
#             self.x = float(line[30:38])
#             self.y = float(line[38:46])
#             self.z = float(line[46:54])
#             if line[54:60]:
#                 self.occupancy = float(line[54:60])
#             if line[60:66]:
#                 self.temperature_factor = float(line[60:66])
#             self.segment_identifier = line[72:76]
#             self.symbol = line[76:78]
#             try:
#                 self.charge = int(line[78:79])
#             except:
#                 self.charge = 0
#         self.residue_key = (self.residue_sequence_number, self.residue_name,
#                             self.chain_identifier)
#         self.atom_key = (self.residue_sequence_number, self.residue_name,
#                          self.chain_identifier, self.atom_name)
#         self.Atom_LO = Atom_LO()
#         self.Atom_LO.atom_serial = self.atom_serial
#         self.Atom_LO.residueKey = self.residue_key
#         self.Atom_LO.repr = self.__repr__()
#         self.v = Vertex((self.x, self.y, self.z), data=self,
#                         unique_id=self.atom_serial)

#     def reinitialize(self):
#         if self.residue:
#             self.residue.num = self.residue_sequence_number
#             self.residue.name = self.residue_name
#             self.residue.chn = self.chain_identifier
#         self.residue_key = (self.residue_sequence_number, self.residue_name,
#                             self.chain_identifier)
#         self.atom_key = (self.residue_sequence_number, self.residue_name,
#                          self.chain_identifier, self.atom_name)
#         self.Atom_LO = Atom_LO()
#         self.Atom_LO.atom_serial = self.atom_serial
#         self.Atom_LO.residueKey = self.residue_key
#         self.Atom_LO.repr = self.__repr__()
#         self.v = Vertex((self.x, self.y, self.z), data=self,
#                         unique_id=self.atom_serial)

#     def copy(self):
#         return Atom(str(self))

#     def __repr__(self, forceHex=0):
#         if self.chain_identifier == "NULL":
#             chain_identifier = " "
#         else:
#             chain_identifier = self.chain_identifier
#         if self.hex_atom_serial and forceHex:
#             stringAtomSerial = "%5s" % hex(self.atom_serial)[2:]
#         else:
#             stringAtomSerial = "%5i" % self.atom_serial
#         if self.hex_residue_sequence_number and forceHex:
#             stringResidueNumber = "%3s" % hex(self.residue_sequence_number)[2:]
#         else:
#             stringResidueNumber = "%3i" % self.residue_sequence_number
#         pdbformat = "%6s%5s %-4s%1s%3s %-1s%4s%1s"
#         if len(chain_identifier) > 1:
#             pdbformat = "%6s%5s %-4s%1s%3s %-2s%3s%1s"
#         pdbformat += " %8.3f%8.3f%8.3f%6.2f%6.2f %-4s%2s%2s\n"
#         return pdbformat % (
#             self.record_name,
#             stringAtomSerial,
#             self.atom_name,
#             self.atom_alternate_location,
#             self.residue_name,
#             chain_identifier,
#             stringResidueNumber,
#             self.residue_insertion_code,
#             self.x,
#             self.y,
#             self.z,
#             self.occupancy,
#             self.temperature_factor,
#             self.segment_identifier,
#             self.symbol,
#             self.charge
#         )


# class CIFAtom:
#     def __init__(self, pdb_file_line="", cif_options=None):
#         cif_options = cif_options or {}
#         self.format = "cif"
#         self.pdb_file_line = pdb_file_line
#         self.cif_options = [k for k, v in cif_options.items() if v]
#         self.pdb_code = " "
#         self.record_name = " "
#         self.atom_serial = 0
#         self.atom_name = " "
#         self.atom_alternate_location = " "
#         self.residue_name = ""
#         self.chain_identifier = " "
#         self.entity_id = " "
#         self.residue_sequence_number = 0
#         self.residue_insertion_code = " "
#         self.x = 0.0
#         self.y = 0.0
#         self.z = 0.0
#         self.occupancy = 0.0
#         self.temperature_factor = 0.0
#         self.sigma_x = 0.0
#         self.sigma_y = 0.0
#         self.sigma_z = 0.0
#         self.sigma_occupancy = 0.0
#         self.sigma_temperature_factor = 0.0
#         self.segment_identifier = " "
#         self.symbol = " "
#         self.charge = 0
#         self.atom_key = " "
#         self.residue_key = " "
#         self.residue = None
#         self.exposed, self.extended = 0, 0
#         self.margin = 0
#         self.core = 0
#         self.cif_dictionary = {}
#         if self.pdb_file_line and cif_options:
#             atom_details = self.pdb_file_line.strip().split()
#             cif_option_keys = list(cif_options)
#             i = 0
#             while i < len(atom_details) and i < len(cif_option_keys):
#                 atom_detail = atom_details[i]
#                 cif_option = cif_option_keys[i]
#                 if cif_option == "_atom_site.group_PDB":
#                     try:
#                         self.record_name = "%-6s" % atom_detail
#                         self.cif_dictionary["_atom_site.group_PDB"] = self.record_name
#                     except:
#                         pass
#                 if cif_option == "_atom_site.id":
#                     try:
#                         self.atom_serial = int(atom_detail)
#                         self.cif_dictionary["_atom_site.id"] = self.atom_serial
#                     except:
#                         pass
#                 if cif_option == "_atom_site.type_symbol":
#                     try:
#                         self.symbol = atom_detail
#                         self.cif_dictionary["_atom_site.type_symbol"] = self.symbol
#                     except:
#                         pass
#                 if cif_option == "_atom_site.label_atom_id":
#                     try:
#                         self.atom_name = atom_detail
#                         self.cif_dictionary["_atom_site.label_atom_id"] = self.atom_name
#                     except:
#                         pass
#                 if cif_option == "_atom_site.label_alt_id":
#                     try:
#                         self.atom_alternate_location = atom_detail
#                         self.cif_dictionary["_atom_site.label_alt_id"] = self.atom_alternate_location
#                     except:
#                         pass
#                 if cif_option == "_atom_site.label_comp_id":
#                     try:
#                         self.residue_name = atom_detail
#                         self.cif_dictionary["_atom_site.label_comp_id"] = self.residue_name
#                     except:
#                         pass
#                 if cif_option == "_atom_site.label_asym_id":
#                     try:
#                         self.chain_identifier = atom_detail
#                         self.cif_dictionary["_atom_site.label_asym_id"] = self.chain_identifier
#                     except:
#                         pass
#                 if cif_option == "_atom_site.label_entity_id":
#                     try:
#                         self.entity_id = atom_detail
#                         self.cif_dictionary["_atom_site.label_entity_id"] = self.entity_id
#                     except:
#                         pass
#                 if cif_option == "_atom_site.label_seq_id":
#                     try:
#                         self.residue_sequence_number = int(atom_detail)
#                         self.cif_dictionary["_atom_site.label_seq_id"] = self.residue_sequence_number
#                     except:
#                         pass
#                 if cif_option == "_atom_site.pdbx_PDB_ins_code":
#                     try:
#                         self.residue_insertion_code = atom_detail
#                         self.cif_dictionary["_atom_site.pdbx_PDB_ins_code"] = self.residue_insertion_code
#                     except:
#                         pass
#                 if cif_option == "_atom_site.Cartn_x":
#                     try:
#                         # Clean atom_detail in case it has extra tokens (e.g. "93.740 1")
#                         token = str(atom_detail).strip().split()[0]
#                         if token != ".":
#                             self.x = float(token)
#                             print(self.x)
#                             self.cif_dictionary["_atom_site.Cartn_x"] = self.x
#                     except Exception:
#                         pass

#                 if cif_option == "_atom_site.Cartn_y":
#                     try:
#                         token = str(atom_detail).strip().split()[0]
#                         if token != ".":
#                             self.y = float(token)
#                             print(self.y)
#                             self.cif_dictionary["_atom_site.Cartn_y"] = self.y
#                     except Exception:
#                         pass

#                 if cif_option == "_atom_site.Cartn_z":
#                     try:
#                         token = str(atom_detail).strip().split()[0]
#                         if token != ".":
#                             self.z = float(token)
#                             print(self.z)
#                             self.cif_dictionary["_atom_site.Cartn_z"] = self.z
#                     except Exception:
#                         pass

#                 if cif_option == "_atom_site.occupancy":
#                     try:
#                         self.occupancy = float(atom_detail)
#                         self.cif_dictionary["_atom_site.occupancy"] = self.occupancy
#                     except:
#                         pass
#                 if cif_option == "_atom_site.B_iso_or_equiv":
#                     try:
#                         self.temperature_factor = float(atom_detail)
#                         self.cif_dictionary["_atom_site.B_iso_or_equiv"] = self.temperature_factor
#                     except:
#                         pass
#                 if cif_option == "_atom_site.Cartn_x_esd":
#                     try:
#                         self.sigma_x = float(atom_detail)
#                         self.cif_dictionary["_atom_site.Cartn_x_esd"] = self.sigma_x
#                     except:
#                         pass
#                 if cif_option == "_atom_site.Cartn_y_esd":
#                     try:
#                         self.sigma_y = float(atom_detail)
#                         self.cif_dictionary["_atom_site.Cartn_y_esd"] = self.sigma_y
#                     except:
#                         pass
#                 if cif_option == "_atom_site.Cartn_z_esd":
#                     try:
#                         self.sigma_z = float(atom_detail)
#                         self.cif_dictionary["_atom_site.Cartn_z_esd"] = self.sigma_z
#                     except:
#                         pass
#                 if cif_option == "_atom_site.occupancy_esd":
#                     try:
#                         self.sigma_occupancy = float(atom_detail)
#                         self.cif_dictionary["_atom_site.occupancy_esd"] = self.sigma_occupancy
#                     except:
#                         pass
#                 if cif_option == "_atom_site.B_iso_or_equiv_esd":
#                     try:
#                         self.sigma_temperature_factor = float(atom_detail)
#                         self.cif_dictionary["_atom_site.B_iso_or_equiv_esd"] = self.sigma_temperature_factor
#                     except:
#                         pass
#                 if cif_option == "_atom_site.pdbx_formal_charge":
#                     try:
#                         self.charge = float(atom_detail)
#                         self.cif_dictionary["_atom_site.pdbx_formal_charge"] = self.charge
#                     except:
#                         pass
#                 i += 1
#         self.residue_key = (self.residue_sequence_number, self.residue_name,
#                             self.chain_identifier)
#         self.atom_key = (self.residue_sequence_number, self.residue_name,
#                          self.chain_identifier, self.atom_name)
#         self.v = Vertex((self.x, self.y, self.z), data=self,
#                         unique_id=self.atom_serial)

#     def reinitialize(self):
#         if self.residue:
#             self.residue.num = self.residue_sequence_number
#             self.residue.name = self.residue_name
#             self.residue.chn = self.chain_identifier
#         self.residue_key = (self.residue_sequence_number, self.residue_name,
#                             self.chain_identifier)
#         self.atom_key = (self.residue_sequence_number, self.residue_name,
#                          self.chain_identifier, self.atom_name)
#         self.v = Vertex((self.x, self.y, self.z), data=self,
#                         unique_id=self.atom_serial)

#     def get_pdb_format(self):
#         atom = Atom()
#         atom.record_name = self.record_name
#         atom.atom_serial = self.atom_serial
#         atom.atom_name = self.atom_name
#         if self.atom_alternate_location == ".":
#             self.atom_alternate_location = ""
#         atom.atom_alternate_location = self.atom_alternate_location
#         atom.residue_name = self.residue_name
#         atom.chain_identifier = self.chain_identifier
#         atom.residue_sequence_number = self.residue_sequence_number
#         if self.residue_insertion_code == "?":
#             self.residue_insertion_code = ""
#         atom.residue_insertion_code = self.residue_insertion_code
#         atom.x = self.x
#         atom.y = self.y
#         atom.z = self.z
#         atom.occupancy = self.occupancy
#         atom.temperature_factor = self.temperature_factor
#         atom.segment_identifier = self.segment_identifier
#         atom.symbol = self.symbol
#         atom.charge = self.charge
#         return str(atom)

#     def __repr__(self):
#         atom_format = ""
#         for cif_option in self.cif_options:
#             if cif_option in self.cif_dictionary:
#                 atom_format += str(self.cif_dictionary[cif_option]) + " "
#         return atom_format + "\n"


class Residue:
    def __init__(self, Atom_instance):
        self.nonSidechainAtomNames = ['N', 'CA', 'C', 'O', 'D', 'H', 'DA']
        self.num = Atom_instance.residue_sequence_number
        self.name = Atom_instance.residue_name
        self.chn = Atom_instance.chain_identifier
        self.key = (self.num, self.name, self.chn)
        self.branched = 0
        self.missingSidechainAtoms = 0
        self.hasAlternateConformation = 0
        self.alt = {}
        self.atoms = {}
        self.atomList = [Atom_instance]
        atom_name = Atom_instance.atom_name.split()[0]
        if getattr(Atom_instance, "atom_alternate_location",
                   getattr(Atom_instance, "alternate_location", "")):
            alt_id = getattr(Atom_instance, "atom_alternate_location",
                             getattr(Atom_instance, "alternate_location", ""))
            self.alt.update({alt_id: {atom_name: Atom_instance}})
        self.atoms.update({atom_name: Atom_instance})
        self.ter = None
        self.bcom = None
        self.scom = None
        self.xtra_info = None

    def addatom(self, Atom_instance):
        atom_name = Atom_instance.atom_name.split()[0]
        self.atoms.update({atom_name: Atom_instance})
        self.atomList.append(Atom_instance)
        alt_id = getattr(Atom_instance, "atom_alternate_location",
                         getattr(Atom_instance, "alternate_location", ""))
        if alt_id:
            if alt_id in self.alt:
                self.alt[alt_id].update({atom_name: Atom_instance})
            else:
                self.alt.update({alt_id: {atom_name: Atom_instance}})

    def change_chain_identifier(self, chain):
        for atom_key in self.atoms:
            a = self.atoms[atom_key]
            a.chain_identifier = chain
        for alt_key in self.alt:
            conf = self.alt[alt_key]
            for atom_name in conf:
                conf[atom_name].chain_identifier = chain
        for atom in self.atomList:
            atom.chain_identifier = chain
        self.chn = chain
        self.key = (self.num, self.name, self.chn)

    def setAltConfBackboneAtoms(self):
        if self.hasAlternateConformation:
            for confKey in self.alt.keys():
                self.alt[confKey].update(self.atoms)

    def setsub(self):
        for atom in self.atoms.values():
            atom.residue = self
        for conf in self.alt.values():
            for atom in conf.values():
                atom.residue = self

    def set_sidechain_representations(self):
        self.ter = self.get_terminal_sidechain_atom()
        if not self.ter:
            self.ter = self.get_reduced_sidechain_atom()
        if not self.ter:
            try:
                self.ter = self.atoms["CA"]
            except:
                self.ter = self.get_reduced_backbone_atom()
        self.bcom = self.get_reduced_backbone_atom()
        self.scom = self.get_reduced_sidechain_atom()

    def get_atoms(self):
        atoms = []
        if self.alt:
            for conf in self.alt.values():
                for atom in conf.values():
                    atoms.append(atom.v)
        else:
            for atom in self.atoms.keys():
                atoms.append(self.atoms[atom].v)
        return atoms

    def get_backbone_atoms(self):
        atoms = []
        if self.alt:
            for conf in self.alt.values():
                for atom in conf.values():
                    if atom.atom_name in self.nonSidechainAtomNames:
                        atoms.append(atom)
        else:
            for atom in self.atoms.keys():
                if atom in self.nonSidechainAtomNames:
                    atoms.append(self.atoms[atom])
        return atoms

    def get_polar_atoms(self):
        atoms = []
        if self.alt:
            for conf in self.alt.values():
                for atom in conf.values():
                    if atom[0] in ["O", "N", "S"]:
                        atoms.append(atom.v)
        else:
            for atom in self.atoms.keys():
                if atom[0] in ["O", "N", "S"]:
                    atoms.append(self.atoms[atom].v)
        return atoms

    def get_sidechain_atoms(self):
        atoms = []
        if 0:
            for conf in self.alt.values():
                for atom in conf.values():
                    if atom.atom_name not in self.nonSidechainAtomNames:
                        if atom not in atoms:
                            atoms.append(atom)
                    elif atom.residue_name == "GLY":
                        try:
                            atoms.append(conf['CA'])
                        except:
                            continue
                if not atoms:
                    try:
                        atoms.append(conf['CA'])
                    except:
                        pass
        else:
            for atom in self.atoms.keys():
                if atom not in self.nonSidechainAtomNames:
                    atoms.append(self.atoms[atom])
                elif self.atoms[atom].residue_name == "GLY":
                    try:
                        atoms.append(self.atoms['CA'])
                    except:
                        continue
            if not atoms:
                try:
                    atoms.append(self.atoms['CA'])
                except:
                    pass
        return atoms

    def get_terminal_sidechain_atom(self):
        if self.name not in NONPOLAR + POLAR + IONIZABLE:
            return None
        terminal_atoms = {
            "GLY": "CA", "ALA": "CB", "ILE": "CD1", "MET": "CE", "CYS": "SG",
            "SER": "OG", "LYS": "NZ", "MLY": "NZ", "MLZ": "NZ", "M3L": "NZ"
        }
        multiple_terminal_atoms = {
            "THR": ("CG2", "OG1"),
            "LEU": ("CD1", "CD2"),
            "PHE": ("CZ",),
            "TYR": ("OH",),
            "PRO": ("CB", "CG", "CD"),
            "VAL": ("CG1", "CG2"),
            "ASN": ("OD1", "ND2"),
            "GLN": ("OE1", "NE2"),
            "TRP": ("CZ3", "CH2"),
            "ARG": ("NH1", "NH2"),
            "ASP": ("OD1", "OD2"),
            "GLU": ("OE1", "OE2"),
            "HIS": ("ND1", "NE2")
        }
        if self.name in terminal_atoms:
            try:
                terAtom = self.atoms[terminal_atoms[self.name]]
                ter = PseudoAtom()
                ter.x = terAtom.x
                ter.y = terAtom.y
                ter.z = terAtom.z
                ter.atom_serial = terAtom.atom_serial + 1
                ter.residue_sequence_number = terAtom.residue_sequence_number
                ter.atom_name = terAtom.atom_name
                ter.residue = self
                ter.residue_name = ter.residue.name
                ter.chain_identifier = terAtom.chain_identifier
                ter.reinitialize()
                self.ter = ter
                return ter
            except:
                self.missingSidechainAtoms = 1
                return None
        else:
            self.branched = 1
            n, x, y, z = 0.0, 0.0, 0.0, 0.0
            atomSerial = 0
            atomName = "TER"
            residueSequenceNumber = 0
            chain_identifier = ""
            for key in multiple_terminal_atoms[self.name]:
                try:
                    atom = self.atoms[key]
                    n += 1
                    x += atom.x
                    y += atom.y
                    z += atom.z
                    if atom.atom_serial > atomSerial:
                        atomSerial = atom.atom_serial
                        atomName = atom.atom_name
                        residueSequenceNumber = atom.residue_sequence_number
                        chain_identifier = atom.chain_identifier
                except:
                    self.missingSidechainAtoms = 1
                    return None
            ter = PseudoAtom()
            ter.atom_name = atomName
            ter.residue = self
            ter.atom_serial = atomSerial + 1
            ter.residue_sequence_number = residueSequenceNumber
            ter.residue_name = ter.residue.name
            ter.chain_identifier = chain_identifier
            ter.x = x / n
            ter.y = y / n
            ter.z = z / n
            ter.reinitialize()
            self.ter = ter
            return ter

    def get_preterminal_sidechain_atom(self):
        if self.alt:
            try:
                return self.atoms["CA"]
            except:
                return None
        if self.name not in NONPOLAR + POLAR + IONIZABLE:
            return None
        preterminal_atoms = {
            "ALA": "CA", "ARG": "CZ", "ASN": "CG", "ASP": "CG", "CYS": "CB",
            "GLN": "CD", "GLU": "CD", "HIS": "CG", "ILE": "CG1", "LEU": "CG",
            "LYS": "CE", "MET": "SD", "PHE": "CB", "PRO": "CA", "SER": "CB",
            "THR": "CB", "TRP": "CB", "TYR": "CZ", "VAL": "CB"
        }
        if self.name in preterminal_atoms:
            try:
                ter = self.atoms[preterminal_atoms[self.name]]
                ter.residue = self
                ter.reinitialize()
                return ter
            except:
                return None

    def get_reduced_backbone_atom(self):
        if self.alt:
            try:
                return self.atoms["CA"].v
            except:
                return None
        if self.name not in NONPOLAR + POLAR + IONIZABLE:
            return None
        atoms = self.get_backbone_atoms()
        n, x, y, z = 0.0, 0.0, 0.0, 0.0
        for atom in atoms:
            n += 1
            x += atom.x
            y += atom.y
            z += atom.z
        com = Atom()
        com.record_name = "ATOM "
        keys = sorted(self.atoms)
        com.atom_serial = self.atoms[keys[0]].atom_serial
        com.atom_name = "COM"
        com.residue_name = self.name
        com.chain_identifier = self.chn
        com.residue_sequence_number = self.num
        com.residue_key = self.key
        if n:
            com.x = x / n
            com.y = y / n
            com.z = z / n
        else:
            com.x = 0.0
            com.y = 0.0
            com.z = 0.0
        com.residue = Residue(com)
        com.residue.name = self.name
        com.residue.num = self.num
        com.residue.chn = self.chn
        com.reinitialize()
        com.residue = self
        return com.v

    def get_reduced_sidechain_atom(self):
        if self.alt:
            try:
                return self.atoms["CA"].v
            except:
                return None
        if self.name not in NONPOLAR + POLAR + IONIZABLE:
            return None
        atoms = self.get_sidechain_atoms()
        if not atoms:
            return None
        n, x, y, z = 0.0, 0.0, 0.0, 0.0
        for atom in atoms:
            n += 1
            x += atom.x
            y += atom.y
            z += atom.z
        com = Atom()
        com.record_name = "ATOM "
        keys = sorted(self.atoms)
        com.atom_serial = self.atoms[keys[0]].atom_serial
        com.atom_name = "COM"
        com.residue_name = self.name
        com.chain_identifier = self.chn
        com.residue_sequence_number = self.num
        com.residue_key = self.key
        if n:
            com.x = x / n
            com.y = y / n
            com.z = z / n
        else:
            com.x = 0.0
            com.y = 0.0
            com.z = 0.0
        com.residue = Residue(com)
        com.residue.name = self.name
        com.residue.num = self.num
        com.residue.chn = self.chn
        com.reinitialize()
        com.residue = self
        return com

    def __repr__(self):
        string = ""
        atoms = {}
        for atom in self.atoms.values():
            atoms[atom.atom_serial] = atom
        sortedKeys = sorted(atoms.keys())
        for sortedKey in sortedKeys:
            string += str(atoms[sortedKey])
        return string


class Chain:
    def __init__(self, pdbfilename, pdbcode, chain_identifier,
                 resdict, peptide=1, het_residues=None):
        self.pdbFileName = pdbfilename
        self.pdbCode = pdbcode
        self.chain_identifier = chain_identifier
        self.residues = resdict
        self.length = len(self.residues)
        self.gap = None
        self.peptide = peptide
        self.het_residues = het_residues
        if self.peptide:
            self.gap = self.maxgap()
        self.nonpolar_sc = {}
        self.nonpolar_polar_sc = {}
        self.nonpolar_ionizable_sc = {}
        self.polar_sc = {}
        self.polar_ionizable_sc = {}
        self.ionizable_sc = {}
        self.active_sc = {}
        for r in self.residues:
            if self.residues[r].name in NONPOLAR:
                self.nonpolar_sc.update({r: self.residues[r]})
            if self.residues[r].name in NONPOLAR + POLAR:
                self.nonpolar_polar_sc.update({r: self.residues[r]})
            if self.residues[r].name in NONPOLAR + IONIZABLE:
                self.nonpolar_ionizable_sc.update({r: self.residues[r]})
            if self.residues[r].name in POLAR:
                self.polar_sc.update({r: self.residues[r]})
            if self.residues[r].name in POLAR + IONIZABLE:
                self.polar_ionizable_sc.update({r: self.residues[r]})
            if self.residues[r].name in IONIZABLE:
                self.ionizable_sc.update({r: self.residues[r]})
            if self.residues[r].name in ACTIVE:
                self.active_sc.update({r: self.residues[r]})

    def maxgap(self):
        keys = sorted(self.residues)
        gap = 1
        resnum = self.residues[keys[0]].num
        for key in keys[1:]:
            current_gap = key[0] - resnum
            if current_gap > gap:
                gap = current_gap
            resnum = key[0]
        return gap

    def get_c_alpha_atoms(self, type="all", return_Vertex=0, return_Atom=1):
        atoms = []
        if type == "all":
            for res in self.residues.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "nonpolar":
            for res in self.nonpolar_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "nonpolar-polar":
            for res in self.nonpolar_polar_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "nonpolar-ionizable":
            for res in self.nonpolar_ionizable_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "polar":
            for res in self.polar_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "polar-ionizable":
            for res in self.polar_ionizable_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "ionizable":
            for res in self.ionizable_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        elif type == "active":
            for res in self.active_sc.values():
                try:
                    atoms.append(res.atoms['CA'])
                except:
                    continue
        if return_Vertex:
            vertices = []
            for atom in atoms:
                vertex = atom.v
                vertex.data = atom.Atom_LO
                vertices.append(vertex)
            return vertices
        else:
            return atoms

    def get_reduced_sidechain_atoms(self, type="all"):
        atoms = []
        if type == "all":
            for res in self.residues.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "nonpolar":
            for res in self.nonpolar_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "nonpolar-polar":
            for res in self.nonpolar_polar_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "nonpolar-ionizable":
            for res in self.nonpolar_ionizable_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "polar":
            for res in self.polar_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "polar-ionizable":
            for res in self.polar_ionizable_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "ionizable":
            for res in self.ionizable_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        elif type == "active":
            for res in self.active_sc.values():
                if res.scom:
                    atoms.append(res.scom)
        return atoms

    def get_terminal_sidechain_atoms(self, type="all",
                                     return_Vertex=0, return_Atom=1):
        atoms = []
        if type == "all":
            for res in self.residues.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "nonpolar":
            for res in self.nonpolar_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "nonpolar-polar":
            for res in self.nonpolar_polar_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "nonpolar-ionizable":
            for res in self.nonpolar_ionizable_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "polar":
            for res in self.polar_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "polar-ionizable":
            for res in self.polar_ionizable_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "ionizable":
            for res in self.ionizable_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        elif type == "active":
            for res in self.active_sc.values():
                if res.ter:
                    atoms.append(res.ter)
        if return_Vertex:
            vertices = []
            for atom in atoms:
                vertex = atom.v
                vertex.data = atom.Atom_LO
                vertices.append(vertex)
            return vertices
        else:
            return atoms

    def get_backbone_atoms(self, type="all"):
        atoms = []
        if type == "all":
            for res in self.residues.values():
                atoms += res.get_backbone_atoms()
        elif type == "nonpolar":
            for res in self.nonpolar_sc.values():
                atoms.append(res.get_backbone_atoms())
        elif type == "nonpolar-polar":
            for res in self.nonpolar_polar_sc.values():
                atoms.append(res.get_backbone_atoms())
        elif type == "nonpolar-ionizable":
            for res in self.nonpolar_ionizable_sc.values():
                atoms.append(res.get_backbone_atoms())
        elif type == "polar":
            for res in self.polar_sc.values():
                atoms.append(res.get_backbone_atoms())
        elif type == "polar-ionizable":
            for res in self.polar_ionizable_sc.values():
                atoms.append(res.get_backbone_atoms())
        elif type == "ionizable":
            for res in self.ionizable_sc.values():
                atoms.append(res.get_backbone_atoms())
        elif type == "active":
            for res in self.active_sc.values():
                atoms.append(res.get_backbone_atoms())
        return atoms

    def get_sidechain_atoms(self, type="all"):
        atoms = []
        if type == "all":
            for res in self.residues.values():
                atoms += res.get_sidechain_atoms()
        elif type == "nonpolar":
            for res in self.nonpolar_sc.values():
                atoms.append(res.get_sidechain_atoms())
        elif type == "nonpolar-polar":
            for res in self.nonpolar_polar_sc.values():
                atoms.append(res.get_sidechain_atoms())
        elif type == "nonpolar-ionizable":
            for res in self.nonpolar_ionizable_sc.values():
                atoms.append(res.get_sidechain_atoms())
        elif type == "polar":
            for res in self.polar_sc.values():
                atoms.append(res.get_sidechain_atoms())
        elif type == "polar-ionizable":
            for res in self.polar_ionizable_sc.values():
                atoms.append(res.get_sidechain_atoms())
        elif type == "ionizable":
            for res in self.ionizable_sc.values():
                atoms.append(res.get_sidechain_atoms())
        elif type == "active":
            for res in self.active_sc.values():
                atoms.append(res.get_sidechain_atoms())
        return atoms

    def general_information(self):
        string = ""
        string += "%-6s%s\n" % ("LCHN", str(self.length))
        string += "%-6s%s\n" % ("GAPS", str(self.gap))
        return string

class PDBfile:
    def __init__(self, pdbFilePath="", pdbFileName="", pdbFileAsString="",
                 twoCharacterChain=0, zip_status=0):
        self.zip = zip_status
        self.pdbFilePath = pdbFilePath
        self.pdbFileName = pdbFileName
        self.format = "pdb"
        self.pdbCode = "NONE"

        # Strip a compression suffix before looking at the extension, or
        # "8hs2.cif.gz" matches neither branch, falls through to the default,
        # and a CIF gets handed to the PDB parser.
        lower_name = pdbFileName.lower()
        if lower_name.endswith(".gz"):
            lower_name = lower_name[:-len(".gz")]
        if lower_name.endswith(".cif") or lower_name.endswith(".mmcif"):
            self.format = "cif"
            self.pdbCode = pdbFileName.split(".cif")[0].split(".mmcif")[0]
        elif lower_name.endswith(".pdb"):
            self.format = "pdb"
            self.pdbCode = pdbFileName.split(".pdb")[0]
        else:
            self.pdbCode = pdbFileName

        self.chain_selected_in_file_name, self.chains = 0, []
        if "-chains." in self.pdbCode:
            self.chain_selected_in_file_name = 1
            chain_string = self.pdbCode.split("-chains.")[1]
            for chain in chain_string.split("."):
                if chain:
                    self.chains.append(chain)

        self.atoms = {}
        self.hetatoms = {}
        self.residues, self.res_bychain = {}, {}
        self.het_residues, self.het_bychain = {}, {}
        self.isNmrStructure = 0
        self.twoCharacterChain = twoCharacterChain

        if self.format == "pdb":
            self._init_from_pdb(pdbFileAsString)
        else:
            self._init_from_cif(pdbFileAsString)

    # PDB path
    def _init_from_pdb(self, pdbFileAsString):
        self.rewritePdbFiles = []
        makeChainA = 0

        # Scan for NULL chain
        if self.zip:
            pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
            for line in pdbfile:
                line = line.decode()
                if line[0:4] == "ATOM" and line[21:22] == " ":
                    makeChainA = 1
                    self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
                    print("PDB file has NULL chain. Rename as chain A....")
                    break
        else:
            pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
            for line in pdbfile:
                if line[0:4] == "ATOM" and line[21:22] == " ":
                    makeChainA = 1
                    self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
                    print("PDB file has NULL chain. Rename as chain A....")
                    break

        # Scan for insertion codes
        if self.zip:
            pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
            for line in pdbfile:
                line = line.decode()
                if line[0:4] == "ATOM" and line[26:27] != " ":
                    self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
                    print("PDB file has residue insertion codes. Must renumber....")
                    break
        else:
            pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
            for line in pdbfile:
                if line[0:4] == "ATOM" and line[26:27] != " ":
                    self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
                    print("PDB file has residue insertion codes. Must renumber....")
                    break

        # Rewrite for NULL chain / insertion codes if needed
        if self.rewritePdbFiles:
            i = 0
            fixedPdbFileLines = ""
            if self.zip:
                pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
                for line in pdbfile:
                    line = line.decode()
                    if line[0:4] == "ATOM":
                        atom = Atom(pdbfileline=line)
                        if makeChainA:
                            atom.chain_identifier = "A"
                        if atom.atom_name == " N":
                            i += 1
                        atom.residue_sequence_number = i
                        atom.residue_insertion_code = ""
                        atom.reinitialize()
                        fixedPdbFileLines += str(atom)
                    else:
                        fixedPdbFileLines += line
                pdbFileRewrite = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "wb")
                pdbFileRewrite.write(fixedPdbFileLines.encode())
                pdbFileRewrite.close()
            else:
                pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
                for line in pdbfile:
                    if line[0:4] == "ATOM":
                        atom = Atom(pdbfileline=line)
                        if makeChainA:
                            atom.chain_identifier = "A"
                        if atom.atom_name == " N":
                            i += 1
                        atom.residue_sequence_number = i
                        atom.residue_insertion_code = ""
                        atom.reinitialize()
                        fixedPdbFileLines += str(atom)
                    else:
                        fixedPdbFileLines += line
                pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName, "w")
                pdbFileRewrite.write(fixedPdbFileLines)
                pdbFileRewrite.close()

        self.rewritePdbFiles = []

        # Scan for hex serials / residue numbers
        if self.zip:
            pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
            for line in pdbfile:
                line = line.decode()
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    atomInstance = Atom(pdbfileline=line)
                    if atomInstance.hex_atom_serial or atomInstance.hex_residue_sequence_number:
                        self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
                        print("PDB file has hexidecimal atom serial or residue numbers. Must renumber....")
                        print("By default, a single-character chain identifier, 5-digit atom serial, and 4-digit residue number is assumed.")
                        break
        else:
            pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
            for line in pdbfile:
                if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                    atomInstance = Atom(pdbfileline=line)
                    if atomInstance.hex_atom_serial or atomInstance.hex_residue_sequence_number:
                        self.rewritePdbFiles.append(self.pdbFilePath + self.pdbFileName)
                        print("PDB file has hexidecimal atom serial or residue numbers. Must renumber....")
                        print("By default, a single-character chain identifier, 5-digit atom serial, and 4-digit residue number is assumed.")
                        break

        # Rewrite for hex IDs if needed
        if self.rewritePdbFiles:
            fixedPdbFileLines = ""
            if self.zip:
                pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
                for line in pdbfile:
                    line = line.decode()
                    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        atom = Atom(pdbfileline=line)
                        fixedPdbFileLines += str(atom)
                    else:
                        fixedPdbFileLines += line
                pdbFileRewrite = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "wb")
                pdbFileRewrite.write(fixedPdbFileLines.encode())
                pdbFileRewrite.close()
            else:
                pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
                for line in pdbfile:
                    if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                        atom = Atom(pdbfileline=line)
                        fixedPdbFileLines += str(atom)
                    else:
                        fixedPdbFileLines += line
                pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName, "w")
                pdbFileRewrite.write(fixedPdbFileLines)
                pdbFileRewrite.close()

        # Scan for heteroatoms mislabelled as ATOM
        self.rewritePdbFilesHeteroAtoms = []
        if self.zip:
            pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
            for line in pdbfile:
                line = line.decode()
                if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
                    self.rewritePdbFilesHeteroAtoms.append(self.pdbFilePath + self.pdbFileName)
                    print("PDB file has incorrectly formatted heteroatoms. Must reformat....")
                    break
        else:
            pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
            for line in pdbfile:
                if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
                    self.rewritePdbFilesHeteroAtoms.append(self.pdbFilePath + self.pdbFileName)
                    print("PDB file has incorrectly formatted heteroatoms. Must reformat....")
                    break

        # Rewrite heteroatoms if needed
        if self.rewritePdbFilesHeteroAtoms:
            fixedPdbFileLines = ""
            if self.zip:
                pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
                for line in pdbfile:
                    line = line.decode()
                    if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
                        atom = Atom(pdbfileline=line)
                        atom.record_name = "HETATM"
                        atom.reinitialize()
                        fixedPdbFileLines += str(atom)
                    else:
                        fixedPdbFileLines += line
                pdbFileRewrite = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "wb")
                pdbFileRewrite.write(fixedPdbFileLines.encode())
                pdbFileRewrite.close()
            else:
                pdbfile = open(self.pdbFilePath + self.pdbFileName, 'r').readlines()
                for line in pdbfile:
                    if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
                        atom = Atom(pdbfileline=line)
                        atom.record_name = "HETATM"
                        atom.reinitialize()
                        fixedPdbFileLines += str(atom)
                    else:
                        fixedPdbFileLines += line
                pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName, "w")
                pdbFileRewrite.write(fixedPdbFileLines)
                pdbFileRewrite.close()

        # Final parse
        self.pdbfile = []
        if self.zip:
            self.pdbfile += gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb").readlines()
        else:
            self.pdbfile += open(self.pdbFilePath + self.pdbFileName, 'r').readlines()

        i, resnum, previousAtom = 0, 0, Atom()
        for line in self.pdbfile:
            if self.zip:
                line = line.decode()
            if "EXPDTA" in line:
                self.isNmrStructure = 1
            if "ENDMDL" in line and self.isNmrStructure:
                break
            if line[0:3] == "TER":
                i = 0
            if line[0:4] == "ATOM":
                atom = Atom(pdbfileline=line, twoCharacterChain=self.twoCharacterChain)
                if atom.chain_identifier == "NULL":
                    atom.chain_identifier = "A"
                    atom.reinitialize()
                previousAtom = atom
                if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
                    self.chains.append(atom.chain_identifier)
                self.atoms.update({atom.atom_serial: atom})
                if atom.residue_key in self.residues:
                    chain = atom.chain_identifier
                    residue = self.residues[atom.residue_key]
                    self.residues[atom.residue_key].addatom(atom)
                    residue.addatom(atom)
                    self.res_bychain[chain].update({atom.residue_key: residue})
                else:
                    chain = atom.chain_identifier
                    residue = Residue(atom)
                    self.residues.update({residue.key: residue})
                    if atom.chain_identifier not in self.res_bychain:
                        self.res_bychain.update({chain: {atom.residue_key: residue}})
                    else:
                        self.res_bychain[atom.chain_identifier].update({atom.residue_key: residue})
                if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
                    self.chains.append(atom.chain_identifier)
            if line[0:6] == "HETATM":
                atom = Atom(pdbfileline=line, twoCharacterChain=self.twoCharacterChain)
                self.hetatoms.update({atom.atom_serial: atom})
                if atom.residue_key in self.het_residues:
                    self.het_residues[atom.residue_key].addatom(atom)
                else:
                    residue = Residue(atom)
                    self.het_residues[residue.key] = residue
                    if atom.chain_identifier not in self.het_bychain:
                        self.het_bychain[atom.chain_identifier] = {}
                    self.het_bychain[atom.chain_identifier][atom.residue_key] = residue
                if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
                    self.chains.append(atom.chain_identifier)

        if self.residues:
            for residue in self.residues.values():
                residue.setAltConfBackboneAtoms()
                residue.setsub()
                residue.set_sidechain_representations()

        self.res_chains = {}
        for chain in self.res_bychain.keys():
            self.res_chains.update({
                chain: Chain(
                    self.pdbFileName,
                    self.pdbCode,
                    chain,
                    self.res_bychain[chain],
                    peptide=1,
                    het_residues=self.het_residues
                )
            })

        self.het_chains = {}
        for chain in self.het_bychain:
            self.het_chains.update({
                chain: Chain(
                    self.pdbFileName,
                    self.pdbCode,
                    chain,
                    self.het_bychain[chain],
                    peptide=0
                )
            })

    # CIF path
    def _init_from_cif(self, pdbFileAsString):
        self.pdbFileLines = []
        self.cif_options_all = {
            "_atom_site.group_PDB": False,
            "_atom_site.id": False,
            "_atom_site.type_symbol": False,
            "_atom_site.label_atom_id": False,
            "_atom_site.label_alt_id": False,
            "_atom_site.label_comp_id": False,
            "_atom_site.label_asym_id": False,
            "_atom_site.label_entity_id": False,
            "_atom_site.label_seq_id": False,
            "_atom_site.pdbx_PDB_ins_code": False,
            "_atom_site.Cartn_x": False,
            "_atom_site.Cartn_y": False,
            "_atom_site.Cartn_z": False,
            "_atom_site.occupancy": False,
            "_atom_site.B_iso_or_equiv": False,
            "_atom_site.Cartn_x_esd": False,
            "_atom_site.Cartn_y_esd": False,
            "_atom_site.Cartn_z_esd": False,
            "_atom_site.occupancy_esd": False,
            "_atom_site.B_iso_or_equiv_esd": False,
            "_atom_site.pdbx_formal_charge": False,
            "_atom_site.auth_seq_id": False,
            "_atom_site.auth_comp_id": False,
            "_atom_site.auth_asym_id": False,
            "_atom_site.auth_atom_id": False,
            "_atom_site.pdbx_PDB_model_num": False
        }

        # Read the file, recording the _atom_site column headers IN THE ORDER THE
        # FILE DECLARES THEM. The atom lines are matched to columns positionally,
        # so the order has to come from the file: cif_options_all is a fixed
        # 26-key superset, and any file that omits a column -- almost every
        # PDB-deposited mmCIF omits the five *_esd fields -- would otherwise have
        # every field after the first gap read out of the wrong column.
        column_order = []

        def note(candidate):
            if candidate.startswith("_atom_site."):
                if candidate in self.cif_options_all:
                    self.cif_options_all[candidate] = True
                if candidate not in column_order:
                    column_order.append(candidate)

        if self.zip and self.pdbFilePath and self.pdbFileName:
            f = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb")
            for line in f:
                line = line.decode()
                self.pdbFileLines.append(line)
                note(line.strip())
        elif pdbFileAsString:
            for line in pdbFileAsString.split("\n"):
                note(line.strip())
                self.pdbFileLines.append(line + "\n")
        elif self.pdbFilePath and self.pdbFileName:
            f = open(self.pdbFilePath + self.pdbFileName, "r")
            for line in f:
                self.pdbFileLines.append(line)
                note(line.strip())

        # What the atom parser walks. Falls back to the fixed order only when a
        # file declares no headers at all.
        self.cif_columns = {k: True for k in column_order} or self.cif_options_all

        # Parse atoms and hetatoms
        for line in self.pdbFileLines:
            if line[0:4] == "ATOM":
                atom = Atom(cif_line=line, cif_options=self.cif_columns)
                self.atoms[atom.atom_serial] = atom
                if atom.residue_key in self.residues:
                    self.residues[atom.residue_key].addatom(atom)
                else:
                    residue = Residue(atom)
                    self.residues[residue.key] = residue
                    if atom.chain_identifier not in self.res_bychain:
                        self.res_bychain[atom.chain_identifier] = {}
                    self.res_bychain[atom.chain_identifier][atom.residue_key] = residue
                if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
                    self.chains.append(atom.chain_identifier)

            if line[0:6] == "HETATM":
                atom = Atom(cif_line=line, cif_options=self.cif_columns)
                self.hetatoms[atom.atom_serial] = atom
                if atom.residue_key in self.het_residues:
                    self.het_residues[atom.residue_key].addatom(atom)
                else:
                    residue = Residue(atom)
                    self.het_residues[residue.key] = residue
                    if atom.chain_identifier not in self.het_bychain:
                        self.het_bychain[atom.chain_identifier] = {}
                    self.het_bychain[atom.chain_identifier][atom.residue_key] = residue
                if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
                    self.chains.append(atom.chain_identifier)

        if self.residues:
            for residue in self.residues.values():
                residue.setAltConfBackboneAtoms()
                residue.setsub()
                residue.set_sidechain_representations()

        self.res_chains = {}
        for chain in self.res_bychain.keys():
            self.res_chains.update({
                chain: Chain(
                    self.pdbFileName,
                    self.pdbCode,
                    chain,
                    self.res_bychain[chain],
                    peptide=1,
                    het_residues=self.het_residues
                )
            })

        self.het_chains = {}
        for chain in self.het_bychain:
            self.het_chains.update({
                chain: Chain(
                    self.pdbFileName,
                    self.pdbCode,
                    chain,
                    self.het_bychain[chain],
                    peptide=0
                )
            })

    def res_atom_hash(self):
        string = ""
        keys = list(self.atoms.keys())
        keys.sort()
        for key in keys:
            atom = self.atoms[key]
            string += "%-6s|%s" % (str(atom.atom_serial), str(atom))
        return string

    def het_atom_hash(self):
        string = ""
        keys = list(self.hetatoms.keys())
        keys.sort()
        for key in keys:
            atom = self.hetatoms[key]
            string += "%-6s|%s" % (str(atom.atom_serial), str(atom))
        return string

    def convert_to_pdb_format(self):
        pdb_formatted_string = ""
        for atom_serial in self.atoms:
            atom = self.atoms[atom_serial]
            pdb_formatted_string += str(atom)
        for atom_serial in self.hetatoms:
            atom = self.hetatoms[atom_serial]
            pdb_formatted_string += str(atom)
        return pdb_formatted_string


# class PDBfile:
#     def __init__(self, pdbFilePath="", pdbFileName="", pdbFileAsString="",
#                  twoCharacterChain=0, zip_status=0):
#         self.zip = zip_status
#         self.pdbFilePath = pdbFilePath
#         self.pdbFileName = pdbFileName
#         self.format = "pdb"
#         self.pdbCode = "NONE"
#         lower_name = pdbFileName.lower()
#         if lower_name.endswith(".cif") or lower_name.endswith(".mmcif"):
#             self.format = "cif"
#             self.pdbCode = pdbFileName.split(".cif")[0].split(".mmcif")[0]
#         elif lower_name.endswith(".pdb"):
#             self.format = "pdb"
#             self.pdbCode = pdbFileName.split(".pdb")[0]
#         else:
#             self.pdbCode = pdbFileName
#         self.chain_selected_in_file_name, self.chains = 0, []
#         if "-chains." in self.pdbCode:
#             self.chain_selected_in_file_name = 1
#             chain_string = self.pdbCode.split("-chains.")[1]
#             for chain in chain_string.split("."):
#                 if chain:
#                     self.chains.append(chain)
#         self.atoms = {}
#         self.hetatoms = {}
#         self.residues, self.res_bychain = {}, {}
#         self.het_residues, self.het_bychain = {}, {}
#         self.isNmrStructure = 0
#         self.twoCharacterChain = twoCharacterChain
#         if self.format == "pdb":
#             self._init_from_pdb(pdbFileAsString)
#         else:
#             self._init_from_cif(pdbFileAsString)

#     # PDB path from original pdbFile_pdb.py
#     def _init_from_pdb(self, pdbFileAsString):
#         self.rewritePdbFiles = []
#         makeChainA = 0
#         if self.zip:
#             pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                     "rb").readlines()
#             for line in pdbfile:
#                 line = line.decode()
#                 if line[0:4] == "ATOM" and line[21:22] == " ":
#                     makeChainA = 1
#                     self.rewritePdbFiles.append(self.pdbFilePath +
#                                                 self.pdbFileName)
#                     print("PDB file has NULL chain. Rename as chain A....")
#                     break
#         else:
#             pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                            'r').readlines()
#             for line in pdbfile:
#                 if line[0:4] == "ATOM" and line[21:22] == " ":
#                     makeChainA = 1
#                     self.rewritePdbFiles.append(self.pdbFilePath +
#                                                 self.pdbFileName)
#                     print("PDB file has NULL chain. Rename as chain A....")
#                     break
#         if self.zip:
#             pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                     "rb").readlines()
#             for line in pdbfile:
#                 line = line.decode()
#                 if line[0:4] == "ATOM" and line[26:27] != " ":
#                     self.rewritePdbFiles.append(self.pdbFilePath +
#                                                 self.pdbFileName)
#                     print("PDB file has residue insertion codes. Must renumber....")
#                     break
#         else:
#             pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                            'r').readlines()
#             for line in pdbfile:
#                 if line[0:4] == "ATOM" and line[26:27] != " ":
#                     self.rewritePdbFiles.append(self.pdbFilePath +
#                                                 self.pdbFileName)
#                     print("PDB file has residue insertion codes. Must renumber....")
#                     break
#         if self.rewritePdbFiles:
#             i = 0
#             fixedPdbFileLines = ""
#             if self.zip:
#                 pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                         "rb").readlines()
#                 for line in pdbfile:
#                     line = line.decode()
#                     if line[0:4] == "ATOM":
#                         atom = Atom(line)
#                         if makeChainA:
#                             atom.chain_identifier = "A"
#                         if atom.atom_name == " N":
#                             i += 1
#                         atom.residue_sequence_number = i
#                         atom.residue_insertion_code = ""
#                         atom.reinitialize()
#                         fixedPdbFileLines += str(atom)
#                     else:
#                         fixedPdbFileLines += line
#                 pdbFileRewrite = gzip.GzipFile(self.pdbFilePath +
#                                                self.pdbFileName, "wb")
#                 pdbFileRewrite.write(fixedPdbFileLines.encode())
#                 pdbFileRewrite.close()
#             else:
#                 pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                                'r').readlines()
#                 for line in pdbfile:
#                     if line[0:4] == "ATOM":
#                         atom = Atom(line)
#                         if makeChainA:
#                             atom.chain_identifier = "A"
#                         if atom.atom_name == " N":
#                             i += 1
#                         atom.residue_sequence_number = i
#                         atom.residue_insertion_code = ""
#                         atom.reinitialize()
#                         fixedPdbFileLines += str(atom)
#                     else:
#                         fixedPdbFileLines += line
#                 pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName,
#                                       "w")
#                 pdbFileRewrite.write(fixedPdbFileLines)
#                 pdbFileRewrite.close()
#         self.rewritePdbFiles = []
#         if self.zip:
#             pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                     "rb").readlines()
#             for line in pdbfile:
#                 line = line.decode()
#                 if line[0:4] == "ATOM" or line[0:6] == "HETATM":
#                     atomInstance = Atom(line)
#                     if atomInstance.hex_atom_serial or atomInstance.hex_residue_sequence_number:
#                         self.rewritePdbFiles.append(self.pdbFilePath +
#                                                     self.pdbFileName)
#                         print("PDB file has hexidecimal atom serial or residue numbers. Must renumber....")
#                         print("By default, a single-character chain identifier, 5-digit atom serial, and 4-digit residue number is assumed.")
#                         break
#         else:
#             pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                            'r').readlines()
#             for line in pdbfile:
#                 if line[0:4] == "ATOM" or line[0:6] == "HETATM":
#                     atomInstance = Atom(line)
#                     if atomInstance.hex_atom_serial or atomInstance.hex_residue_sequence_number:
#                         self.rewritePdbFiles.append(self.pdbFilePath +
#                                                     self.pdbFileName)
#                         print("PDB file has hexidecimal atom serial or residue numbers. Must renumber....")
#                         print("By default, a single-character chain identifier, 5-digit atom serial, and 4-digit residue number is assumed.")
#                         break
#         if self.rewritePdbFiles:
#             fixedPdbFileLines = ""
#             if self.zip:
#                 pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                         "rb").readlines()
#                 for line in pdbfile:
#                     line = line.decode()
#                     if line[0:4] == "ATOM" or line[0:6] == "HETATM":
#                         atom = Atom(line)
#                         fixedPdbFileLines += str(atom)
#                     else:
#                         fixedPdbFileLines += line
#                 pdbFileRewrite = gzip.GzipFile(self.pdbFilePath +
#                                                self.pdbFileName, "wb")
#                 pdbFileRewrite.write(fixedPdbFileLines.encode())
#                 pdbFileRewrite.close()
#             else:
#                 pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                                'r').readlines()
#                 for line in pdbfile:
#                     if line[0:4] == "ATOM" or line[0:6] == "HETATM":
#                         atom = Atom(line)
#                         fixedPdbFileLines += str(atom)
#                     else:
#                         fixedPdbFileLines += line
#                 pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName,
#                                       "w")
#                 pdbFileRewrite.write(fixedPdbFileLines)
#                 pdbFileRewrite.close()
#         self.rewritePdbFilesHeteroAtoms = []
#         if self.zip:
#             pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                     "rb").readlines()
#             for line in pdbfile:
#                 line = line.decode()
#                 if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
#                     self.rewritePdbFilesHeteroAtoms.append(self.pdbFilePath +
#                                                            self.pdbFileName)
#                     print("PDB file has incorrectly formatted heteroatoms. Must reformat....")
#                     break
#         else:
#             pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                            'r').readlines()
#             for line in pdbfile:
#                 if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
#                     self.rewritePdbFilesHeteroAtoms.append(self.pdbFilePath +
#                                                            self.pdbFileName)
#                     print("PDB file has incorrectly formatted heteroatoms. Must reformat....")
#                     break
#         if self.rewritePdbFilesHeteroAtoms:
#             fixedPdbFileLines = ""
#             if self.zip:
#                 pdbfile = gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                         "rb").readlines()
#                 for line in pdbfile:
#                     line = line.decode()
#                     if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
#                         atom = Atom(line)
#                         atom.record_name = "HETATM"
#                         atom.reinitialize()
#                         fixedPdbFileLines += str(atom)
#                     else:
#                         fixedPdbFileLines += line
#                 pdbFileRewrite = gzip.GzipFile(self.pdbFilePath +
#                                                self.pdbFileName, "wb")
#                 pdbFileRewrite.write(fixedPdbFileLines.encode())
#                 pdbFileRewrite.close()
#             else:
#                 pdbfile = open(self.pdbFilePath + self.pdbFileName,
#                                'r').readlines()
#                 for line in pdbfile:
#                     if line[0:4] == "ATOM" and line[17:20] not in ALL_SIDECHAINS:
#                         atom = Atom(line)
#                         atom.record_name = "HETATM"
#                         atom.reinitialize()
#                         fixedPdbFileLines += str(atom)
#                     else:
#                         fixedPdbFileLines += line
#                 pdbFileRewrite = open(self.pdbFilePath + self.pdbFileName,
#                                       "w")
#                 pdbFileRewrite.write(fixedPdbFileLines)
#                 pdbFileRewrite.close()
#         self.pdbfile = []
#         if self.zip:
#             self.pdbfile += gzip.GzipFile(self.pdbFilePath + self.pdbFileName,
#                                           "rb").readlines()
#         else:
#             self.pdbfile += open(self.pdbFilePath + self.pdbFileName,
#                                  'r').readlines()
#         i, resnum, previousAtom = 0, 0, Atom()
#         for line in self.pdbfile:
#             if self.zip:
#                 line = line.decode()
#             if "EXPDTA" in line:
#                 self.isNmrStructure = 1
#             if "ENDMDL" in line and self.isNmrStructure:
#                 break
#             if line[0:3] == "TER":
#                 i = 0
#             if line[0:4] == "ATOM":
#                 atom = Atom(line, twoCharacterChain=self.twoCharacterChain)
#                 if atom.chain_identifier == "NULL":
#                     atom.chain_identifier = "A"
#                     atom.reinitialize()
#                 previousAtom = atom
#                 if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
#                     self.chains.append(atom.chain_identifier)
#                 self.atoms.update({atom.atom_serial: atom})
#                 if atom.residue_key in self.residues:
#                     chain = atom.chain_identifier
#                     residue = self.residues[atom.residue_key]
#                     self.residues[atom.residue_key].addatom(atom)
#                     residue.addatom(atom)
#                     self.res_bychain[chain].update({atom.residue_key: residue})
#                 else:
#                     chain = atom.chain_identifier
#                     residue = Residue(atom)
#                     self.residues.update({residue.key: residue})
#                     if atom.chain_identifier not in self.res_bychain:
#                         self.res_bychain.update({chain: {atom.residue_key: residue}})
#                     else:
#                         self.res_bychain[atom.chain_identifier].update({atom.residue_key: residue})
#                 if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
#                     self.chains.append(atom.chain_identifier)
#             if line[0:6] == "HETATM":
#                 atom = Atom(line, twoCharacterChain=self.twoCharacterChain)
#                 self.hetatoms.update({atom.atom_serial: atom})
#                 if atom.residue_key in self.het_residues:
#                     self.het_residues[atom.residue_key].addatom(atom)
#                 else:
#                     residue = Residue(atom)
#                     self.het_residues[residue.key] = residue
#                     if atom.chain_identifier not in self.het_bychain:
#                         self.het_bychain[atom.chain_identifier] = {}
#                     self.het_bychain[atom.chain_identifier][atom.residue_key] = residue
#                 if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
#                     self.chains.append(atom.chain_identifier)
#         if self.residues:
#             for residue in self.residues.values():
#                 residue.setAltConfBackboneAtoms()
#                 residue.setsub()
#                 residue.set_sidechain_representations()
#         self.res_chains = {}
#         for chain in self.res_bychain.keys():
#             self.res_chains.update({
#                 chain: Chain(
#                     self.pdbFileName,
#                     self.pdbCode,
#                     chain,
#                     self.res_bychain[chain],
#                     peptide=1,
#                     het_residues=self.het_residues
#                 )
#             })
#         self.het_chains = {}
#         for chain in self.het_bychain:
#             self.het_chains.update({
#                 chain: Chain(
#                     self.pdbFileName,
#                     self.pdbCode,
#                     chain,
#                     self.het_bychain[chain],
#                     peptide=0
#                 )
#             })

#     # CIF path from original pdbFile_cif.py
#     def _init_from_cif(self, pdbFileAsString):
#         self.pdbFileLines = []
#         self.cif_options_all = {
#             "_atom_site.group_PDB": False,
#             "_atom_site.id": False,
#             "_atom_site.type_symbol": False,
#             "_atom_site.label_atom_id": False,
#             "_atom_site.label_alt_id": False,
#             "_atom_site.label_comp_id": False,
#             "_atom_site.label_asym_id": False,
#             "_atom_site.label_entity_id": False,
#             "_atom_site.label_seq_id": False,
#             "_atom_site.pdbx_PDB_ins_code": False,
#             "_atom_site.Cartn_x": False,
#             "_atom_site.Cartn_y": False,
#             "_atom_site.Cartn_z": False,
#             "_atom_site.occupancy": False,
#             "_atom_site.B_iso_or_equiv": False,
#             "_atom_site.Cartn_x_esd": False,
#             "_atom_site.Cartn_y_esd": False,
#             "_atom_site.Cartn_z_esd": False,
#             "_atom_site.occupancy_esd": False,
#             "_atom_site.B_iso_or_equiv_esd": False,
#             "_atom_site.pdbx_formal_charge": False,
#             "_atom_site.auth_seq_id": False,
#             "_atom_site.auth_comp_id": False,
#             "_atom_site.auth_asym_id": False,
#             "_atom_site.auth_atom_id": False,
#             "_atom_site.pdbx_PDB_model_num": False
#         }
#         if self.zip and self.pdbFilePath and self.pdbFileName:
#             f = gzip.GzipFile(self.pdbFilePath + self.pdbFileName, "rb")
#             for line in f:
#                 line = line.decode()
#                 self.pdbFileLines.append(line)
#                 possible_cif_option = line.strip()
#                 if possible_cif_option in self.cif_options_all:
#                     self.cif_options_all[possible_cif_option] = True
#         elif pdbFileAsString:
#             for line in pdbFileAsString.split("\n"):
#                 possible_cif_option = line
#                 if possible_cif_option in self.cif_options_all:
#                     self.cif_options_all[possible_cif_option] = True
#                 self.pdbFileLines.append(line + "\n")
#         elif self.pdbFilePath and self.pdbFileName:
#             f = open(self.pdbFilePath + self.pdbFileName, "r")
#             for line in f:
#                 self.pdbFileLines.append(line)
#                 possible_cif_option = line.strip()
#                 if possible_cif_option in self.cif_options_all:
#                     self.cif_options_all[possible_cif_option] = True
#         for line in self.pdbFileLines:
#             if line[0:4] == "ATOM":
#                 atom = CIFAtom(line, self.cif_options_all)
#                 self.atoms[atom.atom_serial] = atom
#                 if atom.residue_key in self.residues:
#                     self.residues[atom.residue_key].addatom(atom)
#                 else:
#                     residue = Residue(atom)
#                     self.residues[residue.key] = residue
#                     if atom.chain_identifier not in self.res_bychain:
#                         self.res_bychain[atom.chain_identifier] = {}
#                     self.res_bychain[atom.chain_identifier][atom.residue_key] = residue
#                 if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
#                     self.chains.append(atom.chain_identifier)
#             if line[0:6] == "HETATM":
#                 atom = CIFAtom(line, self.cif_options_all)
#                 self.atoms[atom.atom_serial] = atom
#                 if atom.residue_key in self.residues:
#                     self.het_residues[atom.residue_key].addatom(atom)
#                 else:
#                     residue = Residue(atom)
#                     self.het_residues[residue.key] = residue
#                     if atom.chain_identifier not in self.het_bychain:
#                         self.het_bychain[atom.chain_identifier] = {}
#                     self.het_bychain[atom.chain_identifier][atom.residue_key] = residue
#                 if not self.chain_selected_in_file_name and atom.chain_identifier not in self.chains:
#                     self.chains.append(atom.chain_identifier)
#         if self.residues:
#             for residue in self.residues.values():
#                 residue.setAltConfBackboneAtoms()
#                 residue.setsub()
#                 residue.set_sidechain_representations()
#         self.res_chains = {}
#         for chain in self.res_bychain.keys():
#             self.res_chains.update({
#                 chain: Chain(
#                     self.pdbFileName,
#                     self.pdbCode,
#                     chain,
#                     self.res_bychain[chain],
#                     peptide=1,
#                     het_residues=self.het_residues
#                 )
#             })
#         self.het_chains = {}
#         for chain in self.het_bychain:
#             self.het_chains.update({
#                 chain: Chain(
#                     self.pdbFileName,
#                     self.pdbCode,
#                     chain,
#                     self.het_bychain[chain],
#                     peptide=0
#                 )
#             })

#     def res_atom_hash(self):
#         string = ""
#         keys = list(self.atoms.keys())
#         keys.sort()
#         for key in keys:
#             atom = self.atoms[key]
#             string += "%-6s|%s" % (str(atom.atom_serial), str(atom))
#         return string

#     def het_atom_hash(self):
#         string = ""
#         keys = list(self.hetatoms.keys())
#         keys.sort()
#         for key in keys:
#             atom = self.hetatoms[key]
#             string += "%-6s|%s" % (str(atom.atom_serial), str(atom))
#         return string

#     def convert_to_pdb_format(self):
#         pdb_formatted_string = ""
#         for atom_serial in self.atoms:
#             atom = self.atoms[atom_serial]
#             if isinstance(atom, CIFAtom):
#                 pdb_formatted_string += atom.get_pdb_format()
#             else:
#                 pdb_formatted_string += str(atom)
#         for atom_serial in self.hetatoms:
#             atom = self.hetatoms[atom_serial]
#             if isinstance(atom, CIFAtom):
#                 pdb_formatted_string += atom.get_pdb_format()
#             else:
#                 pdb_formatted_string += str(atom)
#         return pdb_formatted_string
