from typing import List, Dict, Union
from rdkit import Chem
from rdkit.Chem import AllChem

class Params:
    ordering = ['NAME', 'IO_STRING', 'TYPE', 'AA', '##' 'ATOM', 'BOND_TYPE', 'CHI', 'CONNECT', 'NBR_ATOM',
                'NBR_RADIUS', 'ICOOR_INTERNAL', 'PDB_ROTAMERS']
    prevent_newlines = True

    def __init__(self, mol :Chem.Mol,
                       resn:Union[str, None]=None,
                       chain:Union[str, None]=None,
                       mtype='POLYMER',
                       aa='UNK'):
        self.mol = mol
        self.resn = resn
        self.chain = chain
        self.type = mtype
        self.aa = aa
        self.root_id = 0
        self.parts = self.make_parts()

    def make_parts(self) -> Dict[List[str]]:
        """
        ``.parts`` is a dictionary of lists which stores the parts of the params file.
        The key order is dictated by ``.ordering`` class attribute.
        '##' key will be placed after AA but before ATOM.
        """
        parts = {k: [] for k in self.ordering}
        parts['NAME'].append(self.resn)
        parts['IO_STRING'].append(f'{self.resn} {self.chain}')
        parts['TYPE'] = self.mtype
        parts['AA'].append(self.aa)
        if self.mol.HasProp("_Name"):
            parts['##'].append(self.mol.GetProp("_Name"))
        for atom in self.mol.GetAtoms():
            d = self.make_atom_descriptor_line(atom)
            parts['ATOM'].append(d)
        parts['BONDS'].extend(self.make_bond_lines())
        return parts

    def write(self, file_name):
        with open(file_name,'w') as w:
            w.write(str(self))

    def __str__(self):
        text = ''
        for section in self.ordering:
            for line in self.parts[section]:
                if self.prevent_newlines:
                    sline = line.replace('\n', '')
                else:
                    sline = line
                text += f'{section} {sline}\n'
        return text

    def make_bond_lines(self) -> List[str]:
        # Crawl down from atom root_idx
        visited_bonds = []
        root = self.mol.GetAtomWithIdx(self.root_idx)
        bonds = root.GetBonds()
        for bond in bonds:
            if bond.GetIdx() in visited_bonds:
                continue
            else:
                visited_bonds.append(bond.GetIdx())

        # self.parts['BOND'].append(d)

    def make_atom_descriptor_line(self, atom) -> str:
        if not self.mol.HasProp('_GasteigerCharge'):
            AllChem.ComputeGasteigerCharges(self.mol)
        descriptors = {'name': self.correct_name(atom),
                       'type': self.get_atom_type(atom),
                       'type2': self.get_atom_type2(atom),
                       'charge': self.mol.GetDoubleProp('_GasteigerCharge')}
        return '{name} {type} {type2} {charge:.2f}'.format_map(descriptors)

    def correct_name(self, atom) -> str:
        info = atom.GetPDBResidueInfo()
        if info is not None:
            return info.GetName()
        else:
            return atom.GetSymbol()

    def get_atom_type(self, atom):
        return '1234'

    def get_atom_type2(self, atom):
        return '1234'
