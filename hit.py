from typing import Union, List, Dict
from covalent_dock import pymol2, pyrosetta

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdmolfiles, rdFMCS, Draw

import molfile_to_params
import os, re, shutil, random
from collections import defaultdict
from warnings import warn

class Hit:
    hits_path = '/well/brc/matteo/Mpro'
    work_path = 'templates'
    if not os.path.exists(work_path):
        os.mkdir(work_path)
    relax_cycles = 15

    def __init__(self, name):
        self._covalent_atomname = None
        self._name_map = None
        self._idx2name = None
        self._name2idx = None
        self.name = name
        self.mol_file = os.path.join(self.hits_path, f'Mpro-{name}_0', f'Mpro-{name}_0.mol')
        self.pdb_file = os.path.join(self.hits_path, f'Mpro-{name}_0', f'Mpro-{name}_0.pdb')
        self.mol2_file = os.path.join(self.hits_path, f'Mpro-{name}_0', f'Mpro-{name}_0.mol2')
        self.bound_file = os.path.join(self.hits_path, f'Mpro-{name}_0', f'Mpro-{name}_0_bound.pdb')
        self.relaxbound_file = os.path.join(self.work_path, f'Mpro-{name}_0_bound.r.pdb')
        self.apo_file = os.path.join(self.hits_path, f'Mpro-{name}_0', f'Mpro-{name}_0_apo.pdb')
        self.params_file = os.path.join(self.work_path, f'{self.name}.params')
        self.con_file = os.path.join(self.work_path, f'{self.name}.con')
        self.mol = self.load()

    def is_covalent(self):
        f = open(self.bound_file).read()
        if 'LINK' in f:
            return True
        else:
            return False

    @property #cached
    def covalent_atomname(self) -> Union[str, None]:
        if self._covalent_atomname is None:
            f = open(self.bound_file).read()
            if 'LINK' in f:
                self._covalent_atomname = ' '+re.search('LINKR?\s+(\S+)(.*)', f).group(1).rjust(3)
            else:
                raise ValueError('Why do you think its covalent?')
        return self._covalent_atomname

    @property
    def covalent_idx(self) -> int:
        if self.covalent_atomname in self.name2idx:
            return self.name2idx[self.covalent_atomname]
        elif len(self.covalent_atomname) != 4:
            return self.name2idx[' '+self.covalent_atomname.rjust(3)]
        elif self.covalent_atomname.strip() in self.name2idx:
            return self.name2idx[self.covalent_atomname.strip()]
        elif self.covalent_atomname.strip() in [x.strip() for x in self.name2idx]:
            print(f'WEIRD ATOM NAME SPACING: >{self.covalent_atomname}<')
            return {k.strip(): v for k, v in self.name2idx.items()}[self.covalent_atomname.strip()]
        else:
            raise ValueError(f'What is >{self.covalent_atomname}<')


    def load(self):
        # bond order and atom names. Why not?
        mol = Chem.MolFromMolFile(self.mol_file)
        pdb = Chem.MolFromPDBFile(self.pdb_file)
        for ma, pa in zip(mol.GetAtoms(), pdb.GetAtoms()):
            ma.SetMonomerInfo(pa.GetPDBResidueInfo())
        AllChem.ComputeGasteigerCharges(mol)
        Chem.AddHs(mol)
        return mol

    @property #cached
    def idx2name(self) -> Dict[int, str]:
        if self._idx2name is None:
            self._idx2name = {atom.GetIdx(): atom.GetPDBResidueInfo().GetName() for atom in self.mol.getAtoms()}
        return self._idx2name

    @property  # cached
    def name2idx(self) -> Dict[str, int]:
        if self._name2idx is None:
            self._name2idx = {atom.GetPDBResidueInfo().GetName(): atom.GetIdx() for atom in self.mol.GetAtoms()}
        return self._name2idx

    def get_icoor_from_ref(self) -> str:
        thio = Chem.MolFromSmiles('CCCS')
        thio = Chem.AddHs(thio)
        Chem.GetSSSR(thio)
        AllChem.EmbedMolecule(thio)
        AllChem.UFFOptimizeMolecule(thio, maxIters=2000)
        conf = thio.GetConformer()
        distance = Chem.rdMolTransforms.GetBondLength(conf, 2, 3)
        angle = Chem.rdMolTransforms.GetAngleDeg(conf,1, 2, 3)
        torsion = Chem.rdMolTransforms.GetDihedralDeg(conf, 0, 1, 2, 3)
        n = [atom for atom in self.mol.GetAtomWithIdx(self.covalent_idx).GetNeighbors() if atom.GetSymbol != 'H'][0]
        m = [atom for atom in n.GetNeighbors() if atom.GetSymbol != 'H' and atom.GetIdx != self.covalent_idx][0]
        return f'ICOOR_INTERNAL   CONN1  {torsion:.3f} {180 - angle:.3f} {distance:.3f} ' + \
               f'{self.covalent_atomname} {n.GetPDBResidueInfo().GetName()} {m.GetPDBResidueInfo().GetName()}\n'

    def relax(self):
        if os.path.exists(self.relaxbound_file):
            warn('Preexisting relaxed')
            return None
        self.parameterise()
        self.create_constraint_file()
        pose = self.load_pose()
        setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
        setup.constraint_file(self.con_file)
        setup.apply(pose)
        self.do_relax(pose)
        pose.dump_pdb(self.relaxbound_file)

    def create_constraint_file(self):
        with open(self.con_file,'w') as w:
            w.write('AtomPair SG 145A NE2 41A HARMONIC 3.8 0.2\n')
            if self.is_covalent():
                w.write(f'AtomPair SG 145A {self.covalent_atomname.strip()} 1B HARMONIC 1.8 0.2\n')
                w.write(f'Angle CB 145A SG 145A {self.covalent_atomname.strip()} 1B HARMONIC 1.71 0.35\n')

    def load_pose(self):
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(self.bound_file)
            pymol.cmd.remove('solvent')
            pymol.cmd.remove('not polymer and not resn LIG')  # should run it through P Curran\s list of artefacts!
            pymol.cmd.alter('resn LIG', 'chain="B"')
            pymol.cmd.alter('resn LIG', 'resi="1"')
            pymol.cmd.sort()
            holo = pymol.cmd.get_pdbstr('*')
        if self.is_covalent():
            holo =  f'LINK         SG  CYS A 145                {self.covalent_atomname} LIG B   1     1555   1555  1.8\n' + holo
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([self.params_file])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, holo)
        ## Force HIE
        r = pose.pdb_info().pdb2pose(res=41, chain='A')
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res='HIS').apply(pose)
        return pose

    def do_relax(self, pose):
        scorefxn = pyrosetta.get_fa_scorefxn()
        score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        scorefxn.set_weight(score_manager.score_type_from_name("atom_pair_constraint"), 20)
        scorefxn.set_weight(score_manager.score_type_from_name("angle_constraint"), 20)
        scorefxn.set_weight(score_manager.score_type_from_name("coordinate_constraint"), 1.0)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, self.relax_cycles)
        relax.constrain_relax_to_start_coords(True)
        relax.constrain_coords(True)
        relax.apply(pose)

    def parameterise(self):
        #Chem.MolToPDBFile(self.mol, flavor=0, removeHs=False)
        #os.system(f'obabel {self.pdb_file} -O {self.pdb_file.replace(".pdb",".mol2")} -h')
        ## make a params.
        molfile_to_params.run(self.mol2_file,
                              conformers_in_one_file=True,
                              name='LIG',
                              keep_names=True,
                              amino_acid=None,
                              chain='B',
                              pdb=os.path.join(self.work_path, self.name),
                              clobber=True)
        ### split into parts
        allparts = defaultdict(list)
        with open(self.params_file) as r:
            for line in r:
                parts = line.split()
                allparts[parts[0]].append(line)
        ### fix partial charges.
        # for li in range(len(allparts['ATOM'])):
        #     line = allparts['ATOM'][li]
        #     rex = re.match('^ATOM\s+(\w+)(\s+\w+\s+\w+\s+)([.\d\-])*', line)
        #     name = rex.group(1)
        #     print(self.name2idx)
        #     gc = self.mol.GetAtomWithIdx(self.name2idx[name]).GetDoubleProp('_GasteigerCharge')
        #     allparts['ATOM'][li] = line.replace(rex.group(2) + rex.group(3), f'{gc:.2f}')
        # ### add connect
        if self.is_covalent():
            allparts['CONNECT'].append(f'CONNECT {self.covalent_atomname}\n')
            allparts['ICOOR_INTERNAL'].append(self.get_icoor_from_ref())
        ordering = ['NAME', 'IO_STRING', 'TYPE', 'AA', 'ATOM', 'BOND_TYPE', 'CHI', 'CONNECT', 'NBR_ATOM',
                    'NBR_RADIUS',
                    'ICOOR_INTERNAL', 'PDB_ROTAMERS']
        with open(self.params_file, 'w') as w:
            for section in ordering:
                for line in allparts[section]:
                    w.write(line)

    def score_ligand_alone(self):
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([self.params_file])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, self.pdb_file)
        return pyrosetta.get_fa_scorefxn()(pose)

    def dummy(self):
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(self.params_file.replace('.params','.pdb'), 'ligand')
            pymol.cmd.fab('ACA', chain='A')
            pymol.cmd.alter('resn LIG', 'chain="B"')
            pymol.cmd.alter('resn LIG', 'resi="1"')
            pymol.cmd.sort()
            for atom in pymol.cmd.get_model('resn LIG').atom:
                pymol.cmd.translate([random.random() for i in range(3)], f'resn LIG and name {atom.name}')
            dummy = pymol.cmd.get_pdbstr('*')
        if self.is_covalent():
            dummy = f'LINK         SG  CYS A   2                {self.covalent_atomname} LIG B   1     1555   1555  1.8\n' + dummy
        pyrosetta.init(extra_options='-load_PDB_components false -no_optH true -relax:jump_move true')
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([f"{self.name}.params"])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, dummy)
        pose.dump_pdb(f'test_{self.name}.0.pdb')
        self.do_relax(pose)
        pose.dump_pdb(f'test_{self.name}.1.pdb')
        #pyrosetta.rosetta.protocols.docking.ConformerSwitchMover().apply(pose)
        docking = pyrosetta.rosetta.protocols.docking.DockMCMProtocol().apply(pose)
        pose.dump_pdb(f'test_{self.name}.2.pdb')


