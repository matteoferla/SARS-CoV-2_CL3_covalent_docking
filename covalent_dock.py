########################################################################################################################

__doc__ = \
    """
    Automated covalent ligand docking via PyRosetta.
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2020 A.D."
__license__ = "MIT"
__version__ = "1"
__citation__ = ""

########################################################################################################################

import pyrosetta

pyrosetta.init(extra_options='-load_PDB_components false -no_optH true -relax:jump_move false -mute all')

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles
import pymol2


from locking_singleton_pymol import GlobalPyMOL
# pair_fit needs to be singleton.

import os

import json
import molfile_to_params
# see https://github.com/matteoferla/mol_to_params.py
# or
# https://www.well.ox.ac.uk/~matteo/molfile_to_params.zip


import re
from collections import defaultdict, namedtuple


class CovDock:
    apo_pdbfilename = 'apo.r.pdb'
    constraint_filename = 'cysBound.cst'

    def __init__(self, smiles: str,name: str='ligand', refine:bool=True):
        self.name = name
        self.smiles = smiles
        # cached property
        self._index_map = None
        self._name_map = None
        # operations
        self.ori_mol = Chem.MolFromSmiles(self.smiles)
        self.thio_mol = self.thiolate()
        # last one because of substitution
        self.CX_idx, self.SX_idx = self.thio_mol.GetSubstructMatches(Chem.MolFromSmiles('C[S-]'))[-1]
        self.dethio_mol = self.dethiolate()
        sdffile = f'{self.name}/{self.name}.sdf'
        self.save_confs(self.dethio_mol, sdffile)
        self.parameterise(sdffile)  # self.pdb_mol gets assigned.
        pdbblock = self.make_initial_pdb()
        self.pose = self.make_pose(pdbblock)
        self.dock_pose()
        if refine:
            self.refine_pose()
        self.pose.dump_pdb(f'{self.name}/docked_{self.name}.pdb')
        self.score = self.calculate_score()

    def thiolate(self) -> Chem.Mol:
        thio = AllChem.ReplaceSubstructs(self.ori_mol,
                                         Chem.MolFromSmiles('C[SiH3]'),
                                         Chem.MolFromSmiles('C[S-]'),
                                         replacementConnectionPoint=0)[0]
        thio.UpdatePropertyCache()
        thio = Chem.AddHs(thio)
        Chem.GetSSSR(thio)
        AllChem.EmbedMultipleConfs(thio, numConfs=100)
        AllChem.UFFOptimizeMoleculeConfs(thio, maxIters=2000)
        AllChem.ComputeGasteigerCharges(thio)
        return thio

    def dethiolate(self) -> Chem.Mol:
        dethio = Chem.EditableMol(self.thio_mol)
        dethio.RemoveAtom(self.SX_idx)
        return dethio.GetMol()

    def save_confs(self, out, file) -> True:
        # Chem.MolToMolFile(out,'test.mol')
        writer = rdmolfiles.SDWriter(file)
        writer.SetKekulize(False)
        for i in range(out.GetNumConformers()):
            writer.write(out, confId=i)
        writer.close()
        return True

    def thio2pdb_name(self, index: int) -> str:
        """
        Given an atom index of thio version. What is the atom name in params/PDB?

        :param index: thio version atom index
        :return: atom name in params
        """
        if self._index_map is None:
            self._fill_maps()
        return self._name_map[index]

    def thio2pdb_index(self, index: int) -> int:
        """
        Given an atom index of thio version. What is the atom index of pdb version

        :param index: atom index of thio version
        :return: atom index of pdb version
        """
        if self._index_map is None:
            self._fill_maps()
        return self._index_map[index]

    def _fill_maps(self):
        self._index_map = {k: v for v, k in enumerate(self.thio_mol.GetSubstructMatch(self.pdb_mol))}
        self._name_map = {k: self.pdb_mol.GetAtomWithIdx(v).GetPDBResidueInfo().GetName() for k, v in
                          self._index_map.items()}

    def get_icoor_from_ref(self) -> str:
        if self._index_map is None:
            self._fill_maps()
        # mol is thio_mol
        cx = self.thio_mol.GetAtomWithIdx(self.CX_idx)
        sx = self.thio_mol.GetAtomWithIdx(self.SX_idx)
        neigh = [a.GetIdx() for a in cx.GetNeighbors() if a.GetSymbol() not in ("H", "S")]
        grandneigh = [a.GetIdx() for a in self.thio_mol.GetAtomWithIdx(neigh[0]).GetNeighbors() if
                      a.GetSymbol() not in ("H",) and a.GetIdx() not in neigh + [self.CX_idx]]
        conf = self.thio_mol.GetConformer()
        angle = Chem.rdMolTransforms.GetAngleDeg(conf, self.CX_idx, self.SX_idx, neigh[0])
        distance = Chem.rdMolTransforms.GetBondLength(conf, self.CX_idx, self.SX_idx)
        torsion = Chem.rdMolTransforms.GetDihedralDeg(conf, self.CX_idx, self.SX_idx, neigh[0], grandneigh[0])
        return f'ICOOR_INTERNAL   CONN1  {torsion:.3f} {180 - angle:.3f} {distance:.3f} ' + \
               f'{self.thio2pdb_name(self.CX_idx)} {self.thio2pdb_name(neigh[0])} {self.thio2pdb_name(grandneigh[0])}\n'

    def fix_pdb(self, pdbfile):
        pdb_mol = Chem.MolFromPDBFile(pdbfile, removeHs=False)
        pdb_mol = AllChem.AssignBondOrdersFromTemplate(self.dethio_mol, pdb_mol)
        return pdb_mol

    def parameterise(self, infile: str) -> None:
        ## preventt append to preexisting
        if os.path.exists(f'{self.name}/{self.name}_conformers.pdb'):
            os.remove(f'{self.name}/{self.name}_conformers.pdb')
        ## make a params.
        molfile_to_params.run(infile,
                              conformers_in_one_file=True,
                              name='LIG',
                              amino_acid=None,
                              chain='B',
                              pdb=f'{self.name}/{self.name}',
                              clobber=True)
        ## modify
        self.pdb_mol = self.fix_pdb(f'{self.name}/{self.name}.pdb')
        ### split into parts
        allparts = defaultdict(list)
        with open(f'{self.name}/{self.name}.params') as r:
            for line in r:
                parts = line.split()
                allparts[parts[0]].append(line)
        ### change names
        name_of_CX = self.thio2pdb_name(self.CX_idx)
        new_name = ' CX '
        ### fix partial charges.
        for li in range(len(allparts['ATOM'])):
            line = allparts['ATOM'][li]
            rex = re.match('^ATOM\s+(\w+)(\s+\w+\s+\w+\s+)([.\d\-])*', line)
            name = rex.group(1)
            for ti, n in self._name_map.items():
                if n == name or n == name.strip() or n.strip() == name:
                    break
            else:
                raise ValueError(f'name {name} not in {list(self._name_map.values())}')
            gc = self.thio_mol.GetAtomWithIdx(ti).GetDoubleProp('_GasteigerCharge')
            allparts['ATOM'][li] = line.replace(rex.group(2) + rex.group(3), f'{gc:.2f}')
        ### add connect
        allparts['CONNECT'].append(f'CONNECT  CX\n')
        allparts['ICOOR_INTERNAL'].append(self.get_icoor_from_ref())
        ordering = ['NAME', 'IO_STRING', 'TYPE', 'AA', 'ATOM', 'BOND_TYPE', 'CHI', 'CONNECT', 'NBR_ATOM', 'NBR_RADIUS',
                    'ICOOR_INTERNAL', 'PDB_ROTAMERS']
        with open(f'{self.name}/{self.name}.params', 'w') as w:
            for section in ordering:
                for line in allparts[section]:
                    w.write(re.sub(name_of_CX + '(?!\w)', new_name, line))

        with open(f'{self.name}/{self.name}_conformers.pdb') as r:
            rotalib = r.read()
        with open(f'{self.name}/{self.name}_conformers.pdb', 'w') as w:
            w.write(re.sub(name_of_CX + '(?!\w)', new_name, rotalib))
        self.pdb_mol.GetAtomWithIdx(self.CX_idx).GetPDBResidueInfo().SetName(new_name)
        Chem.MolToPDBFile(self.pdb_mol, f'{self.name}/{self.name}.pdb')
        self._name_map[self.CX_idx] = new_name
        return None

    def make_initial_pdb(self):
        ## make a version with beta carbon on ligand.
        methio = AllChem.ReplaceSubstructs(self.ori_mol, Chem.MolFromSmiles('C[SiH3]'), Chem.MolFromSmiles('CSBr'),
                                           replacementConnectionPoint=0)[0]
        AllChem.EmbedMolecule(methio)
        cx, s, cb = methio.GetSubstructMatch(Chem.MolFromSmiles('CSBr'))
        methio.GetAtomWithIdx(s).SetMonomerInfo(
            Chem.rdchem.AtomPDBResidueInfo(atomName=' SX ', isHeteroAtom=True, residueName='UNL', residueNumber=1))
        methio.GetAtomWithIdx(cx).SetMonomerInfo(
            Chem.rdchem.AtomPDBResidueInfo(atomName=' CX ', isHeteroAtom=True, residueName='UNL', residueNumber=1))
        methio.GetAtomWithIdx(cb).SetMonomerInfo(
            Chem.rdchem.AtomPDBResidueInfo(atomName=' CB ', isHeteroAtom=True, residueName='UNL', residueNumber=1))
        ligpdb = Chem.MolToPDBBlock(methio)
        ## make a pdb
        # This step cannot be done with Pymol2 parallel.
        with GlobalPyMOL() as pymol:
            pymol.cmd.delete('*')
            pymol.cmd.load(self.apo_pdbfilename, 'apo')
            pymol.cmd.read_pdbstr(ligpdb, 'ligand')
            pymol.cmd.alter('resn UNL', 'chain="B"')
            pymol.cmd.alter('resn UNL', 'resn="LIG"')
            pymol.cmd.sort()
            pymol.cmd.pair_fit('resn LIG and name SX', 'resi 145 & name SG', 'resn LIG and name CB', 'resi 145 & name CB')
            pymol.cmd.remove('resn LIG and name SX')
            pymol.cmd.remove('resn LIG and name CB')
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return 'LINK         SG  CYS A 145                 CX  LIG B   1     1555   1555  1.8\n' + pdbblock

    def make_pose(self, pdbblock) -> pyrosetta.Pose:
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([f"{self.name}/{self.name}.params"])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdbblock)
        return pose

    def dock_pose(self):
        ## constraints
        def add_weights(scorefxn):
            stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
            scorefxn.set_weight(stm.score_type_from_name("atom_pair_constraint"), 20)
            scorefxn.set_weight(stm.score_type_from_name("angle_constraint"), 20)

        setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
        setup.constraint_file(self.constraint_filename)
        setup.apply(self.pose)
        scorefxn = pyrosetta.get_fa_scorefxn()
        add_weights(scorefxn)
        ### First relax
        movemap = pyrosetta.MoveMap()
        v = self.get_ligand_selector().apply(self.pose)
        n = self.get_neighbour_selector().apply(self.pose)
        movemap.set_bb(allow_bb=v)
        movemap.set_chi(allow_chi=n)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 1)
        relax.set_movemap(movemap)
        print(f'Settings set: {self.name}')
        relax.apply(self.pose)
        print(f'FastRelax: {self.name}')
        self.pose.dump_pdb(f'{self.name}/min1_{self.name}.pdb')
        with open(f'{self.name}/min1_{self.name}.json', 'w') as w:
            json.dump(self.calculate_score(), w)
        ### Second relax
        movemap.set_bb(allow_bb=n)
        relax.set_movemap(movemap)
        relax.apply(self.pose)
        print(f'FastRelax 2: {self.name}')
        self.pose.dump_pdb(f'{self.name}/min2_{self.name}.pdb')
        with open(f'{self.name}/min2_{self.name}.json', 'w') as w:
            json.dump(self.calculate_score(), w)
        ### Docking
        pyrosetta.rosetta.protocols.docking.setup_foldtree(self.pose, 'A_B', pyrosetta.Vector1([1]))
        scorefxn = pyrosetta.create_score_function('ligand')
        add_weights(scorefxn)
        restrict_to_focus = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),
                                                                                              self.get_neighbour_selector(), True)
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())
        tf.push_back(restrict_to_focus)
        docking = pyrosetta.rosetta.protocols.docking.DockMCMProtocol()
        docking.set_task_factory(tf)
        docking.set_ignore_default_task(True)
        #docking.set_move_map()
        docking.set_scorefxn(scorefxn)
        docking.apply(self.pose)
        print(f'Dock: {self.name}')
        with open(f'{self.name}/docked_{self.name}.json', 'w') as w:
            json.dump(self.calculate_score(), w)

    def get_ligand_selector(self):
        ligand_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
        ligand_selector.set_residue_name3('LIG')
        cys_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        cys_selector.set_index(self.pose.pdb_info().pdb2pose(chain='A', res=145))
        and_selector = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector()
        and_selector.add_residue_selector(cys_selector)
        and_selector.add_residue_selector(ligand_selector)
        return and_selector

    def get_neighbour_selector(self):
        and_selector = self.get_ligand_selector()
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        return NeighborhoodResidueSelector(and_selector, distance=5, include_focus_in_subset=True)


    def get_neighbour_vector(self):
        return self.get_neighbour_selector.apply(self.pose)

    def refine_pose(self):
        n_vector =  self.get_neighbour_vector()
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=n_vector)
        movemap.set_chi(allow_chi=n_vector)

        scorefxn = pyrosetta.get_fa_scorefxn()
        print(scorefxn)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
        relax.set_movemap(movemap)
        relax.apply(self.pose)

    def calculate_score(self):
        split_pose = pyrosetta.Pose()
        split_pose.assign(self.pose)
        cys_pos = split_pose.pdb_info().pdb2pose(chain='A', res=145)
        lig_pos = split_pose.pdb_info().pdb2pose(chain='B', res=1)
        # RESCON: 305 LIG n-conn= 1 n-poly= 0 n-nonpoly= 1 conn# 1 22 145 3
        split_pose.conformation().sever_chemical_bond(seqpos1=cys_pos, res1_resconn_index=3, seqpos2=lig_pos, res2_resconn_index=1)
        xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
        xyz.x = 500.0
        xyz.y = 0.0
        xyz.z = 0.0
        for a in range(1, split_pose.residue(lig_pos).natoms() + 1):
            split_pose.residue(lig_pos).set_xyz(a, split_pose.residue(lig_pos).xyz(a) + xyz)
        scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015")
        x = scorefxn(split_pose)
        b = scorefxn(self.pose)
        alone = self.score_ligand_alone()
        data = self.pose.energies().residue_total_energies_array()  # structured numpy array
        i = lig_pos - 1  ##pose numbering is fortran style. while python is C++
        ligand_data = {data.dtype.names[j]: data[i][j] for j in range(len(data.dtype))}
        return {'xyz_unbound': x,
                'bound': b,
                'apriori_ligand': alone,
                'ligand_data': ligand_data,
                'xyz_difference': b - x}


    @classmethod
    def create_apo(cls, template_pdbfilename: str):
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(template_pdbfilename)
            pymol.cmd.remove('solvent')
            pymol.cmd.remove('resn DMS')  # should run it through P Curran\s list of artefacts!
            apo = pymol.cmd.get_pdbstr('not resn LIG')
        pose = pyrosetta.Pose()
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, apo)
        scorefxn = pyrosetta.get_fa_scorefxn()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
        relax.apply(pose)
        pose.dump_pdb(cls.apo_pdbfilename)
        
    @classmethod
    def create_constraintfile(cls):
        with open(cls.constraint_filename, 'w') as w:
            w.write('AtomPair SG 145A CX 1B HARMONIC 1.8 0.2\n')
            w.write('Angle CB 145A SG 145A CX 1B HARMONIC 1.71 0.35\n')

    def score_ligand_alone(self):
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([f'{self.name}/{self.name}.params'])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, f'{self.name}/{self.name}.pdb')
        return pyrosetta.get_fa_scorefxn()(pose)
        

    
    


