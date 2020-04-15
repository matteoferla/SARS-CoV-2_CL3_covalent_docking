from typing import Union, List, Dict
from covalent_dock import CovDock, pymol2, GlobalPyMOL, pyrosetta
from hit import Hit

from collections import namedtuple

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdmolfiles, rdFMCS, Draw

import molfile_to_params
import os, re, shutil, random, json
from collections import defaultdict
from warnings import warn

from fragmenstein import Fragmenstein
from egor import Egor

class OverCov(CovDock):
    apo_pdbfilename = 'apo.r.pdb'
    constraint_filename = 'cysBound.cst'
    hits_path = '/Users/matteo/Coding/rosettaOps/Mpro'
    work_path = 'output'
    if not os.path.exists(work_path):
        os.mkdir(work_path)
    placeholder = '*' #'C[SiH3]'

    def __init__(self, smiles: str, hits: List, name: str = 'ligand', refine: bool = True):
        # entry attributes
        self.smiles = smiles
        self.name = name
        path = os.path.join(self.work_path, self.name)
        if not os.path.exists(path):
            os.mkdir(path)
        if len(hits) == 0:
            warn('using regular CovDock as without hits.')
            super().__init__(smiles, name, refine=True)
        else:
            print(f'Analysing: {name}')
            self.hits = [Hit(hit) for hit in hits]
            # cached property
            self._index_map = None
            self._name_map = None
            self._best_hit = None
            self.notebook = {}
            # operations
            self.best_hit.relax()
            self.ori_mol = Chem.MolFromSmiles(self.smiles)
            self.thio_mol = self.thiolate()
            # S has the highest index because of substitution, but prioritise the anteconnection atom that is a carbon
            x = []
            for element in ('C', 'N', 'O'):
                x.extend(self.thio_mol.GetSubstructMatches(Chem.MolFromSmiles(element + 'C[S-]')))
            self.CY_idx, self.CX_idx, self.SX_idx = sorted(x, key=lambda v: -v[2])[0]
            self.dethio_mol = self.dethiolate()
            ## align ligand
            cid = self.align_probe_to_target()
            aligned_file = self.write_probe(cid)
            self.save_confs(self.dethio_mol, aligned_file)
            print(f'SMILES converted: {name}')
            self.parameterise(aligned_file)  # self.pdb_mol gets assigned.
            print(f'Parameterised: {name}')
            #pdbblock = self.make_placed_pdb()
            self.fragmenstein = self.make_fragmenstein()
            pdbblock = self.place_fragmenstein()
            open(f'{self.work_path}/{self.name}/pre_{self.name}.pdb','w').write(pdbblock)
            self.make_overlap_image()
            self.pose = self.make_pose(pdbblock)
            print(f'PyRosetta loaded: {name}')
            self.egor = self.call_egor()
            print(f'EM: {name}')
            # self.dock_pose()
            # print(f'Docked: {name}')
            # # if refine:
            # #     self.refine_pose()
            # self.pose.dump_pdb(f'{self.work_path}/{self.name}/docked_{self.name}.pdb')
            print(f'Docked saved: {name}')
            self.snap_shot()
            print(f'Snapped: {name}')
            self.score = self.calculate_score()
            json.dump(self.notebook, open(f'{self.work_path}/{self.name}/{self.name}.json', 'w'))
            print(f'Done: {name}')

    @classmethod
    def from_fragmenstein_ouput(cls, name, path, best_hit):
        if not os.path.exists(path):
            raise FileNotFoundError('No folder')
        elif not os.path.exists(os.path.join(path, f'{name}.followup.mol')):
            raise FileNotFoundError('No placed molecule')
        elif not os.path.exists(os.path.join(path, f'{name}.thiol.mol')):
            raise FileNotFoundError('No thio molecule')
        else:
            followup = Chem.MolFromMolFile(os.path.join(path, f'{name}.followup.mol'))
            thio = Chem.MolFromMolFile(os.path.join(path, f'{name}.thiol.mol'))
            cls.from_mols(cls,
                          name=name,
                          followup=followup,
                          thio=thio)

    @classmethod
    def from_mols(cls, name: str, followup: Chem.Mol, thio: Chem.Mol, best_hit):
        self = cls.__new__(cls)
        self.smiles = Chem.MolToSmiles(followup, kekuleSmiles=False)
        self.name = name
        # cached property
        self._index_map = None
        self._name_map = None
        self._best_hit = best_hit
        self.notebook = {}
        # operations
        self.ori_mol = followup
        Chem.AddHs(self.ori_mol)
        self.thio_mol = thio
        Chem.AddHs(self.thio_mol)
        # Dethiolate manually. All this for a proton that will be out of place anyway.
        self.CX_idx, self.SX_idx = self.thio_mol.GetSubstructMatches(Chem.MolFromSmiles('C[S-]'))[-1]
        Chem.GetSSSR(thio)
        AllChem.EmbedMultipleConfs(thio, numConfs=100)
        AllChem.UFFOptimizeMoleculeConfs(thio, maxIters=2000)
        AllChem.ComputeGasteigerCharges(thio)
        AllChem.MMFFOptimizeMolecule(thio)
        dethio = Chem.EditableMol(thio)
        dethio.RemoveAtom(self.SX_idx)
        self.dethio_mol = dethio.GetMol()
        ## align ligand
        aligned_file = 'temp.mol'
        self.save_confs(self.dethio_mol, aligned_file)
        print(f'SMILES converted: {name}')
        self.parameterise(aligned_file)  # self.pdb_mol gets assigned.
        print(f'Parameterised: {name}')
        # pdbblock = self.make_placed_pdb()
        for ma, pa in zip(self.dethio_mol.GetAtoms(), self.pdb_mol.GetAtoms()):
            assert ma.GetSymbol() == pa.GetSymbol(), f'The indices do not align! {ma.GetIdx()}:{ma.GetSymbol()} vs. {pa.GetIdx()}:{pa.GetSymbol()}'
            ma.SetMonomerInfo(pa.GetPDBResidueInfo())
        with GlobalPyMOL() as pymol:
            pymol.cmd.delete('*')
            pymol.cmd.load(self.best_hit.relaxbound_file, 'apo')
            # fix drift
            pymol.cmd.load(self.best_hit.bound_file, 'ref')
            pymol.cmd.align('apo', 'ref')
            pymol.cmd.delete('ref')
            pymol.cmd.remove('resn LIG')
            # distort positions
            pymol.cmd.read_pdbstr(Chem.MolToPDBBlock(self.fragmenstein.positioned_mol), 'scaffold')
            pymol.cmd.save(f'{self.work_path}/{self.name}/{self.name}.scaffold.pdb')
            pymol.cmd.delete('scaffold')
            pymol.cmd.read_pdbstr(Chem.MolToPDBBlock(self.fragmenstein.positioned_mol), 'ligand')
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return 'LINK         SG  CYS A 145                 CX  LIG B   1     1555   1555  1.8\n' + pdbblock


        open(f'{self.work_path}/{self.name}/pre_{self.name}.pdb', 'w').write(pdbblock)
        self.make_overlap_image()
        self.pose = self.make_pose(pdbblock)
        print(f'PyRosetta loaded: {name}')
        self.egor = self.call_egor()
        print(f'EM: {name}')
        self.dock_pose()
        print(f'Docked: {name}')
        # if refine:
        #     self.refine_pose()
        self.snap_shot()
        print(f'Snapped: {name}')
        self.score = self.calculate_score()
        json.dump(self.notebook, open(f'{self.work_path}/{self.name}/{self.name}.json', 'w'))
        print(f'Done: {name}')





        @classmethod
        def reanimate(cls,
                      mol: Chem.Mol,
                      hits: List[Hit],
                      constraint_file: str,
                      ligand_residue: Union[str, int, Tuple[int, str], pyrosetta.Vector1],
                      key_residues: Union[None, Sequence[Union[int, str, Tuple[int, str]]], pyrosetta.Vector1] = None
                      ):
            fragmenstein = Fragmenstein(mol, hits)
            fragmenstein.positioned_mol

    def align_probe_to_target(self) -> int:  # implace
        """
        This aligns the probe to the first molecule in order to get most positions okay.

        :param probe: modified inplace
        :return: index of best conformer
        """
        ### find what is common
        common = self._get_common(self.best_hit.mol)
        ### Align them
        overlap_target = self.best_hit.mol.GetSubstructMatch(common)
        overlap_probe = self.dethio_mol.GetSubstructMatch(common)
        atomMap = [(probe_at, target_at) for probe_at, target_at in zip(overlap_probe, overlap_target)]
        rmss = [rdMolAlign.AlignMol(self.dethio_mol,
                                    self.best_hit.mol,
                                    prbCid=i,
                                    atomMap=atomMap,
                                    maxIters=500) for i in range(self.dethio_mol.GetNumConformers())]
        # print(rmss)
        best_i = rmss.index(min(rmss))
        return best_i

    def _get_common(self, target: Chem.Mol) -> Chem.Mol:
        res = rdFMCS.FindMCS([self.dethio_mol, target]  # ,
                                  # matchValences=True, threshold=0.1
                                  # atomCompare=Chem.rdFMCS.AtomCompare.CompareElements #,
                                  # bondCompare=Chem.rdFMCS.BondCompare.CompareAny #CompareOrder
                                  )
        return Chem.MolFromSmarts(res.smartsString)

    def write_probe(self, cid):
        aligned_file = f'{self.work_path}/{self.name}/{self.name}.aligned.mol'
        # Chem.MolToMolFile(probe, aligned_file, confId=cid)
        writer = rdmolfiles.SDWriter(aligned_file)
        writer.SetKekulize(False)
        writer.write(self.dethio_mol, confId=cid)
        for i in range(self.dethio_mol.GetNumConformers()):
            if i == cid:
                continue
            writer.write(self.dethio_mol, confId=i)
        writer.close()
        return aligned_file

    def get_hit2probe_map(self, hit_mol) -> Dict:
        overlap_hit, overlap_probe = self._get_overlaps(hit_mol)
        return {hit_at: probe_at for probe_at, hit_at in zip(overlap_probe, overlap_hit)}

    def get_probe2hit_map(self, hit: Hit) -> Dict:
        overlap_hit, overlap_probe = self._get_overlaps(hit)
        print(overlap_hit, overlap_probe)
        return {probe_at: hit_at for probe_at, hit_at in zip(overlap_probe, overlap_hit)}

    def _get_overlaps(self, hit: Hit):
        common = self._get_common(hit.mol)
        overlap_hit = hit.mol.GetSubstructMatch(common)
        overlap_probe = self.dethio_mol.GetSubstructMatch(common)
        return overlap_hit, overlap_probe

    def get_fudge_positions(self) -> Dict[int, List[int]]:
        """
        A fudged position is just the average of the atom position that overlap with two or more fragemnets.
        This has terrible and unexpected consequences.

        :return:
        """
        overlapx = defaultdict(list)
        overlapy = defaultdict(list)
        overlapz = defaultdict(list)
        for hit in self.hits:
            probe2hit_map = self.get_probe2hit_map(hit)
            hconf = hit.mol.GetConformers()[0]
            for pa, ha in probe2hit_map.items():
                overlapx[pa].append(hconf.GetAtomPosition(ha).x)
                overlapy[pa].append(hconf.GetAtomPosition(ha).y)
                overlapz[pa].append(hconf.GetAtomPosition(ha).z)
        mean = lambda a: sum(a)/len(a)
        print(overlapx)
        print(overlapy)
        print(overlapz)
        return {p: [mean(overlapx[p]), mean(overlapy[p]), mean(overlapz[p])] for p in overlapx.keys()}

    def make_overlap_image(self):
        # crappiest deepcopy
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(self.dethio_mol))
        AllChem.Compute2DCoords(mol)
        for hit in self.hits:
            overlap_hit, overlap_probe = self._get_overlaps(hit)
            drawer = Draw.rdMolDraw2D.MolDraw2DSVG(400, 300)
            Draw.rdDepictor.Compute2DCoords(mol)
            drawer.DrawMolecule(mol, highlightAtoms=overlap_probe)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            open(f'{self.work_path}/{self.name}/{self.name}@{hit.name}.svg', 'w').write(svg)



    def make_fragmenstein(self):
        ff = Fragmenstein(self.dethio_mol, [h.mol for h in self.hits])
        Chem.MolToMolFile(ff.scaffold, f'{self.work_path}/{self.name}/scaffold.fragmenstein.mol', kekulize=False)
        Chem.MolToMolFile(ff.chimera, f'{self.work_path}/{self.name}/chimera.fragmenstein.mol', kekulize=False)
        Chem.MolToMolFile(ff.positioned_mol, f'{self.work_path}/{self.name}/{self.name}.fragmenstein.mol', kekulize=False)
        return ff

    def place_fragmenstein(self):
        for ma, pa in zip(self.fragmenstein.positioned_mol.GetAtoms(), self.pdb_mol.GetAtoms()):
            assert ma.GetSymbol() == pa.GetSymbol(), f'The indices do not align! {ma.GetIdx()}:{ma.GetSymbol()} vs. {pa.GetIdx()}:{pa.GetSymbol()}'
            ma.SetMonomerInfo(pa.GetPDBResidueInfo())
        Chem.AddHs(self.fragmenstein.positioned_mol)
        with GlobalPyMOL() as pymol:
            pymol.cmd.delete('*')
            pymol.cmd.load(self.best_hit.relaxbound_file, 'apo')
            # fix drift
            pymol.cmd.load(self.best_hit.bound_file, 'ref')
            pymol.cmd.align('apo', 'ref')
            pymol.cmd.delete('ref')
            pymol.cmd.remove('resn LIG')
            # distort positions
            pymol.cmd.read_pdbstr(Chem.MolToPDBBlock(self.fragmenstein.positioned_mol), 'scaffold')
            pymol.cmd.save(f'{self.work_path}/{self.name}/{self.name}.scaffold.pdb')
            pymol.cmd.delete('scaffold')
            pymol.cmd.read_pdbstr(Chem.MolToPDBBlock(self.fragmenstein.positioned_mol), 'ligand')
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return 'LINK         SG  CYS A 145                 CX  LIG B   1     1555   1555  1.8\n' + pdbblock

    def call_egor(self):
        egor = Egor(pose=self.pose,
                    constraint_file=self.constraint_filename,
                    ligand_residue='LIG',
                    key_residues=['145A']
                    )
        self.notebook['pre-min'] = egor.ligand_score()
        mover = self.get_FastRelax(10)
        mover.apply(self.pose)
        self.notebook['post-min'] = egor.ligand_score()
        self.pose.dump_pdb(f'{self.work_path}/{self.name}/min1_{self.name}.pdb')
        return egor

    def make_placed_pdb(self):
        with GlobalPyMOL() as pymol:
            pymol.cmd.delete('*')
            # pymol.cmd.load(self.best_hit.relaxbound_file, 'apo')
            # # fix drift
            # pymol.cmd.load(self.best_hit.bound_file, 'ref')
            # pymol.cmd.align('apo', 'ref')
            # pymol.cmd.delete('ref')
            # pymol.cmd.remove('resn LIG')

            pymol.cmd.load(self.best_hit.apo_file, 'apo')
            # distort positions
            pymol.cmd.load(f'{self.work_path}/{self.name}/{self.name}.pdb', 'ligand')
            # fudge_positions = self.get_fudge_positions()
            # print(fudge_positions)
            # for i, atom in enumerate(pymol.cmd.get_model('resn LIG').atom):
            #     if atom.symbol == 'H':
            #         pymol.cmd.remove(f'resn LIG and name {atom.name}')
            #     elif i in fudge_positions:
            #         pymol.cmd.translate([fudge_positions[i][ax] - atom.coord[ax] for ax in range(3)], f'resn LIG and name {atom.name}', camera=0)
            #     else:
            #         pass #novel atom
            pdbblock = pymol.cmd.get_pdbstr('*')
            pymol.cmd.delete('*')
        return 'LINK         SG  CYS A 145                 CX  LIG B   1     1555   1555  1.8\n' + pdbblock

    # @property
    # def apo_pdbfilename(self):
    #     relaxed_filename = self.best_hit.relaxbound_file
    #     if not os.path.exists(relaxed_filename):
    #         self.best_hit.relax()
    #     original_filename = self.best_hit.apo_file
    #     relaxed_filename = os.path.join(self.hits_path, f'Mpro-{self.best_hit.name}_0', f'Mpro-{self.best_hit.name}_0_apo.r.pdb')
    #     if not os.path.exists(relaxed_filename):
    #         with pymol2.PyMOL() as pymol:
    #             pymol.cmd.load(original_filename)
    #             pymol.cmd.remove('solvent')
    #             pymol.cmd.remove('not polymer')  # should run it through P Curran\s list of artefacts!
    #             apo = pymol.cmd.get_pdbstr('*')
    #         pose = pyrosetta.Pose()
    #         pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, apo)
    #         scorefxn = pyrosetta.get_fa_scorefxn()
    #         relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
    #         relax.apply(pose)
    #         pose.dump_pdb(relaxed_filename)
    #     else:
    #         pass # already exists.
    #     return relaxed_filename

    def snap_shot(self):
        with pymol2.PyMOL() as pymol:
            pymol.cmd.bg_color('white')
            pymol.cmd.load(f'{self.work_path}/{self.name}/pre_{self.name}.pdb', 'pre')
            for k in ('pre','min'): #('min1','min2', 'docked'):
                pymol.cmd.load(f'{self.work_path}/{self.name}/{k}_{self.name}.pdb', k)
                pymol.cmd.align(k, 'pre')
            for hit in self.hits:
                pymol.cmd.load(hit.bound_file)
            # pymol.cmd.read_pdbstr(Chem.MolToPDBBlock(self.fragmenstein.positioned_mol), 'scaffold')
            pymol.cmd.remove('solvent')
            pymol.cmd.remove('resn DMS')
            pymol.cmd.show('sticks', 'resi 145')
            pymol.cmd.show('lines', 'byres resn LIG around 4')
            pymol.cmd.zoom('resi 145 or resn LIG')
            pymol.cmd.save(f'{self.work_path}/{self.name}/{self.name}_protein.pse')

    @property
    def best_hit(self) -> Hit:
        #cached.
        if self._best_hit is not None:
            return self._best_hit
        best_d = 99999
        best_hit = -1
        with pymol2.PyMOL() as pymol:
            for hit in self.hits:
                pymol.cmd.load(hit.bound_file)
                d = min([pymol.cmd.distance('resi 145 and name SG',f'resn LIG and name {atom.name}') for atom in pymol.cmd.get_model('resn LIG').atom])
                if d < best_d:
                    best_hit = hit
                    best_d = d
                pymol.cmd.delete('*')
        print('Best hit', best_hit.name)
        self._best_hit = best_hit
        return best_hit

if __name__ == '__main__':
    Hit.hits_path = '../Mpro'
    OverCov.hits_path = '../Mpro'
    # # h.relax()
    # c = OverCov(name='572_ACL',
    #         smiles='Cc1ncnc(-c2ccc(CN3CCN(C(=O)C[SiH3])CC3)cc2F)c1C',
    #         hits=['x0692', 'x0770', 'x0995'],
    #         refine=False)
    c = OverCov(name='2_ACL',
                        smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C[SiH3])CC1',
                        hits=('x0692', 'x0305', 'x1249'),
                        refine=False)
    print(c.score)
