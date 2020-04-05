from typing import Union, List, Dict
from covalent_dock import CovDock, pymol2, GlobalPyMOL, pyrosetta
from hit import Hit

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdmolfiles, rdFMCS, Draw

import molfile_to_params
import os, re, shutil, random
from collections import defaultdict
from warnings import warn


class OverCov(CovDock):
    hits_path = '/Users/matteo/Coding/rosettaOps/Mpro'

    def __init__(self, smiles: str, hits: List, name: str = 'ligand', refine: bool = True):
        # entry attributes
        self.smiles = smiles
        self.name = name
        if not os.path.exists(self.name):
            os.mkdir(self.name)
        if len(hits) == 0:
            warn('using regular CovDock as without hits.')
            super().__init__(smiles, name, refine=True)
        else:
            self.hits = [Hit(hit) for hit in hits]
            # cached property
            self._index_map = None
            self._name_map = None
            self._best_hit = None
            # operations
            self.ori_mol = Chem.MolFromSmiles(self.smiles)
            self.thio_mol = self.thiolate()
            # last one because of substitution
            self.CX_idx, self.SX_idx = self.thio_mol.GetSubstructMatches(Chem.MolFromSmiles('C[S-]'))[-1]
            self.dethio_mol = self.dethiolate()
            ## align ligand
            cid = self.align_probe_to_target()
            aligned_file = self.write_probe(cid)
            self.save_confs(self.dethio_mol, aligned_file)
            self.parameterise(aligned_file)  # self.pdb_mol gets assigned.
            pdbblock = self.make_placed_pdb()
            open(f'{self.name}/pre_{self.name}.pdb','w').write(pdbblock)
            self.make_overlap_image()
            self.pose = self.make_pose(pdbblock)
            self.dock_pose()
            if refine:
                self.refine_pose()
            self.pose.dump_pdb(f'{self.name}/holo_{self.name}.pdb')
            self.snap_shot()
            self.score = self.calculate_score()

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

    def align_probe_to_targets(self) -> int:  # implace
        common = self._get_common(self.best_hit.mol)
        ### Align them
        overlap_target = self.best_hit.mol.GetSubstructMatch(common)
        for o, overlap_probe in enumerate(self.dethio_mol.GetSubstructMatches(common)):
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
        aligned_file = f'{self.name}/{self.name}.aligned.mol'
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
            open(f'{self.name}/{self.name}@{hit.name}.svg', 'w').write(svg)

    def make_placed_pdb(self):
        with GlobalPyMOL() as pymol:
            pymol.cmd.delete('*')
            self.best_hit.relax()
            pymol.cmd.load(self.best_hit.relaxbound_file, 'apo')
            # fix drift
            pymol.cmd.load(self.best_hit.bound_file, 'ref')
            pymol.cmd.align('apo', 'ref')
            pymol.cmd.delete('ref')
            pymol.cmd.remove('resn LIG')
            # distort positions
            pymol.cmd.load(f'{self.name}/{self.name}.pdb', 'ligand')
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
            pymol.cmd.load(f'{self.name}/{self.name}.pdb', 'pre')
            #pymol.cmd.load(f'pre_{self.name}.pdb', 'fudged')
            pymol.cmd.load(f'holo_{self.name}.pdb', 'fixed')
            for hit in self.hits:
                pymol.cmd.load(hit.bound_file)
            pymol.cmd.remove('solvent')
            pymol.cmd.show('sticks', 'resi 145')
            pymol.cmd.zoom('resi 145 or resn LIG')
            pymol.cmd.save(f'{self.name}/{self.name}.pse')

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
    h = Hit('x0692')
    print(h.score_ligand_alone())
    # # h.relax()
    # c = OverCov(name='JAN-GHE-fd8-1',
    #             smiles='CCNc1ncc(C#N)cc1CN1CCN(C(=O)C[SiH3])CC1',
    #             hits=['x0305', 'x0692', 'x1249'],
    #             refine=False)
    # print(c.score)
