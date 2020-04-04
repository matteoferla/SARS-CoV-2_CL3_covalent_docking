# SARS-CoV-2_CL3_covalent_docking

Pyrosetta pipeline to do covalent docking of ligands given a SMILES string with a '[SiH3]' as the linking atom*.

	
## Make a template

	CovDock.create_apo('Mpro-x0705_0.pdb')

3 FastRelax cycles as Rosetta wants stuff energy minimised.

## Dock

	CovDock(name='64_NCL', smiles='CC(N)NC(CC(N)=O)C(=O)NC1CCC(c2nc3cc(-c4ccc(S(N)(=O)=O)cc4)ccc3s2)CN1C(=O)C[SiH3]')

Requires a constrain file `CovDock.constraint_filename`, which can be generated by `CovDock.create_constraintfile()`.

What it does is make a thiol version in rdkit, then a dethiol version which gets parameterised, however the output gets hacked: Gasteiger values are added and the connection add.
The atom name that connects to C145.SG is changed to CX.
The constraint is set to x20. The values were taken from the thioether of methionine.
Then it uses `DockMCMProtocol` mover to dock it better and finally a quick 3x fastrelax of the 6 Å neighbourhood —this can be disabled with `refine=False`.

It uses my modded `molfile_to_params` module ([runs as a module and py3 ported](https://github.com/matteoferla/mol_to_params.py)), which actually is the source of half my troubles due to the fact RDKit does not output Mol2 files.
The module is not in this folder as the original code has licencing issues.
The codebase for it is big, but that is because it does everything from scratch, so a less blackbox modules could be written building on RDKit...

One thing to note is that PyMOL singleton has to be used for `pair_fit`. The global singleton is from [another project](https://raw.githubusercontent.com/matteoferla/MichelaNGLo-transpiler/master/michelanglo_transpiler/locking_singleton_pymol.py) as it is basically a way to not worry about it being a singleton as it uses a lock.	

The `.pose` attribute is the PyRosetta pose of the docked ligand.

## Installation

Script written for Python 3.7. Requires the whole circus: `pyrosetta`, `pymol`, `rdkit`.

## Footnote

I was given the files like this, but I am at a total loss when it comes to dummy atoms (RDKit and Obabel croak on R or X, while * matches everything in RDKit) like "why is R not okay albeit not a leggit SMILES"? It really confuses me!
I use `Br` as my dummy atom, 
