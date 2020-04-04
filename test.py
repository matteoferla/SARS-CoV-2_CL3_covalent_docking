from covalent_dock import CovDock

def main():
    original_filename = 'Mpro-x0072.pdb'
    print(f'Create template file {ConDock.apo_pdbfilename} from {original_filename}')
    CovDock.create_apo(original_filename)
    print('dock ligand')
    CovDock(name='64_NCL', smiles='CC(N)NC(CC(N)=O)C(=O)NC1CCC(c2nc3cc(-c4ccc(S(N)(=O)=O)cc4)ccc3s2)CN1C(=O)C[SiH3]')
    print('Done')


if __name__ == '__main__':
    main()