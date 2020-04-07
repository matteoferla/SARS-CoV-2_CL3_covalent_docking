from multiprocessing import Pool
import csv, re, unicodedata, requests, os
from warnings import warn


def slack_me(msg):
    """
    Send message to a slack webhook
    :param msg:
    :return:
    """
    # sanitise.
    msg = unicodedata.normalize('NFKD',msg).encode('ascii','ignore').decode('ascii')
    msg = re.sub('[^\w\s\-.,;?!@#()\[\]]','', msg)
    r = requests.post(url=os.environ['SLACK_WEBHOOK'],
                      headers={'Content-type': 'application/json'},
                      data=f"{{'text': '{msg}'}}")
    if r.status_code == 200 and r.content == b'ok':
        return True
    else:
        return False

def g(h):
    try:
        from hit import Hit
        Hit(h).relax()
    except Exception as err:
        print(f'{err.__class__.__name__}: {str(err)}')

def f(d):
    print('**********************************************')
    print('**********************************************')
    print('**********************************************')
    print('**********************************************')
    try:
        from substitute import OverCov
        print(dict(name=d['name'], hits=d['hits'], smiles=d['silane_smiles']))
        c = OverCov(name=d['name'], hits=d['hits'], smiles=d['silane_smiles'])
        s = c.score
        x = {k: str(v) for k,v in {**s, **d}.items()}
        print(x)
        with open(f'scores_{dataset}.csv','a') as w:
            dw =csv.DictWriter(w, fieldnames=['n_id',
                                              'original_name',
                                              'name',
                                              'hits',
                                              'silane_smiles',
                                              'reacted_smiles',
                                              'xyz_unbound',
                                              'bound',
                                            'apo',
                                            'apriori_ligand',
                                            'ligand_data',
                                            'xyz_difference',
                                            'apo_difference'])

            dw.writerow(x)
        #slack_me(d['name']+' - COMPLETE')
    except Exception as err:
        warn(f'FAILED! {err.__class__.__name}: {str(err)}')
        print(x)
        with open(f'scores_{dataset}.csv', 'a') as w:
            dw = csv.DictWriter(w, fieldnames=['n_id',
                                               'original_name',
                                               'name',
                                               'hits',
                                               'silane_smiles',
                                               'reacted_smiles',
                                               'xyz_unbound',
                                               'bound',
                                               'apo',
                                               'apriori_ligand',
                                               'ligand_data',
                                               'xyz_difference',
                                               'apo_difference'])

            dw.writerow(d)
        #slack_me(d['name']+' - FAILED')
        pass

def get_smiles(file):
    master = {}
    for line in open('all_submissions.smi'):
        nid, smiles, oriname, hits = line.strip().split('\t')
        #1	Cc1ccc(OCC(=O)N2CCN(CCCCC(=O)Nc3cnccc3C)CC2)cc1	PET-SGC-fed-1	x0107
        master[int(nid)] = {'reacted_smiles': smiles,
                            'n_id': int(nid),
                            'original_name': oriname.replace(' ', '_'),
                            'hits': hits.split(',')}
    data = []
    for line in open(file):
        #CCNc1ncc(CN)cc1CN1CCN(C(=O)C[SiH3])CC1	2_ACL
        smiles, name = line.split()
        #nid = int(name.split('_')[0]) Nir's ID is not longer used.
        data.append({**master[nid], 'silane_smiles': smiles, 'name': name.replace(' ', '_')})
    return data

if __name__ == '__main__' and 1==0:
    cores = int(os.environ['cores'])
    path = '/well/brc/matteo/Mpro'
    hits = [folder.replace('_0', '').replace('Mpro-','') for folder in os.listdir(path) if os.path.isdir(os.path.join(path, folder))]
    print(hits)
    with Pool(cores) as p:
        print(p.map(g, hits))


if __name__ == '__main__':
    cores = int(os.environ['cores'])
    dataset='new'
    with open(f'scores_{dataset}.csv', 'w') as w:
        dw = csv.DictWriter(w, fieldnames=['n_id',
                                              'original_name',
                                              'name',
                                              'hits',
                                              'silane_smiles',
                                              'reacted_smiles',
                                              'xyz_unbound',
                                              'bound',
                                            'apo',
                                            'apriori_ligand',
                                            'ligand_data',
                                            'xyz_difference',
                                            'apo_difference'])
        dw.writeheader()
    data = get_smiles(file=f'new_SiH3.smi')
    with Pool(cores) as p:
        print(p.map(f, data))
