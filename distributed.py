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

def get_master(master_file) -> list:
    master = {}
    for line in open(master_file):
        if line.strip == '':
            continue
        nid, smiles, oriname, hits = line.strip().split('\t')
        #1	Cc1ccc(OCC(=O)N2CCN(CCCCC(=O)Nc3cnccc3C)CC2)cc1	PET-SGC-fed-1	x0107
        master[oriname.replace(' ', '_')] = {'reacted_smiles': smiles,
                            'n_id': int(nid),
                            'original_name': oriname.replace(' ', '_'),
                            'hits': hits.split(',')}
    return master

def f(d): #'name', 'hits', 'smiles' keys.
    try:
        if os.path.exists(os.path.join('output', name, f'{name}.json')):
            return 'Already'
        else:
            from substitute import OverCov, Hit
            Hit.hits_path = '../Mpro'
            OverCov.hits_path = '../Mpro'
            #OverCov.placeholder = '[SiH3]'  # '*' is default
            OverCov(**d)
            return 'Done'

    except Exception as err:
        issue = f'{err.__class__.__name__}: {err}'.replace('\n', '')
        with open('error.txt', 'a') as w:
            w.write(f'{name}{issue}\n')
        print(issue)
        return err.__class__.__name__


if __name__ == '__main__':
    cores = int(os.environ['cores'])

    master = get_master('all_submissions.smi')
    databall = []
    for line in open('new_SiH3.smi'):
        if not line.strip():
            continue
        smiles, name = line.strip().split('\t')
        pre = re.match('^(.*)_\w+$', name).group(1)
        if pre not in master:
            pass
        else:
            hits = master[pre]['hits']
            databall.append({'name': name, 'hits': hits, 'smiles': smiles.replace('[SiH3]', '*')})
    with Pool(cores) as p:
        print(p.map(f, databall))
    slack_me('All DONE!')
