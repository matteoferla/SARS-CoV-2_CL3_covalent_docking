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


###########################

def g(h):
    try:
        from hit import Hit
        Hit(h).relax()
        return f'DONE {h}'
    except Exception as err:
        print(f'{err.__class__.__name__}: {str(err)}')
        return f'FAIL {h} {err.__class__.__name__}'

if __name__ == '__main__':
    cores = int(os.environ['cores'])
    path = '/well/brc/matteo/Mpro'
    hits = [folder.replace('_0', '').replace('Mpro-','') for folder in os.listdir(path) if os.path.isdir(os.path.join(path, folder))]
    print(hits)
    with Pool(cores) as p:
        print(p.map(g, hits))
    slack_me('All DONE!')