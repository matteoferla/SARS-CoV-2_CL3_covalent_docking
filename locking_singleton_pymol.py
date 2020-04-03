########################################################################################################################

__doc__ = \
    """
    Transpiler module for Michelanglo --- ``GlobalPyMOL`` class
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "3"
__citation__ = "Ferla et al. (2020) MichelaNGLo:  sculpting  protein  views on web  pages without coding. Bioinformatics"

########################################################################################################################

import os
from pprint import PrettyPrinter
pprint = PrettyPrinter().pprint
import pymol2

from threading import Lock

###############################################################
# Hardcoded variable.
# override with pymol.cmd.set('fetch_path', tempfolder)
tempfolder = os.getcwd() + '/michelanglo_app/temp'


###############################################################

class GlobalPyMOL:
    """
    SingletonPyMOL that works like the parallel but that waits for the other thread to release it.
    Use with a context manager!
    """
    pymol = pymol2.SingletonPyMOL()
    pymol.start()
    pymol.cmd.set('fetch_path', tempfolder)
    pylock = Lock()

    def __init__(self):
        pass

    def __enter__(self):
        if not self.pylock.acquire(timeout=60):
            # something hung up.
            self.pymol.cmd.remove('*')
            self.pymol.cmd.delete('*')
            self.pylock.release() #pointless roundtrip to be safe.
            self.pylock.acquire()
            return self.pymol
        else:
            self.pymol.cmd.delete('*')
            return self.pymol

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.pymol.cmd.delete('*')
        self.pylock.release()

    def kill(self):
        ## the assumption is that it died.
        self.pymol.cmd.reinitialize() #dangerous for other threads.
        if self.pylock.locked():
            self.pylock.release()
