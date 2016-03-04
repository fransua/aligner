#!/usr/bin/python
"""
30 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "1.0"

from subprocess import Popen, PIPE 
from os         import path, environ, access, X_OK, pathsep
from sys        import stderr
from sequences  import Sequences


class Aligner():
    def __init__(self, tool):
        self.tool = tool
        self.exe  = None
        self.rank = 0
        self.run  = None
        self.check_aligner()

    def check_aligner(self):
        """
        checks if aligners are installed.
        Forgets about the ones that are not installed.
        """
        def is_exe(fpath):
            return path.isfile(fpath) and access(fpath, X_OK)
        if self.tool in ALIGNERS:
            found = False
            ali_path, _ = path.split(ALIGNERS[self.tool]['bin'])
            if not ali_path:
                for loc in environ["PATH"].split(pathsep):
                    exe_file = path.join(loc, ALIGNERS[self.tool]['bin'])
                    if is_exe(exe_file):
                        self.exe = exe_file
                        found=True
                        break
            else:
                self.exe = ALIGNERS[self.tool]['bin']
                found=True
        if not found:
            raise('ERROR: %s executable not found in path' % self.tool)


def __run_muscle(path):
    return Popen([BINARIES['muscle']['bin'],
                  '-quiet', #'-stable',
                  '-maxiters' , '999',
                  '-maxhours' , '24 ',
                  '-maxtrees' , '100',
                  '-in'       , path,
                  '-out'      , path + '_muscle',
                  ], stdout=PIPE, stderr=PIPE)


def __run_mafft(path):
    return Popen ('%s --auto %s > %s_mafft' % (BINARIES['mafft']['bin'], path, path),
                  shell=True, stdout=PIPE, stderr=PIPE)


def __run_dialign(path):
    return Popen ('%s /usr/share/dialign-tx %s %s_dialign' \
                  % (BINARIES['dialign']['bin'], path, path),
                  shell=True, stdout=PIPE, stderr=PIPE)


def __run_probcons(path):
    return Popen ('%s %s > %s_probcons' \
                  % (BINARIES['probcons']['bin'], path, path),
                  shell=True, stdout=PIPE, stderr=PIPE)


ALIGNERS = { 'muscle'  : {'fun': __run_muscle  ,
                          'ord': 2             ,
                          'bin': 'muscle'      },
             'mafft'   : {'fun': __run_mafft   ,
                          'ord': 3             ,
                          'bin': 'mafft'       },
             'probcons': {'fun': __run_probcons,
                          'ord': 0             ,
                          'bin': 'probcons'    },
             'dialign' : {'fun': __run_dialign ,
                          'ord': 1             ,
                          'bin': 'dialign-tx'  }}
