#!/usr/bin/python
"""
30 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from itertools               import groupby
from utils.seq_utils import translate
from sys                     import stderr


def parse_mcoffee_score(path, sequences):
    """
    """
    def is_header(line):
        return line=='\n'
    with open (path) as score:
        for is_header, lines in groupby (score, key=is_header):
            if is_header: continue
            for line in lines:
                if line.startswith ('T-COFFEE,'): break
                if line.startswith ('cons  '): continue
                name, sco = line.split ()
                sco = ['0' if n.isalpha() else n for n in sco]
                sequences[name].setdefault ('score', []).extend (sco)


def parse_mcoffee_aln(path, sequences):
    """
    """
    def is_header(line):
        return line=='\n'
    with open (path) as score:
        for is_header, lines in groupby (score, key=is_header):
            if is_header: continue
            for line in lines:
                if line.startswith ('CLUSTAL'): break
                if line.startswith ('    '): continue
                name, seq = line.split ()
                sequences[name].setdefault ('aa_ali', []).extend (seq)
