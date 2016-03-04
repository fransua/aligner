#!/usr/bin/python
"""
11 Jul 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from itertools       import groupby
from re              import sub
from re              import compile as compil
from sys             import stdout, stderr
from utils.seq_utils import translate, get_genetic_code
from subprocess      import Popen, PIPE
from time            import sleep

class Sequences():
    def __init__(self, genetic_code=None, path=None, aa=False, tmp_dir='tmp'):
        self.genetic_code = None if aa else get_genetic_code (genetic_code)
        self.__sequences = {}
        self.headers = []
        self.items = []
        self.tmp_dir = tmp_dir
        if path:
            self.read(path)
        self.alignments = {}
        
    def __getitem__(self, name):
        return self.__sequences[name]

    def __setitem__(self, name, value):
        self.__sequences[name] = value

    def __delitem__(self, name):
        del (self.__sequences[name])
        self.headers.pop(self.headers.index(name))

    def __iter__(self):
        for seq in self.__sequences:
            yield seq

    def __repr__(self):
        seqs = ''
        for head in self.headers:
            seqs += '>%s |%s\n' % (head, self[head]['descr'])
            for item in self.items:
                seqs += '    %s length: %s\n' % (item, len(self[head][item]))
        return seqs

    def __str__(self):
        seqs = ''
        for head in self.headers:
            seqs += '>%s |%s\n' % (head, self[head]['descr'])
            for item in self.items:
                seqs += '%s: %s\n' % (item, self[head][item])
        return seqs

    def align(self, tools=None, nprocs=1, refresh=0.1):
        i = 0
        alignments = []
        procs = {}
        try:
            while alignments or procs:
                if alignments and len (procs)<nprocs[0]:
                    line = alignments.pop(0)
                    i += 1
                    procs[i] = {'p': Popen(line, shell=True,
                                           stderr=PIPE, stdout=PIPE),
                                'cmd': line}
                for p in procs:
                    if procs[p]['p'].poll() is None:
                        continue
                    if procs[p]['p'].returncode == -9:
                        print ' WAHOOO!!! this was killed:'
                        print procs[p]
                        return
                    alignments[p] = {'cmd': procs[p]['cmd']}
                    alignments[p]['out'], alignments[p]['err'] = procs[p]['p'].communicate()
                    del (procs[p])
                    break
                sleep(refresh)
        except Exception as e:
            print 'ERROR at', i
            print e
        print '\n\nall jobs done...'

        
    def write(self, outfile=None, item='seq', reverse=False, width=60, descr=False):
        """
        Write sequence object to file in fasta format
        
        :argument None outfile: path to outfile, if None than, print to stdout
        :argument seq item: what to put in place of sequence
        :argument False reverse: wether to reverse or not the sequence
        :argument 60 width: number of sites per line when printing sequence
        :argument False descr: put description of sequence also, not recommended if you are not sure how the aligner will read it.
        
        """
        if outfile:
            out = open (outfile, 'w')
        else:
            out = stdout
        wsub = compil('([A-Za-z-]{'+str(width)+'})')
        for elt in self:
            if descr:
                out.write ('>%s |%s\n' % (elt, self[elt]['descr']))
            else:
                out.write ('>%s\n' % (elt))
            seq = self[elt][item][::-1] if reverse else self[elt][item]
            seq = seq if type(seq) is str else ''.join(seq)
            out.write ('%s\n' % (sub(wsub, '\\1\n', seq)))
        if outfile:
            out.close()
        
    
    def read(self, path, store='seq'):
        """
        """
        def check_header(line):
            return line.startswith('>')
        def get_sep (line):
            if ' |' in line:
                return ' |'
            elif '\t|' in line:
                return '\t|'
            elif '|' in line:
                return '|'
            elif '\t' in line:
                return '\t'
            elif ' ' in line:
                return ' '
            else:
                return '\n'
        sep = None
        with open(path) as fasta:
            for is_header, lines in groupby (fasta, key=check_header):
                text = ''.join (line.strip() for line in lines)
                if is_header:
                    if not sep:
                        sep = get_sep(text)
                    if sep == '\n':
                        head, descr = text[1:], ''
                    else:
                        head, descr = text[1:].split(sep, 1)
                    self [head] = {'descr': descr}
                else:
                    seq = text.replace('\n', '')
                    if self.genetic_code:
                        if len (text)%3 != 0:
                            seq = seq[:-(len(seq)%3)]
                        try:
                            self[head]['prot']  = translate(text, self.genetic_code,
                                                           stop=True)
                            self[head]['codon'] = [seq[i:i+3] for i in xrange(0,len(seq), 3)]
                        except KeyError: # in case sequence == 'Sequence unavailable'
                            del (self[head])
                            print >> stderr, 'No sequence found for ' + head
                            print >> stderr, text
                            continue
                    self [head][store] = seq
                    self.headers.append(head)
        for item in self[head]:
            if item == 'descr': continue
            self.items.append(item)
        
    