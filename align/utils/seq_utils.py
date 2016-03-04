#!/usr/bin/python
"""
30 Sep 2011


"""

__author__  = "Francois-Jose Serra"
__email__   = "francois@barrabin.org"
__licence__ = "GPLv3"
__version__ = "0.0"

from itertools import product
from sys import stderr


def translate(sequence, gencode, stop=False):
    '''
    little function to translate DNA to protein...
    '''
    #dictionary with the genetic code
    proteinseq = ''
    #loop to read DNA sequence in codons, 3 nucleotides at a time
    l_seq = len (sequence)
    for n in xrange (0, l_seq-l_seq%3, 3):
        #checking to see if the dictionary has the key
        try:
            proteinseq += gencode[sequence[n:n+3]]
        except KeyError:
            newcod = [AMBIG[nt] for nt in sequence[n:n+3]]
            aa = ''
            for nt1, nt2, nt3 in product (*newcod):
                if nt1+nt2+nt3 in gencode:
                    aa  = gencode[nt1+nt2+nt3]
                    if aa != '*': break
            if not aa:
                aa = 'X'
            proteinseq += aa
    #return protein sequence
    if '*' in proteinseq[:-1]:
        stderr.write ('Found 1 STOP codon found inside.\n')
    if stop:
        if proteinseq.endswith('*'):
            return proteinseq[:-1]
        else:
            return proteinseq
    else:
        return proteinseq


def write_fasta(sequences, outfile, what='seq', lon=60):
    """
    """
    out = open (outfile, 'w')
    for elt in sequences:
        out.write ('>%s |%s\n' % (elt, sequences[elt]['descr']))
        seq = sequences[elt][what]
        seq = seq if type(seq) is str else ''.join(seq)
        out.write ('%s\n' % ('\n'.join ([seq[s:s+lon] for s in xrange (0, len (seq), lon)])))
    out.close()


def write_rfasta(sequences, outfile, what='seq', rev=False):
    """
    """
    out = open (outfile, 'w')
    for elt in sequences:
        out.write ('>%s\n%s\n' % (elt,
                                  sequences[elt][what][::(-1 if rev else 1)]))
    out.close()


def get_genetic_code (code):
    '''
    gencode, choose between:
        * [std] Standard
        * [vmt] Vertebrate Mitochondrial
        * [ymt] Yeast Mitochondrial
        * [mmt] Mold Mitochondrial, Protozoan Mitochondrial, Coelenterate Mitochondrial, Mycoplasma and Spiroplasma
        * [imt] Invertebrate Mitochondrial
        * [cnc] Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear
        * [emi] Echinoderm Mitochondrial and Flatworm Mitochondrial
        * [enu] Euplotid Nuclear
        * [bpp] Bacterial and Plant Plastid
        * [ayn] Alternative Yeast Nuclear
        * [ami] Ascidian Mitochondrial
        * [afm] Alternative Flatworm Mitochondrial
        * [bma] Blepharisma Macronuclear
        * [cmi] Chlorophycean Mitochondrial
        * [tmi] Trematode Mitochondrial
        * [som] Scenedesmus obliquus Mitochondrial
        * [thm] Thraustochytrium Mitochondrial
    '''
    gencode = {
        'std' :{
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T',
            'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L',
            'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',
            'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D',
            'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F',
            'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', '---':'-'},
        'vmt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D', 'GAC':'D',
            'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
            'TAA':'*', 'TAG':'*', 'AGA':'*', 'AGG':'*', '---':'-'},
        'ymt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'T', 'CTC':'T', 'CTA':'T', 'CTG':'T',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'mmt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'imt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'cnc' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAA':'Q', 'TAG':'Q',
            'TGT':'C', 'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L',
            'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H',
            'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R',
            'CGG':'R', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K',
            'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G',
            'GGC':'G', 'GGA':'G', 'GGG':'G', 'TGA':'*', '---':'-'},
        'emi' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'N', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'enu' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'bpp' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P',
            'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q',
            'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I',
            'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T',
            'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S',
            'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V',
            'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D',
            'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G',
            'GGG':'G', 'TAA':'*', 'TAG':'*', 'TGA':'*', '---':'-'},
        'ayn' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'S', 'CCT':'P',
            'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q',
            'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I',
            'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T',
            'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S',
            'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V',
            'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D',
            'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G',
            'GGG':'G', 'TAA':'*', 'TAG':'*', 'TGA':'*', '---':'-'},
        'ami' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'G', 'AGG':'G', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'afm' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAA':'Y', 'TGT':'C',
            'TGC':'C', 'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L',
            'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H',
            'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R',
            'CGG':'R', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'N',
            'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G',
            'GGC':'G', 'GGA':'G', 'GGG':'G', 'TAG':'*', '---':'-'},
        'bma' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAG':'Q', 'TGT':'C',
            'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'cmi' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAG':'L', 'TGT':'C',
            'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'tmi' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'N', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'som' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAG':'L', 'TGT':'C', 'TGC':'C',
            'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P',
            'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q',
            'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I',
            'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T',
            'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S',
            'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V',
            'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D',
            'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G',
            'GGG':'G', 'TCA':'*', 'TAA':'*', 'TGA':'*', '---':'-'},
        'thm' :{
            'TTT':'F', 'TTC':'F', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C', 'TGG':'W',
            'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P',
            'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
            'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I', 'ATC':'I',
            'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
            'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S',
            'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D', 'GAC':'D',
            'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
            'TTA':'*', 'TAA':'*', 'TAG':'*', 'TGA':'*', '---':'-'},
    }
    return gencode [code]
    
AMBIG = {'Y':['A', 'G'], 'R':['C', 'T'], 'M':['G', 'T'], 'K':['A', 'C'], \
         'S':['G', 'C'],'W':['A', 'T'], 'V':['C', 'G', 'T'], \
         'H':['A', 'G', 'T'], 'D':['A', 'C', 'T'], 'B':['A', 'C', 'G'], \
         'N':['A', 'C', 'G', 'T'], 'A': ['A'], 'C': ['C'], 'T': ['T'], 'G': ['G'] }

