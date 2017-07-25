"""
@author: Lei Dai
updated: 04/19/2017

summary: this script generates a list of all possible mutants in NS5A aa 18-103 by NNK saturation mutagenesis
N: A/G/T/C, K: G/T
86*32=2752 nt genotypes
86*19=1634 single aa mutants

input: 
./reference/NS5A_all: nucleotide sequence of the entire NS5A mutated region

output:
./reference/NS5Amutlist
format: nt mutation; amino acid substitution

"""

import os
import sys
import math

#compare two codons to identify mutations
def makeID (codon,mutcodon,pos):
  ID = []
  for i in range(0,3):
    if codon[i]!= mutcodon[i]:
      ID.append(codon[i]+str(pos+1+i)+mutcodon[i])
  if len(ID) != 0:
    return ('-'.join(ID))
  else:
    return 'WT' 

#dictionary: DNA codon table
#'_': stop codon
dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','_'] #this variable is not used
bp = ['A','T','C','G']

#%% WT:DNA sequence of NS5A mutated region: amino acid 18-103
NS5A_bp=''
refWT = open('/mnt/hgfs/analysis_sunlab/analysis_HCV_NS5A/reference/NS5A_all','r') 
for line in refWT.xreadlines():
  if '>' in line: continue
  else:
    NS5A_bp= NS5A_bp+line.rstrip()
refWT.close()
print NS5A_bp

#%% WT:codon of mutated region
NS5A_codon={}
aapos=18 #starting amino acid position
for pos in range(0,len(NS5A_bp)):
  if pos%3 == 0:
    ID = aapos+ math.floor(pos/3)
    NS5A_codon[ID]= NS5A_bp[pos:pos+3] #NS5A_codon[0:18] is empty
print NS5A_codon

#%%
outfile = open('/mnt/hgfs/analysis_sunlab/analysis_HCV_NS5A/reference/NS5Amutlist','w')
for pos in range(0,len(NS5A_bp),3):
    AApos= aapos+math.floor(pos/3)
    WTaa= dnamap[NS5A_codon[AApos]]
    codon = list(NS5A_codon[AApos])
    for nt1 in bp: #N
        pos1 = nt1
        for nt2 in bp: #N
            pos2 = nt2
            for nt3 in ['G','T']: #K
                pos3 = nt3
                Mutcodon = [pos1,pos2,pos3] #list
                ID = makeID(codon,Mutcodon,pos)  
                Mutaa= dnamap[''.join(Mutcodon)] #transfrom from list to string
                out=[ID,WTaa,str(int(AApos)),Mutaa]
                outfile.write('\t'.join(out)+'\n')
outfile.close()
