"""
@author: Lei Dai
updated: 04/18/2017

summary: this script maps sequencing reads to HCV NS5A reference sequence
reads with > 5 mismatch to reference genome are filtered

input: 
./parse
sequencing reads (output of parse.py)
./reference
BarCode: barcode in sequencing library preparation
NS5A.fa: reference sequence

output: 
./mapping_5mismatch/*_map
output: see funtion offsetcheck

"""

import os
import sys
import glob
import string
import operator 
from itertools import imap
import time
import numpy as np
from multiprocessing import Pool

#convert a sequence to reverse complementary strand
def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') 
  rcseq = seq.translate(complements)[::-1] #[::-1] string reverse. 5' to 3'
  return rcseq

#calculate Hamming distance
def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(imap(operator.ne, str1, str2)) #interative map

def offsetcheck(R1_seq, R2_seq, refseqs): 
#INPUT 
#R_seq: read; refseqs: dictionary of references corresponding to RGN1-5
#OUTPUT
#region: RGN1-5; Roffset: starting point of the region; strand: read is F/R 
#"bad": if no match is found

  mismatch = 5; # number of mismatch allowed
  regionfound = 0
  regionmatch = 0  
  
  for i in range(0,len(R1_seq)-54+1):  
    #54bp=18aa. library1-4: 17aa, library 5: 18aa. 
    #3bp from primer (constant region) are added to library1-4 reference at the 3' end (i.e. bp #52-54)
    Fseq_1 = R1_seq[i:i+54] #forward sequence
    Rseq_1 = rc(Fseq_1) #reverse sequence
    for region in refseqs.keys(): #key: RNG1-5
      refseq  = refseqs[region]
      Fhdist_1  = hamming(refseq,Fseq_1)
      Rhdist_1  = hamming(refseq,Rseq_1)
      if Fhdist_1 <= mismatch:  
        strand_1 = 'F'
        offset_1 = i
        regionfound = 1
        break
      if Rhdist_1 <= mismatch: 
        strand_1 = 'R'
        offset_1 = i
        regionfound = 1
        break         
    if regionfound == 1: #region of R1 is found. map R2 on the same region  
        for j in range(0,len(R2_seq)-54+1):      
            Fseq_2 = R2_seq[j:j+54] 
            Rseq_2 = rc(Fseq_2)  
            if strand_1 == 'F': hdist_2 = hamming(refseq,Rseq_2)
            else: hdist_2 = hamming(refseq,Fseq_2)        
            if hdist_2 <= mismatch:
                offset_2 = j
                regionmatch = 1 #R1 and R2 reads mapped to the same region
                break 
    if regionfound == 1 and regionmatch == 1:        
      return [region, strand_1, offset_1, offset_2]     
      
  return 'bad'


# compare to the reference sequence: call mutation
def mapping(R1parse, R2parse, outfile, refseqs): 
#R1parse, R2parse: specify input file paths (FASTQ file, R1/R2 are paired-ends)
#mfile: m2file: specify output file paths
#refseqs: input 
  R1file = open(R1parse,'r')
  R2file = open(R2parse,'r')
  outfile  = open(outfile,'w')

  #READ BARCODE FILE
  #output: barcodes; pops
  reffolder = '/mnt/hgfs/analysis_sunlab/analysis_HCV_NS5A/reference/'
  reffilename = reffolder +'BarCode'
  reffile   = open(reffilename,'r') # Barcode in sequencing library preparation: demultiplex 
  barcodes = {}
  pops     = []
  for line in reffile.xreadlines():
      line = line.rstrip().rsplit("\t")
      barcodes[line[0]] = line[1] #line[0]: 3bp barcode; line[1]: experimental condition
      pops.append(line[1])
  reffile.close()      
  
  #read sequence from parsed data, map to reference sequence
  for line in R1file.xreadlines():
    R1record = line.rstrip().rsplit("\t") #check the output format of "readrecord" function
    R1_ID    = R1record[0] 
    R1_bc    = R1record[1][0:3] #3-bp barcode + 'T' (dA tailing)
    R1_seq   = R1record[1][4::]

    R2record = R2file.readline().rstrip().rsplit("\t") 
    R2_ID    = R2record[0]
    R2_bc    = R2record[1][0:3]
    R2_seq   = R2record[1][4::]
    
    #QUALITY CONTROL#
    #1) barcodes of paired-end reads should match
    #2) de-multiplex: discard reads of other samples (HCV NS5A sample ~40% of one HiSeq lane)    
    #3) discard if the read length is shorter than 54
    if R1_ID != R2_ID or R1_bc != R2_bc or (R1_bc not in barcodes.keys()) or len(R1_seq) < 54 or len(R2_seq) < 54:
        outfile.write('bad' +"\n")
        continue    
    
    #EXTRACT OFFSET INFO
    R_info   = offsetcheck(R1_seq,R2_seq,refseqs) 
    if R_info != 'bad':
        outfile.write(R_info[0]+"\t" +R_info[1]+"\t" +str(R_info[2])+"\t" + str(R_info[3]) +"\n")
    else:
        outfile.write('bad' +"\n")
        
# process file i
def processfile(i):
    print "--- processing file %d ---" %i
    start_time = time.time()
    
    #READ IN REFERENCE SEQUENCE
    reffolder = '/mnt/hgfs/analysis_sunlab/analysis_HCV_NS5A/reference/'
    reffilename = reffolder +'NS5A.fa'
    reffile = open(reffilename,'r') 
    refseqs = {}
    for line in reffile.xreadlines():
      if '#' in line: continue #skip comments
      if '>' in line: #label:RGN1-5
        ID = line.rstrip().replace('>','')
      else:
        refseqs[ID] = line.rstrip() #DNA sequence corresponding to RGN1-5    
    
    #specify input file: parsed reads
    filenames = sorted(glob.glob('/home/lei/Data/Project_RS05262015/parse/*_R1'))    
    R1parse = filenames[i-1]
    fileID = R1parse.rsplit('/')[-1]
    fileID = fileID [:-3] #remove "_R1"
    filenames = sorted(glob.glob('/home/lei/Data/Project_RS05262015/parse/*_R2'))    
    R2parse = filenames[i-1]
    #output file
    folder='/home/lei/Data/mapping_5mismatch/'
    outfile  = folder +fileID+'_map'      
    mapping(R1parse,R2parse,outfile,refseqs)

    print "---file %d finished time: %s seconds ---" %(i, time.time()-start_time) 

#######MAIN##############################     
#input: file_start, file_end, # of parallel process
def main():
    start_time = time.time()
    file_start  = int(sys.argv[1])
    file_end    = int(sys.argv[2])
    file_total = 43
    assert(file_start>=1)
    assert(file_end<=file_total)   

    #create output folder    
    folder='/home/lei/Data/Project_RS05262015/mapping_5mismatch/'
    if not os.path.exists(folder):
        os.makedirs(folder)

    poolnum = int(sys.argv[3])    
    p = Pool(poolnum)  # running separate processes
    p.map(processfile, range(file_start, file_end+1)) 
    print "---total time: %s seconds ---" % (time.time() - start_time)

if __name__ == '__main__':
  main()
