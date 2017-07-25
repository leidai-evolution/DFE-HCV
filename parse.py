"""
@author: Lei Dai
updated: 04/18/2017

summary: this script parses Illumina Paired-End data (FASTQ)

input: ./raw/*.fastq
Illumina HiSeq PE 100

output: ./parse/*_R1, ./parse/*_R2
format: see function readrecord

"""

import os
import sys
import glob
import string
from Bio import SeqIO
import time
from multiprocessing import Pool

#conver list to string
def list2string(l):
  newl = [] 
  for i in l:
    newl.append(str(i))
  return newl

#parse FASTQ file using SeqIO
def readrecord(filename,outfilename):
#OUTPUT format: ID \t sequence \t quality score
  outfile = open(outfilename,'w')
  for record in SeqIO.parse(filename, "fastq"):
    ID   = record.id
    seq  = record.seq
    qual = list2string(record.letter_annotations["phred_quality"])
    outfile.write(str(ID)+"\t"+str(seq)+"\t"+'-'.join(qual)+"\n")    
  outfile.close() 

#parse each file: raw sequencing data
def parsefile(i):    
    print("--- processing file %d ---" %i)  
    start_time = time.time()
    
    #specify file path
    filenames = sorted(glob.glob('/home/lei/Data/Project_RS05262015/raw/*_R1_*.fastq')) 
    folder='/home/lei/Data/'
    filename=filenames[i-1]    
    R1raw = filename
    R2raw = filename.replace('_R1','_R2') 
    fileID = filename.rsplit('/')[-1].rsplit('.')[0] 
    R1parse  = folder+'parse/'+fileID+'_R1' 
    R2parse  = folder+'parse/'+fileID+'_R2'
    readrecord(R1raw,R1parse)
    readrecord(R2raw,R2parse)
      
    print filename
    print("--- %s seconds ---" % (time.time() - start_time)) 

    
##########MAIN##########
#input: file_start, file_end
def main():
    file_start  = int(sys.argv[1])
    file_end    = int(sys.argv[2])
    file_total = 43
    assert(file_start>=1)
    assert(file_end<=file_total)   
    p = Pool(4)  # running separate processes
    p.map(parsefile, range(file_start, file_end+1))            
   
if __name__ == '__main__':
  main()
  