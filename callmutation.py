"""
@author: Lei Dai
updated: 04/18/2017

summary: this script identifies mutations by comparing the mapped reads to NS5A reference sequence

input:
parse
./parse (output of parse.py)
mapping
./mapping_5mismatch (output of mapping.py)
reference
./reference/NS5A.fa

output:
./mutation_5mismatch/allmut
format: see function call mutation. region + condition(barcode) + mutation

"""

import os
import sys
import glob
import string
import time
from multiprocessing import Pool

#convert a sequence to reverse complementary strand
def rc(seq):
  complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB') 
  rcseq = seq.translate(complements)[::-1] #[::-1] string reverse. 5' to 3'
  return rcseq

def callmutation(R1parse,R2parse,mapfile,refseqs,outfile):
    R1file = open(R1parse,'r')
    R2file = open(R2parse,'r')
    mapfile = open(mapfile,'r')
    outfile = open(outfile,'w')
 
    for line in mapfile.xreadlines():           
        R_info=line.rstrip().rsplit("\t")    
        R1record = R1file.readline().rstrip().rsplit("\t") #check the output format of readrecord function        
        R2record = R2file.readline().rstrip().rsplit("\t") 
        
        #filter
        if R_info[0] == 'bad' : continue

        WT_Amp    = R_info[0] #RGN1-5, specify which region
        refseq    = refseqs[WT_Amp] 
        R1_strand = R_info[1] #strand
        R1_offset = int(R_info[2]) #offset
        R2_offset = int(R_info[3])    
   
        R1_bc    = R1record[1][0:3]
        R1_seq   = R1record[1][4::]
        R1_qual  = R1record[2].rsplit('-')[4::]  
        R2_seq   = R2record[1][4::]
        R2_qual  = R2record[2].rsplit('-')[4::]         
        #offset is determined based on 5'-3' direction (see offsetcheck funtion)
        R1_seq    = R1_seq[R1_offset:R1_offset+54]
        R2_seq    = R2_seq[R2_offset:R2_offset+54]
        R1_qual   = R1_qual[R1_offset:R1_offset+54]
        R2_qual   = R2_qual[R2_offset:R2_offset+54]    
        #reverse read
        if R1_strand == 'R': #R1 can be either forward or reverse reads
            R1_seq = rc(R1_seq) 
            R1_qual.reverse()
        else:
            R2_seq = rc(R2_seq)
            R2_qual.reverse()    
        #region1-4 only have 51bp
        if WT_Amp != 'RGN5': 
            region_length=51          
            R1_seq = R1_seq[0:region_length] #the entire 54 bp should be reversed before truncation
            R2_seq = R2_seq[0:region_length]
        else:
            region_length=54

        #CALL MUTATION
        Muts = [] #mutations 
        for n in range(0,region_length):
            if R1_seq[n] != refseq[n] and R1_seq[n] == R2_seq[n] and int(R1_qual[n]) >= 30 and int(R2_qual[n]) >= 30:  
                #R1==R2: mutation is called only if it is found on both of the paired-end reads    
                Mut = refseq[n]+str(n+1)+R1_seq[n] 
                Muts.append(Mut)              
        #set to WT, if no mutations are found        
        if len(Muts) == 0: 
            Muts = ['WT']
       
        outfile.write(WT_Amp+"\t" +R1_bc[0:3]+"\t" +'-'.join(Muts) +"\n") #region + condition(barcode) + mutation
    
    R1file.close()
    R2file.close()
    mapfile.close()
    outfile.close()

# process file i
def annotatefile(i):
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
    
    #specify input file
    #reads
#    filenames = sorted(glob.glob('/home/lei/Data/parse/*_R1')) 
    filenames = sorted(glob.glob('/mnt/hgfs/sequencing/Project_RS05262015/parse/*_R1'))     
    R1parse = filenames[i-1]
    fileID = R1parse.rsplit('/')[-1]
    fileID = fileID [:-3] #remove "_R1"
#    filenames = sorted(glob.glob('/home/lei/Data/parse/*_R2'))   
    filenames = sorted(glob.glob('/mnt/hgfs/sequencing/Project_RS05262015/parse/*_R2')) 
    R2parse = filenames[i-1]    
    
    #mapping info
    mapfile  = '/home/lei/Data/Project_RS05262015/mapping_5mismatch/'+fileID+'_map' 
    
    #output file
    folder ='/home/lei/Data/Project_RS05262015/mutation_5mismatch/'
    outfile  = folder+fileID+'_mut'        
    callmutation(R1parse,R2parse,mapfile,refseqs,outfile)

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
    poolnum = int(sys.argv[3])
    #create output folder
    folder = '/home/lei/Data/Project_RS05262015/mutation_5mismatch/'
    if not os.path.exists(folder):
       os.makedirs(folder) 
    
    p = Pool(poolnum)  # running separate processes
    p.map(annotatefile, range(file_start, file_end+1))  
    
    #combine all files
    os.system('cat '+folder+'*_mut > '+folder+'allmut')
    print "---total time: %s seconds ---" % (time.time() - start_time)

if __name__ == '__main__':
  main()