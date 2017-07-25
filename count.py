"""
@author: Lei Dai
updated: 04/18/2017

summary: this script counts the number of reads for each mutant in the library

input:
mutation
./mutation_5mismatch/allmut (output of callmutation.py)
reference
./reference/
NS5A.fa: reference genome
offset_bp: start bp position, NS5A region 1-5 (5 sub libraries)
offset_aa: start and end aa position, NS5A region 1-5 
BarCode: mapping between barcode and sample
NS5AMutTable(output of NS5AMuttable.py): enumerate all possible NNK mutants in the NS5A mutated region

output:
./count_5mismatch
depth: sequencing depth for each sub-library at each condition
AGenotypes: all mutants
SGenotypes: single mutants (nt level)

NS5A_FD: table of mutant counts -> downstream analysis
format (see header): aa position; genotype; WTaa-position#-mutantaa; # of bp changes;
mutant count; WT count for each region/amplicon; total read counts for each region/amplicon

"""

import os
import sys
import glob
import time

#this function can be saved as a separate file
def list2string(l):
  newl = []
  for i in l:
    newl.append(str(i))
  return newl

# change the amplicon position to whole NS5A mutregion position: bp
def adjustmutpos(muts, offset):
  if muts == 'WT': return muts
  else:
    muts  = muts.rsplit('-')
    mlist = []
    for m in muts:
      pos = str(int(m[1:-1])+offset-1) #m[0]:WT bp; m[-1]:mutant bp
      mlist.append(m[0]+pos+m[-1])
    return '-'.join(mlist)

##########MAIN###########
def main():
    start_time = time.time()
    
    ##specify input and output files
    #input folder
    infolder = '/home/lei/Data/Project_RS05262015/mutation_5mismatch/'
    #reference folder
    reffolder = '/mnt/hgfs/analysis_sunlab/analysis_HCV_NS5A/reference/'
    #output folder
#    outfolder='/home/lei/Data/Project_RS05262015/count_5mismatch/'
    outfolder='/mnt/hgfs/analysis_sunlab/analysis_HCV_NS5A/output/count_5mismatch/'
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    
    #READ OFFSET FILE: base pair position
    #output: offsetH
    filename = reffolder + 'offset_bp'
    infile  = open(filename,'r')
    offsetH = {}
    for line in infile.xreadlines():
      line = line.rstrip().rsplit("\t")
      offsetH[line[0]] = [int(line[1])] #line[0]: region; line[1]: offset
    infile.close()
    #READ IN OFFSET FILE: amino acid position
    #output: offsetA
    filename = reffolder + 'offset_aa'
    infile  = open(filename,'r')
    offsetA = {}
    for line in infile.xreadlines():
      line = line.rstrip().rsplit("\t")
      offsetA[line[0]] = [int(line[1]),int(line[2])] #line[0]: region; line[1]: start; line[2]: end
    infile.close()
    #READ BARCODE FILE
    #output: barcodes; pops
    filename = reffolder + 'BarCode'
    infile   = open(filename,'r') # Barcode in sequencing library preparation: demultiplex 
    barcodes = {}
    pops     = []
    for line in infile.xreadlines():
      line = line.rstrip().rsplit("\t")
      barcodes[line[0]] = line[1] #line[0]: 3bp barcode; line[1]: experimental condition
      pops.append(line[1])
    infile.close()
    
    #Initialize: GENOTYPE AND DEPTH 
    GenotypeH = {}
    depthH    = {}
    WTcount   = {}
    for pop in pops:
      GenotypeH[pop] = {} 
      depthH[pop]    = {}
      WTcount[pop]   = {}
      for amp in offsetH.keys():
        depthH[pop][amp] = 0
        WTcount[pop][amp] = 0
    
    ##Count the number of reads for each mutant
    #read: all mutations
    #file format: region; barcode (which condition); mutation(bp, position is relative to each region) 
    filename =infolder +'allmut'
    infile = open(filename,'r')    
    #output files
    filename = outfolder + 'depth'
    outfileD = open(filename,'w') 
    filename = outfolder + 'AGenotypes'#all mutants. not used in downstream analysis.
    outfileA = open(filename,'w') 
    filename = outfolder + 'SGenotypes' #single mutants (nt level). not used in downstream analysis.
    outfileS = open(filename,'w') 
    filename = outfolder + 'NS5A_FD'
    outfileC = open(filename,'w') #summary of counts. input for downstream analysis   
       
    for line in infile.xreadlines():
      line = line.rstrip().rsplit("\t")
      amp    = line[0]
      offset = offsetH[amp][0]
      pop    = barcodes[line[1]]
      muts   = adjustmutpos(line[2],offset)
      depthH[pop][amp] += 1
      if muts == 'WT':
        WTcount[pop][amp] += 1
      else:
        if GenotypeH[pop].has_key(muts): 
          GenotypeH[pop][muts] += 1
        else: 
          GenotypeH[pop][muts] = 1
    infile.close()
    
    ##Genotypes: list of all genotypes found in reads
    Genotypes = []
    for pop in pops:
      Genotypes.extend(GenotypeH[pop].keys()) 
    Genotypes = list(set(Genotypes))
    
    #WTCOUNT: total number of WT reads
    #DEPTH: total number of reads (succefully mapped)
    outfileD.write('pop'+"\t"+'amp'+"\t"+'WTcount'+"\t"+'Depth'+"\n") #header
    for pop in pops: #condition
      for amp in offsetH.keys(): #region
        outfileD.write("\t".join([pop, amp, str(WTcount[pop][amp]), str(depthH[pop][amp])])+"\n")
    outfileD.close()
    
    #GENOTYPES
    header   = 'Genotype'+"\t"+"\t".join(pops)
    outfileA.write(header+"\n")
    outfileS.write(header+"\n")
    for G in Genotypes:
      counts = []
      for pop in pops:
        if GenotypeH[pop].has_key(G): 
          counts.append(GenotypeH[pop][G])
        else: 
          counts.append(0)
      out = G+"\t"+"\t".join(list2string(counts))
      outfileA.write(out+"\n")
      if '-' not in G:
        outfileS.write(out+"\n")
    outfileA.close()
    outfileS.close()
    
    #IDENTIFY AA mutants
    #header
    header   = ['AApos','Genotype','AAChange','DNAdist','Synonymous']
    for pop in pops:
      header.extend([pop, pop+'_WT',pop+'_Dep']) 
    outfileC.write("\t".join(header)+"\n")
    
    #READ PROTEIN LEVEL DATA 
    filename = reffolder + 'NS5Amutlist'
    infile = open(filename,'r') 
    #this file is generated by NS5Amutlist.py
    #enumerate all possible NNK mutants in the NS5A mutated region
    #format: nt mutations; WT aa; aa position; mutant aa
    Rframe = {}
    for line in infile.xreadlines():
      array = line.rstrip().rsplit("\t")
      ID = array[0]
      WTaa1 = array[1]
      aapos = array[2]
      Mutaa1 = array[3]
      Rframe[ID] = ''.join([WTaa1,aapos,Mutaa1])
    infile.close()
    
    # all mutation, using '-' to seperate single & multiple bp mutation
    for G in Genotypes:
        if not Rframe.has_key(G):continue #exclude mutants not supposed to be in the library
        F1 = Rframe[G]
        aapos = F1[1:-1] 
        #count the # of bp changes
        dnadist=len(G.rsplit('-'))
        #check if the mutation is synonymous
        if F1[0] == F1[-1]: 
            synonymous = 1
        else: 
            synonymous = 0
        out = [str(aapos),G,F1,str(dnadist),str(synonymous)] 

        #output format (see header): aa position; genotype; WTaa-position#-mutantaa; # of bp changes;        
        #mutant count; WT count for each region/amplicon; total read counts for each region/amplicon
        for pop in pops:            
            for amp in offsetA.keys(): #find which region the mutation corresponds to
                if int(aapos) >= int(offsetA[amp][0]) and int(aapos) < int(offsetA[amp][1]):
                    wtc = WTcount[pop][amp]
                    dep = depthH[pop][amp]                     
            if GenotypeH[pop].has_key(G): 
                 out.append(str(GenotypeH[pop][G])) #mutant count 
            else: 
                 out.append(str(0)) #mutant count=0              
            out.append(str(wtc)) #WT count
            out.append(str(dep)) #depth
        outfileC.write("\t".join(out)+"\n")
    outfileC.close()
    
    print("---finished time: %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
  main()
