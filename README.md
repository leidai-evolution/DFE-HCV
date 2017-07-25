# DFE-HCV

Analysis: Illumina sequencing data
C:\Users\LDai\Dropbox\Programming\analysis_sunlab\analysis_HCV_NS5A
•	parse.py:  parse Illumina Paired-End data (FASTQ)
•	mapping.py: map sequencing reads to HCV NS5A reference sequence
•	callmuation.py: identify mutations by comparing the mapped reads to NS5A reference sequence
•	count.py: count the number of reads for each mutant in the library
•	NS5Amutlist.py:  generate a list of NNK mutants in the library
•	runall.sh: shell script. Pipeline
Note: check folder destinations before running the scripts.
Input: FASTQ files
Output (used in downstream analysis)
Folder: .\output\count_5mismatch
Filename: NS5A_FD 
Format: aa position; aa mutation; nt mutation; nt distance from WT; synonymous; read count (mutant; WT; depth) for all 6 conditions (input, transfection, selected at DCV 0, 10, 40 ,100pM)

Analysis: fitness
analysis_RF.m
-calculate mutant frequency
Filter low frequency mutants: missing mutants; lethal mutants
-calculate fitness
-group nt mutations with the same aa substitution
identify 1/2/3 nt mutations; synonymous mutations

analysis_RF_HQdata.m
-calculate fitness based on Hangfei Qi’s independent experiment

analysis_comparefitness.m
-compare fitness data sets
-estimate experimental errors of fitness

analysis_RF2.m
-plot DFE
-Output fitness data: fit DBFE
-Deposit fitness data in table format: supplementary data set (single codon substitutions; single amino acid substitutions)

analysis_fitdist.m
-plot
empirical CDF of beneficial single aa mutations at 4 conditions
compare fitted CDF by Generalized Pareto Distribution and exponential distribution (output of R codes)
empirical CDF of beneficial single nt mutations

analysis_geneticcode
-plot DFE of 1,2,3 nt mutations

----------------
Fit DFE: R codes adapted from Biesel et al 2007
MLEpar.R
-main function: handles input and output
-normalize fitness values of beneficial mutations by the threshold, which is determined by synonymous mutations
gpdmle.R
-perform Maximum Likelihood Estimate (MLE) of parameters of GPD and exponential distribution
-perform likelihood ratio test
-bootstrap of test statistics 
evalrtfunc_modified.R
-fixed two bugs from the original codes

----------
Analysis: resistance
analysis_RF3.m
-plot fitness profile of validated resistance mutants: fitness vs [drug]

analysis_compareresistance.m
-calculate W
-compare data sets

analysis_fitIC50.m
-fit dose response curves of validated mutants
-infer IC50 from fitness data
-plot resistance vs fitness

----------------------
Bioinformatics
analysis_alignment.m
-calculate Shannon entropy from MSA data

analysis_entropy.m
-plot entropy vs fitness

analysis_position_plot.m
-plot fitness/entropy vs amino acid site
-output: blue/red color for PyMOL visualization

analysis_RSA.m
-RSA calculated by CalSASA.py
-plot RSA vs fitness

analysis_ddg.m 
-stability ddG simulations: pyrosetta ./ddg folder, benchmark_ddg.m
-plot ddG vs fitness
