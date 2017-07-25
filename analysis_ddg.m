%% downstream analysis: ddG vs fitness
%updated: July 2017
clear;
clc;

% import data from PyRosetta simulations: all mutations
% folder = 'C:\Users\LDai\Dropbox\Programming\pyrosetta\PyRosetta.monolith.ubuntu.release-103\ddg\';
folder='./data_structure/';
filename = strcat(folder,'3fqm_ddg.csv'); %repack 1 residue 
ddg_all_1res = readtable(filename,'Delimiter',',','ReadVariableNames',true);
filename = strcat(folder,'3fqm_ddg_corrected.csv'); %repack 8A 
ddg_all_8A = readtable(filename,'Delimiter',',','ReadVariableNames',true);
% filename = strcat(folder,'3fqm_ddg_background.csv'); %mutate background: this didn't work

%group ddg by position
%the list is ordered by position: 19*72=1368
n_res=103-32+1;
for i=1:n_res
    ind_start=(i-1)*19+1;
    ind_end=i*19;
    ddg_temp=ddg_all_1res.ddg(ind_start:ind_end);
    ddg_mean(i)=mean(ddg_temp);
    ddg_median(i)=median(ddg_temp);
    ddg_temp=ddg_all_8A.ddg(ind_start:ind_end);
    ddg_mean_8A(i)=mean(ddg_temp);
    ddg_median_8A(i)=median(ddg_temp);
end

%% import fitness data, group by position
%copied from: analysis_alignment_NS5A_entropy.m
%fitness: RF is calculated in analysis_RF.m. %grouped into from analysis_RF2.m
inputfile='./data/downstream_5mismatch_v5.mat'; 
load(inputfile);
inputfile_fitness='./data/rf_aa_v4.mat';
load(inputfile_fitness);

%fitness profile
fitness=rf_aa{1}; %in the absense of drug

for i=1:86
    WT_pos=find(~cellfun(@isempty,strfind(AAalphabet,WTsequence{i})));
    mut_pos=setdiff(1:20,WT_pos);
    fitness_pos=fitness(i,mut_pos);
    fitness_pos(isnan(fitness_pos))=[];
    %median,mean fitness; fraction of lethal mutants
    fitness_median(i)=median(fitness_pos);
    fitness_mean(i)=mean(fitness_pos);
    fitness_lethal(i)=nnz(find(fitness_pos==0))/length(fitness_pos);

end

%% plot: ddG vs fitness
%ddg_median vs. fitness_median

ddg_plot=ddg_median;
ind_start=32-18+1;
ind_end=103-18+1;
% fitness_plot=fitness_lethal(ind_start:ind_end);
fitness_plot=fitness_median(ind_start:ind_end);

plot(ddg_plot,fitness_plot,'ko','markerfacecolor','k','markersize',6);
xlabel('\Delta\DeltaG');
% ylabel('Proportion of lethal mutants');
ylabel('Relative fitness');
[pho pvalue]=corr(ddg_plot',fitness_plot','type','spearman')
title(strcat('Spearman correlation=',num2str(round(pho*100)/100)));
set(gca,'xlim',[-3 13],'ylim',[-0.05 1.1],'fontsize',16);
%plot regression line
P = polyfit(ddg_plot',fitness_plot',1);
yfit = P(1)*ddg_plot'+P(2);
hold on;
plot(ddg_plot',yfit,'r-');
box off;

%% compare different repacking
plot(ddg_median,ddg_median_8A,'+')
[corr_temp,p_temp]=corr(ddg_median',ddg_median_8A')

%% analysis: stability-function
threshold_destab=1;
threshold_lethal=0.2;
%find deleterious mutations that do not destabilize the protein
ind_lethal=find(fitness_plot<threshold_lethal);
ind_destab=find(ddg_plot<threshold_destab);
ind_functional=intersect(ind_lethal,ind_destab);
% threshold_conserved=0.01;
% ind_conserved=find(entropy_plot<threshold_conserved);
% ind_functional=intersect(ind_conserved,intersect(ind_lethal,ind_destab));

res_functional=ind_functional+32-1;
functional=[];
for i=1:length(res_functional)
    functional.res(i)=res_functional(i);
    functional.aa(i)=WTsequence{res_functional(i)-18+1};
end
functional



