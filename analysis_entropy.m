%% compare entropy to fitness cost
%updated: 06/24/2015
clear;
clc

%input: pre-processed FASTA file in analysis_alignment_NS5A.m
inputfile1='./data_LANL/genotypeall_NS5AD1.mat';
load(inputfile1);
%fitness: RF is calculated in analysis_RF.m. %grouped into from analysis_RF2.m
inputfile2='./data/downstream_5mismatch_v5.mat'; 
load(inputfile2);
inputfile_fitness='./data/rf_aa_v4.mat';
load(inputfile_fitness);

%% fitness profile
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

%% plot: fitness vs entropy 
figure(1);
plot(entropy,fitness_mean,'ko','markerfacecolor','k','markersize',6);
xlabel('Shannon entropy');
ylabel('Relative fitness');
[pho pvalue]=corr(entropy,fitness_mean','type','spearman')
title(strcat('Spearman correlation=',num2str(round(pho*100)/100)));
set(gca,'xlim',[-0.05 2.1],'ylim',[-0.05 1.1]); %tune xlim
set(gca,'fontsize',16);

%plot regression line
P = polyfit(entropy,fitness_mean',1);
yfit = P(1)*entropy+P(2);
hold on;
plot(entropy,yfit,'r-');
box off;

%% plot entropy vs fraction of lethal
figure(2);
plot(entropy,fitness_lethal,'ko','markerfacecolor','k','markersize',6);
xlabel('Sequence diversity (entropy)');
ylabel('Proportion of lethal mutants');
[pho pvalue]=corr(entropy,fitness_lethal','type','spearman')
title(strcat('Spearman correlation=',num2str(round(pho*100)/100)));


