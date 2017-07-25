%% downstream analysis of fitness data
%updated: July 2017, Lei Dai

%inputs: processed data of analysis_RF.m
clear;
clc;
inputfile='./data/downstream_5mismatch_v5.mat'; 
load(inputfile);

%outputs
%plot DFE of single aa mutants at 4 conditions
%save fitness values -> fit DFE by R codes
%-synonymous mutations
%-single aa mutations 
%-single nt mutations
%save workspace: rf_aa_v4.mat
%deposite fitness data in table format

%% plot DFE of single a.a. mutants
%group nt mutatnts into aa mutants
%transform rf_aa matrix to vector, filter lethal/missing mutants and WT entries
drug_concentration=[0;10;40;100];

aamut_total=86*19;

%fraction of lethal
for k=1:4
    rf_aa_vec{k}=[]; %WT is not included
    mutation_vec{k}={};
    fraction_lethal_condition(k)=nnz(find(rf_aa{k}==0))/aamut_total;
%     fraction_ben_condition(k)=nnz(find(rf_aa{k}>1))/aamut_total; %beneficial threshold:
end
fraction_lethal_condition

for i=1:length(WTsequence)
    for j=1:length(AAalphabet)
        for k=1:4
            if ~strcmp(mutation{i,j},'WT') && rf_aa{k}(i,j)>0 %not WT, not lethal
                mutation_vec{k}=[mutation_vec{k}; mutation{i,j}];               
                rf_aa_vec{k}=[rf_aa_vec{k};rf_aa{k}(i,j)];
            end
        end
    end
end

%% selection coefficient
%s=log(rf)
for k=1:4
    s_aa_vec{k}=log(rf_aa_vec{k});
end

%% plot histogram: DFE, [DCV]=10, 40, 100pM
for i=2:4
    subplot(3,1,i-1);
    interval=0;
    color='b';
    scale=0; %linear scale
    
    %non-synonymous mutation
    hist_custom(s_aa_vec{i},interval,color,scale);
    set(gca,'fontsize',10,'xlim',[-6 6],'ylim',[0 120]);
    xlabel('log(Relative Fitness)');
    ylabel('Count');
    title(strcat('[DCV]=',num2str(drug_concentration(i)),'pM'));

%draw lines: beneficial threshold
%synonymous mutations
s_silent=log(rf_silent(rf_silent>0));
threshold_ben=2*std(s_silent);
threshold_del=-2*std(s_silent);
% line([threshold_del threshold_del],[0 200],'color','k','linewidth',1);
line([threshold_ben threshold_ben],[0 200],'color','k','linewidth',2);

end

%% plot histogram: DFE, [DCV]=0pM + synonymous mutation
%single aa substitutions
interval=0;
color='b';
scale=0; %linear scale
subplot(2,1,1);
hist_custom(s_aa_vec{1},interval,color,scale);
set(gca,'fontsize',12,'xlim',[-6 6],'ylim',[0 120]);
xlabel('log(Relative Fitness)');
ylabel('Count');
title('single amino acid substitutions');
box off
hold on;

%synonymous
color='r';
interval=0.1;
scale=0;
subplot(2,1,2);
hist_custom(s_silent,interval,color,scale);
set(gca,'fontsize',12,'xlim',[-6 6]);
xlabel('log(Relative Fitness)');
ylabel('Count');
title('synonymous substitutions');
box off

%% plot fraction of beneficial mutations vs [DCV]
num_ben=[];
for k=1:4
    num_ben(k)=length(find(s_aa_vec{k}>threshold_ben));
end
plot(drug_concentration,num_ben,'-ko','markersize',12);
set(gca,'fontsize',12,'xlim',[-5 105],'ylim',[0 140]);
xlabel('[DCV] (pM)');
ylabel('Number of beneficial mutations');

%% save fitness data for fitting DFE of beneficial mutations
%1) fitness of single aa mutations
%transform rf_aa to vector: filter WT, keep lethal/missing
mutation_vec_all={};
for k=1:4
    rf_aa_vec_all{k}=[]; %no need to include WT
end
for i=1:length(WTsequence)
    for j=1:length(AAalphabet)
        if ~strcmp(mutation{i,j},'WT') %not WT
            mutation_vec_all=[mutation_vec_all; mutation{i,j}];
            for k=1:4
                rf_aa_vec_all{k}=[rf_aa_vec_all{k};rf_aa{k}(i,j)];
            end
        end
    end
end

%output: rf_aa_vec into table: aa change, rf at different drug concentration
mutation_table=cell2table(mutation_vec_all,'VariableNames',{'mutation'});
rf_table=table(rf_aa_vec_all{1},rf_aa_vec_all{2},rf_aa_vec_all{3},rf_aa_vec_all{4},'VariableNames',{'rf_0','rf_10','rf_40','rf_100'});
out_table = [mutation_table rf_table];

%save: all single aa mutations 
outfile='./data/fitness_singleaa.txt';
writetable(out_table,outfile,'Delimiter','\t');

%% save fitness of synonymous mutations
index_silent=find(NS5AFD.Synonymous == 1);
rf_silent=rf_filter(index_silent,1);
mutation_silent=NS5AFD.AAChange(index_silent);
out_table_silent=table(mutation_silent,rf_silent,'VariableNames',{'mutation','rf_silent'});

outfile='./data/fitness_synonymous.txt';
writetable(out_table_silent,outfile,'Delimiter','\t');

%% save fitness of nonsynonymous point mutations (single nt mutations)
%index_nt_nonsyn{1}: %from analysis_RF.m
mutation_pointmut=NS5AFD.Genotype(index_nt_nonsyn{1});
for i=1:4
    rf_pointmut{i}=rf_filter(index_nt_nonsyn{1},i);
end
out_table_pointmut=table(mutation_pointmut,rf_pointmut{1},rf_pointmut{2},rf_pointmut{3},rf_pointmut{4},'VariableNames',{'mutation','rf_0','rf_10','rf_40','rf_100'});

outfile='./data/fitness_pointmut.txt';
writetable(out_table_pointmut,outfile,'Delimiter','\t');

%% output: save rf_aa matrix, nonsynonymous substitutions grouped into single aa substituions
% save('./data/rf_aa_v4.mat','rf_aa','rf_aa_vec','rf_aa_vec_all','mutation','mutation_vec','mutation_vec_all');

%% deposit fitness as supplementary data
%format: table format
%1) single amino acid subsitutiton: mutation, relative fitness (0,10,40,100pM): fitness_singleaa.txt
%2) single codon subsitutiton: mutation, relative fitness (0,10,40,100pM): 
out_table_codon=table(NS5AFD.Genotype,NS5AFD.AAChange,NS5AFD.AApos,rf_filter(:,1),rf_filter(:,2),rf_filter(:,3),rf_filter(:,4),'VariableNames',{'Genotype','AAChange','AApos','rf_0','rf_10','rf_40','rf_100'});
out_table_codon = sortrows(out_table_codon,'AApos','ascend'); 

outfile='./data/fitness_singlecodon.txt';
writetable(out_table_codon,outfile,'Delimiter','\t');










