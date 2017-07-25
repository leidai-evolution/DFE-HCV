%% analyze genetic code: HCV NS5A fitness data 
%updated: 07/12/2017
%outputs
%figure:
%-plot DFE of non-synonymous mutations (single codon substitutions, not grouped into aa mutants)
%-plot DFE of synonymous mutations (single codon substitutions, not grouped into aa mutants)
%-nonsense mutations, histogram
%-count fraction of lethal mutations for each category; count fraction of beneficial/neutral/deleterious

%figure
%violin plot of 1,2,3-nt non-synonymous mutations
%fraction of lethal/beneficial 1,2,3-nt non-synonymous mutations

% load data: output of analysis_RF.m
clear;
clc;
inputfile='./data/downstream_5mismatch_v5.mat'; 
load(inputfile);
inputfile2='./data/rf_aa_v4.mat';
load(inputfile2);

%% plot DFE: histogram, all non-synonymous mutations 
figure(1);
fitness_nt=[];

%exclude lethal mutations in DFE plot
% rf_0=rf_filter(:,1);
index_nonlethal=find(rf_0>0); 
index_plot=[];
for i=1:3 %non-synonymous mutations
    index_plot=[index_plot;intersect(index_nt_nonsyn{i},index_nonlethal)];
end
s_nonsyn=log(rf_0(index_plot));
%synonymous mutations
s_silent=log(rf_silent(rf_silent>0));

%non-synonymous
color='b';
interval=0.3;
scale=0;
subplot(2,1,1);
hist_custom(s_nonsyn,interval,color,scale);
hist_custom(rf_aa_vec{i},interval,color,scale);
set(gca,'fontsize',12,'xlim',[-5 2]);
% xlabel('log(Relative Fitness)');
ylabel('Count');
title('Non-synonymous');
box off
hold on;

%draw lines
% threshold_ben=2*std(s_silent);
% threshold_del=-2*std(s_silent);
% % line([threshold_del threshold_del],[0 200],'color','k','linewidth',1);
% line([threshold_ben threshold_ben],[0 200],'color','k','linewidth',2);

%synonymous
color='r';
interval=0.1;
scale=0;
subplot(2,1,2);
hist_custom(s_silent,interval,color,scale);
set(gca,'fontsize',12,'xlim',[-5 2]);
xlabel('log(Relative Fitness)');
ylabel('Count');
title('Synonymous');
box off
% line([threshold_del threshold_del],[0 40],'color','k','linewidth',1);
line([threshold_ben threshold_ben],[0 40],'color','k','linewidth',2);

%% DFE at the DNA level: 1,2,3-nt non-synonymous mutations, non-lethal 
%set 1) include synonymous and nonsense mutations: index_nt{i}
%set 2) exclude synonymous and nonsense mutations: index_nt_nonsyn{i}
figure(1);
interval=0;
color='b';
scale=1;
fitness_nt=[];

%exclude lethal mutations in DFE plot
rf_0=rf_filter(:,1);
index_nonlethal=find(rf_0>0); 
for i=1:3
    subplot(3,1,i);    
    %set 1) include synonymous and nonsense
%     index_plot=intersect(index_nt{i},index_nonlethal); 
    %set 2) non-synonymous only
    index_plot=intersect(index_nt_nonsyn{i},index_nonlethal);
    fitness_nt{i}=rf_0(index_plot); 
    hist_custom(fitness_nt{i},interval,color,scale);
    set(gca,'fontsize',12,'xlim',[-5 1]);
    xlabel('log(Relative Fitness)');
    ylabel('Count');
    title(strcat(num2str(i),'-nucleotide mutations, n=',num2str(length(fitness_nt{i}))));
end

%% fraction of lethal mutations 
index_lethal=find(rf_0==0); %find lethal mutations
for i=1:3
    %set 1) include 
%     index_lethal_nt{i}=intersect(index_nt{i},index_lethal); 
% frac_lethal(i)=length(index_nt_lethal{i})/length(index_nt{i});
    %set 2) exclude 
    index_nt_lethal{i}=intersect(index_nt_nonsyn{i},index_lethal);
    frac_lethal(i)=length(index_nt_lethal{i})/length(index_nt_nonsyn{i});
end
%plot
% subplot(2,1,1);
bar(frac_lethal,'b','EdgeColor','b');
set(gca,'ylim',[0 0.8]);
ylabel('Fraction of lethal mutations');
% xlabel('Number of nucleotide substitutions');
box off

%% fraction of beneficial mutations 
%test if there is enrichment of beneficial mutations in point mutations
index_beneficial=find(log(rf_0)>threshold_ben); %set to s>2*std(silent)

index_beneficial_nt=[];
for i=1:3
    %set 2)
    index_beneficial_nt{i}=intersect(index_nt_nonsyn{i},index_beneficial);
    length(index_beneficial_nt{i})
    length(index_nt_nonsyn{i})
    frac_beneficial(i)=length(index_beneficial_nt{i})/length(index_nt_nonsyn{i});
end

%all mutations
frac_beneficial_allnt=length(index_beneficial)/length(rf_0)

%plot
% subplot(2,1,2);
bar(frac_beneficial,'m','EdgeColor','m');
set(gca,'ylim',[0 0.05]);
ylabel('Fraction of beneficial mutations');
xlabel('Number of nucleotide substitutions');
box off


