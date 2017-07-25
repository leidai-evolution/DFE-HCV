%% Compare drug resistnace score W
%demonstrate reproducbility/consistency

%compare drug resistance score
%correlation between conditions: selection at 10, 40, 100pM
%compare to Hangfei's data: selection at 20pM

%save worksapce

%% input: fitness data
clear;
clc;
%fitness data
% inputfile='./data/downstream_5mismatch_v5.mat'; 
inputfile='./data/rf_aa_v4.mat'; %from analysis_RF2.m
load(inputfile);

%% calculate drug resistance score W
%group by aa
for i=1:3
    W_aa{i}=rf_aa{i+1}./rf_aa{1};
end

%not grouped by aa
% for i=1:3
%     W{i}=rf_filter(:,i+1)./rf_filter(:,1);
% end

%% compare W(fold change in RF) between different drug concentrations: RF10/RF0,RF40/RF0,RF100/RF0
%for convenience
W_10=W_aa{1};
W_40=W_aa{2};
W_100=W_aa{3};

cutoff_W=1e-3;
% cutoff_W=1;
upper_W=1e3; %range of plot

ind_plot=find(W_10>cutoff_W & W_10<Inf & W_40>cutoff_W & W_40<Inf & W_100>cutoff_W & W_100<Inf);

subplot(1,2,1);
plot(W_40(ind_plot),W_10(ind_plot),'+');
corr_W=corr(W_10(ind_plot),W_40(ind_plot),'type','Spearman');
title(strcat('Spearman correlation=',num2str(corr_W)));
% corr_W=corr(log(W_10(ind_plot)),log(W_40(ind_plot)),'type','Pearson');
% title(strcat('Pearson correlation=',num2str(corr_W)));
set(gca,'xscale','log','yscale','log');
set(gca,'xlim',[cutoff_W upper_W],'ylim',[cutoff_W upper_W]);
xlabel('W([DCV]=40pM)');
ylabel('W([DCV]=10pM)');

%compare W40 and W100
subplot(1,2,2);
plot(W_40(ind_plot),W_100(ind_plot),'+');
corr_W=corr(W_40(ind_plot),W_100(ind_plot),'type','Spearman');
title(strcat('Spearman correlation=',num2str(corr_W)));
% corr_W=corr(log(W_40(ind_plot)),log(W_100(ind_plot)),'type','Pearson');
% title(strcat('Pearson correlation=',num2str(corr_W)));
set(gca,'xscale','log','yscale','log');
set(gca,'xlim',[cutoff_W upper_W],'ylim',[cutoff_W upper_W]);
xlabel('W([DCV]=40pM)');
ylabel('W([DCV]=100pM)');

%compare W10 and W100
% subplot(1,3,2);
% plot(W_10(ind_plot),W_100(ind_plot),'+');
% corr_W=corr(W_10(ind_plot),W_100(ind_plot),'type','Spearman');
% title(strcat('Spearman correlation=',num2str(corr_W)));
% % corr_W=corr(log(W_10(ind_plot)),log(W_100(ind_plot)),'type','Pearson');
% % title(strcat('Pearson correlation=',num2str(corr_W)));
% set(gca,'xscale','log','yscale','log');
% set(gca,'xlim',[cutoff_W upper_W],'ylim',[cutoff_W upper_W]);
% xlabel('W_{10}');
% ylabel('W_{100}');

%% Hangfei's data: analysis_RF_HQdata.m
% clear
% clc
inputfile2='./data/HCV_HQdata_RF.mat'; 
load(inputfile2);

%calculate W: Hangfei Qi's data
W_HQ=RF_BMS./RF_control; %DCV=20pM

%% compare drug resistance 
cutoff_W=1e-3;

W_HQ_vec=[];
W_LD_vec=[];
W_LD=W_aa{1}; %DCV=10pM
for i=1:length(WTsequence)
    for j=1:length(AAalphabet)
        if ~strcmp(mutation{i,j},'WT') && W_HQ(i,j)>cutoff_W && W_LD(i,j)>cutoff_W  
            W_HQ_vec=[W_HQ_vec; W_HQ(i,j)];           
            W_LD_vec=[W_LD_vec; W_LD(i,j)];
        end
    end
end

plot(W_LD_vec,W_HQ_vec,'+');
set(gca,'xscale','log','yscale','log');
corr_W=corr(W_LD_vec,W_HQ_vec,'type','Spearman');
title(strcat('Spearman correlation=',num2str(corr_W)));
xlabel('W([DCV]=10pM)');
ylabel('W([DCV]=20pM),Qi et al 2014');

%% save workspace
save('./data/W_aa.mat');
