%% Compare fitness data
%demonstrate reproducbility/consistency
%estimate experimental errors of RF estimates

%updated: July 2017
%compare fitness
%1) Hangfei's data: no drug
%2) round1/trasfection vs transfection/plasmid

%% input: fitness data
clear;
clc;
%fitness data, version 5.0: analysis_RF.m
inputfile='./data/downstream_5mismatch_v5.mat'; 
load(inputfile);

%%
%Hangfei's data: analysis_RF_HQdata.m
inputfile2='./data/HCV_HQdata_RF.mat'; 
load(inputfile2);

%% compare data sets: selection(one round)/transfection
%rank correlation between this data set and Hangfei's data set
cutoff_rf=0; %non-lethal mutations only

count=0;
RF_rep1=RF_control; %Hangfei Qi data
RF_rep2=rf_aa{1}; %this experiemnt
RF_rep1_vec=[];
RF_rep2_vec=[];

for i=1:length(WTsequence)
    for j=1:length(AAalphabet)
        if ~strcmp(mutation{i,j},'WT')  && ~isnan(RF_rep1(i,j)) && ~isnan(RF_rep2(i,j)) && RF_rep1(i,j)>cutoff_rf && RF_rep2(i,j)>cutoff_rf
            RF_rep1_vec=[RF_rep1_vec; RF_rep1(i,j)];
            RF_rep2_vec=[RF_rep2_vec; RF_rep2(i,j)]; 
            count=count+1;
            compare_aa(count,:)=[i,j];
        end
    end
end

%selection coefficient
s_rep1=log(RF_rep1_vec);
s_rep2=log(RF_rep2_vec);

%plot
plot(s_rep1,s_rep2,'+');
hold on;
x = -5:2; 
P = polyfit(s_rep1,s_rep2,1)
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');

% set(gca,'xlim',[-2 1],'ylim',[-2.5 1]);
% set(gca,'xscale','linear','yscale','linear');
corr_s=corr(s_rep1,s_rep2,'type','spearman')
legend(strcat('Spearman correlation=',num2str(corr_s)));
title('[DCV]=0 pM');
xlabel('log(Relative Fitness), Qi et al 2014');
ylabel('log(Relative Fitness)');

%%
% std: transform previous data to current data, then plot histogram
s_rep1_transform=s_rep1*P(1)+P(2);
s_dif_transform=s_rep2-s_rep1_transform;
s_dif_mean=mean(s_dif_transform)
s_dif_std=std(s_dif_transform)
histogram(s_dif_transform);
set(gca,'xlim',[-2.5 2.5]);
xlabel('\Delta s');
% legend(strcat('mean=',num2str(s_dif_mean),',std=',num2str(s_dif_std)));
legend(strcat('standard deviation=',num2str(s_dif_std)));

%% compare RF_no drug: transfection/plasmid vs. infection/transfection
%expect to be correlated
%rf_tran can be filter by a cutoff on wt_input
cutoff_rf=0;
ind_plot=find(rf_tran>cutoff_rf & rf_0>cutoff_rf);
s_tran=log(rf_tran(ind_plot));
s_sel=log(rf_0(ind_plot));
plot(s_tran,s_sel,'+');
corr_0=corr(s_tran,s_sel,'type','Spearman');
% set(gca,'xscale','log','yscale','log');
legend(strcat('Spearman correlation=',num2str(corr_0)));
xlabel('log(Relative Fitness),Transfection/Plasmid');
ylabel('log(Relative fitness)');
title('[DCV]=0 pM');

hold on;
x = -4:2; 
P = polyfit(s_tran,s_sel,1)
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');

