%% downstream analysis: spectrum of beneficial mutations under drug selection
%plot
%empirical CDF of beneficial single aa mutations at 4 conditions
%compare fitted CDF by Generalized Pareto Distribution and exponential distribution (output of R codes)
%empirical CDF of beneficial single nt mutations 

%load data
clear;
clc;
inputfile='./data/downstream_5mismatch_v5.mat';
load(inputfile);
inputfile2='./data/rf_aa_v4.mat'; %from analysis_RF2.m
load(inputfile2);

%% plot DFE of beneficial single aa mutations, without fitting
%plot CDF instead of histogram
%set threshold of beneficial mutations
%synonymous mutations
s_silent=log(rf_silent(rf_silent>0));
threshold_ben=2*std(s_silent);

for k=1:4
    rf_aa_ben{k}=[];
end
for k=1:4
    temp=log(rf_aa_vec{k}); %s=log(RF)
    s_aa_ben{k}=temp(temp>threshold_ben);
    cdfplot(s_aa_ben{k}); 
    hold on;
end
set(gca,'fontsize',12,'xlim',[0 6]);
xlabel('log(Relative Fitness)');
ylabel('Cumulative Distribution Function');
% legend(strcat('no drug, #=',num2str(length(rf_aa_ben{1}))),strcat('10pM [Daclatasvir],#=',num2str(length(rf_aa_ben{2}))),strcat('40pM [Daclatasvir],#=',num2str(length(rf_aa_ben{3}))),strcat('100pM [Daclatasvir],#=',num2str(length(rf_aa_ben{4}))));
drug_concentration=[0;10;40;100];
legend(num2str(drug_concentration));
% title('Spectrum of beneficial mutations');
box off;

%% transform to selection coefficients, normalized by the minimum beneficial mutation above threshold
for k=1:4
    s_shift{k}=s_aa_ben{k}-min(s_aa_ben{k}); %shift to minimum s
end

%% load: fitted parameters
parfile='./data/GPDfit_singleaa_log.txt'; %from MLEpar.R
MLE = readtable(parfile,'delimiter','\t');

%% plot CDF: both empirical and fitted
%single aa mutations
figure(1);
for k=1:4
    %CDF for data
    cdfplot(s_shift{k});
    hold on;
end
for k=1:4
        limit=max(s_shift{k});
    %full model
    tau_full=MLE.tau_full(k);
    kappa_full=MLE.kappa_full(k);
    %CDF for GPD
    [F, x]=GPDcdf(tau_full,kappa_full,limit);
    
    plot(x,F,'-k');
    hold on;
end
box off
set(gca,'fontsize',12);
xlabel('selection coefficient');
ylabel('Cumulative Distribution Function');
legend('0pM','10pM','40pM','100pM','Maximum Likelihood Fit');
% title('DFE of beneficial mutations');

%% compare fits of full model and reduced model
%single aa mutations
figure(2);
for k=1:4
    subplot(2,2,k);
    cdfplot(s_shift{k});
    hold on;
    hold on;
    limit=max(s_shift{k});
    %full model
    tau_full=MLE.tau_full(k);
    kappa_full=MLE.kappa_full(k);
    %CDF for GPD
    [F, x]=GPDcdf(tau_full,kappa_full,limit);
    %reduced model (exponential)
    tau_reduced=MLE.tau_reduced(k);
    kappa_reduced=MLE.kappa_reduced(k);
    [F_exp, x_exp]=GPDcdf(tau_reduced,kappa_reduced,limit);
    
    plot(x,F,'--k');
    hold on;
    plot(x_exp,F_exp,'-.');
    hold on;
    set(gca,'fontsize',12,'xlim',[0 6]);
    xlabel('selection coefficient');
    title(num2str(drug_concentration(k)));
end
legend('Empirical','Maximum likelihood fit: GPD','Maximum likelihood fit: exponential');

%% single nt mutations
%plot DFE of beneficial single nt mutations, without fitting
%plot CDF instead of histogram
for i=1:4
    temp=log(rf_filter(index_nt_nonsyn{1},i));
    s_pointmut_ben{i}=temp(temp>threshold_ben); 
    cdfplot(s_pointmut_ben{i}); 
    hold on;
end
set(gca,'fontsize',12,'xlim',[0 4]);
xlabel('log(relative fitness)');
ylabel('Cumulative Distribution Function');
drug_concentration=[0;10;40;100];
legend(num2str(drug_concentration));
% title('Spectrum of beneficial point mutations');
box off;


