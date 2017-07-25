%% plot RSA vs fitness
%updated, 08/02/2016

clear
clc
%RSA calculated by CalSASA.py
%import data: ./3fqm_rsa
% folder = 'C:\Users\LDai\Dropbox\Programming\structure\solvent accessibility\';
folder='./data_structure/';
filename = strcat(folder,'3fqm_rsa.csv');
ns5a_rsa= readtable(filename,'Delimiter',',','ReadVariableNames',true);

%% import fitness data, group by position
%consider merge this code w/ analysis_ddg_NS5A

%copied from: analysis_alignment_NS5A_entropy.m
%fitness: RF is calculated in analysis_RF.m. %grouped into single aa substitution from analysis_RF2.m
inputfile='./data/downstream_5mismatch_v5.mat';
load(inputfile);
inputfile_fitness='./data/rf_aa_v4.mat';
load(inputfile_fitness);

% fitness profile
fitness=rf_aa{1}; %in the absense of drug

for i=1:86
    WT_pos=find(~cellfun(@isempty,strfind(AAalphabet,WTsequence{i})));
    mut_pos=setdiff(1:20,WT_pos);
    fitness_pos=fitness(i,mut_pos);
    fitness_pos(isnan(fitness_pos))=[]; %remove NaN: missing variants
    %median,mean fitness; fraction of lethal mutants
    fitness_median(i)=median(fitness_pos);
    fitness_mean(i)=mean(fitness_pos);
    fitness_lethal(i)=nnz(find(fitness_pos==0))/length(fitness_pos);

end

%%
n_res=103-32+1;
ind_start=32-18+1;
ind_end=103-18+1;

fitness_plot=fitness_median(ind_start:ind_end); %use median
threshold_lethal=0.2;
ind_lethal=find(fitness_plot<threshold_lethal);
ind_nonlethal=find(fitness_plot>=threshold_lethal);

rsa_plot=ns5a_rsa.RSA(1:n_res);

%
nbin=5;
subplot(2,1,1);
histogram(rsa_plot(ind_lethal),nbin);
% interval=0.2;
% color='b';
% scale=0;
% hist_custom(rsa_plot(ind_lethal),interval,color,scale);
set(gca,'xlim',[0 1]);
ylabel('Count');
% xlabel('Relative solvent accessibility');
title('Median fitness<0.2');
set(gca,'fontsize',16);
box off;

subplot(2,1,2);
histogram(rsa_plot(ind_nonlethal),nbin);
set(gca,'xlim',[0 1]);
ylabel('Count');
xlabel('Relative solvent accessibility');
title('Median fitness>=0.2');
set(gca,'fontsize',16);
box off;

%% CDF
figure(2)
cdfplot(rsa_plot(ind_lethal));
hold on;
cdfplot(rsa_plot(ind_nonlethal));
legend('sites: median fitness<0.2','sites: median fitness>=0.2');
set(gca,'fontsize',12,'xlim',[0 1]);
xlabel('Relative Solvent Accessibility');

%% plot RSA vs fitness
figure(3)
plot(rsa_plot,fitness_mean(ind_start:ind_end),'ko','markerfacecolor','k','markersize',6);
[corr_temp p_temp]=corr(rsa_plot,fitness_mean(ind_start:ind_end)','type','spearman')
set(gca,'fontsize',12,'xlim',[-0.05 1.05],'ylim',[-0.05 1.25]);
xlabel('Relative Solvent Accessibility');
ylabel('Relative Fitness');

%plot regression line
P = polyfit(rsa_plot,fitness_mean(ind_start:ind_end)',1);
yfit = P(1)*rsa_plot+P(2);
hold on;
plot(rsa_plot,yfit,'r-');
box off;


