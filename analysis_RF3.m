%% downstream analysis of fitness data: drug resistance
%updated: July 2017
%outputs
%fitness profile of validated resistance mutants: fitness vs [drug]

%% load data
% clear;
% clc;
out_table=readtable('./data/fitness_singleaa.txt','Delimiter','\t');

%% plot fitness change vs [drug]: validated resistant mutants
% mut_valid={'T24G','T24S','F28M','F28C','L31I','T54L','C92S','Y93F','Y93W'}; %validate resistant mutants
mut_valid={'T24G','T24S','F28M','F28C','L31I','T54L','C92S','Y93F','Y93W','R56T'}; %include the senstivie mutant R56T

%find rf of validated mutants
rf_mutant=[];
for i=1:length(mut_valid)
    ind_mut=find(ismember(out_table.mutation, mut_valid{i})==1);
    rf_mutant(i,:)=table2array(out_table(ind_mut,2:5));
end

drug_plot=[0,10,40,100];
c = distinguishable_colors(length(mut_valid));
for i=1:length(mut_valid)
    plot(drug_plot,log(rf_mutant(i,:)),'-o','linewidth',2,'markersize',8,'color',c(i,:));
    hold on;
end
set(gca,'fontsize',12);
xlabel('[DCV] (pM)');
ylabel('log(Relative Fitness)');
legend(mut_valid);
% ylim([-2 2.4]);
xlim([-5 105]);
line([0 100],[0 0],'color','k','linewidth',2,'linestyle','--');
box off;



