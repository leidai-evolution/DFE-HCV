%% plot ddg/fitness/entropy vs position
clear
clc

%% entropy: in vivo fitness
%input: pre-processed FASTA file in analysis_alignment_NS5A.m
inputfile1='./data_LANL/genotypeall_NS5AD1.mat';
load(inputfile1);
% entropy_plot=entropy(ind_start:ind_end);
entropy_plot=entropy;

%% import data from PyRosetta simulations: all mutations
% folder = 'C:\Users\LDai\Dropbox\Programming\pyrosetta\PyRosetta.monolith.ubuntu.release-103\ddg\';
folder='./data_structure/';
filename = strcat(folder,'3fqm_ddg.csv'); %repack 1 residue
% filename = strcat(folder,'3fqm_ddg_corrected.csv'); %repack 8A
ddg_all = readtable(filename,'Delimiter',',','ReadVariableNames',true);

%group ddg by position
%the list is ordered by position: 19*72=1368
n_res=103-32+1;
for i=1:n_res
    ind_start=(i-1)*19+1;
    ind_end=i*19;
    ddg_temp=ddg_all.ddg(ind_start:ind_end);
    ddg_mean(i)=mean(ddg_temp);
    ddg_median(i)=median(ddg_temp);
end

%% import fitness data, group by position
%copied from: analysis_alignment_NS5A_entropy.m
%fitness: RF is calculated in analysis_RF.m. %grouped into from analysis_RF2.m
% inputfile2='./downstream_june 2015/data/downstream2_5mismatch.mat'; 
inputfile='./data/downstream_5mismatch_v5.mat'; 
load(inputfile);
inputfile_fitness='./data/rf_aa_v4.mat';
load(inputfile_fitness);

%% fitness profile
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

%% RF vs site, entropy vs site
site=18:103;
[ax, h1, h2] = plotyy(site, fitness_median, site, entropy_plot, 'plot'); %median fitness
h1.Marker='o';
% h2.Marker = '+';

% Add title and x axis label
xlabel('Amino acid site');
% Use the axis handles to set the labels of the y axes
set(get(ax(1), 'Ylabel'), 'String', 'Relative fitness');
set(get(ax(2), 'Ylabel'), 'String', 'Shannon entropy');
set(gca,'ylim',[0 1.2],'fontsize',16);


%% output: classify positions by average fitness ->PyMol
%average fitness <0.2: red
%else: green
%output format: residue #, color
ind_start=32-18+1;
ind_end=103-18+1;
fitness_plot=fitness_mean(ind_start:ind_end);
for i=1:length(fitness_plot)
    if fitness_plot(i)<0.2
%         color_plot{i}='red';
        color_plot{i}='blue';
    else 
%         color_plot{i}='green';
        color_plot{i}='red';
    end
end
site=32:103;
output_table=table(site',color_plot');
% outfile='./data_structure/3fqm_color.txt';
outfile='./data_structure/3fqm_color_v2.txt';
writetable(output_table,outfile,'Delimiter','\t');

%% archive: stem plot of ddg/entropy/fitness vs site
ddg_plot=ddg_median;
% ind_start=32-18+1;
% ind_end=103-18+1;
% fitness_plot=fitness_lethal(ind_start:ind_end);
% fitness_plot=fitness_mean(ind_start:ind_end);
fitness_plot=fitness_mean;

site=32:103;
figure(2);
subplot(3,1,1);
%ddg
stem(site,ddg_median);
ylabel('\Delta\DeltaG');
set(gca,'xlim',[18 103]);
set(gca,'fontsize',16);
box off

site=18:103;
subplot(3,1,2);
%in vitro fitness
stem(site,fitness_plot);
ylabel('Relative fitness');
set(gca,'xlim',[18 103],'ylim',[-0.1 1.1]);
set(gca,'fontsize',16);
box off

site=18:103;
subplot(3,1,3);
stem(site,entropy_plot);
ylabel('Shannon entropy');
set(gca,'xlim',[18 103],'ylim',[-0.1 2]);
set(gca,'fontsize',16);
xlabel('Amino acid site');
box off