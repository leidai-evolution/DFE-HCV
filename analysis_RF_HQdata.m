%% HCV fitness data
clear
clc

%reference
%WTsequence: 86aa
%AAalphabet: a list of 20aa. the order is the same listed as Hangfei's raw data
load data/HCV_NS5A_reference

%raw data: /data_Hangfei/NS5A_AllroundsDrug-20121025_July2017_edit by LD
%paper: Qi et al, PLOS Pathogens 2014
%frequency relative to wt
%P1_control(transfection): freq_control_P1 
%P2_control(one round infection): freq_control_P2
%P1_drug (transfection): freq_BMS_P1
%P2_drug (one round infection)_freq_BMS_P2
load data/HCV_HQdata_freq

%% calculate RF
%selection: no drug
RF_control=freq_control_P2./freq_control_P1; %frequency is relative to wt
%selection: Daclatasvir 20pM 
RF_BMS=freq_BMS_P2./freq_BMS_P1;

%% find mutants absent in transfection pool: set fitness to NaN
index_control_filter=find(freq_control_P1==0);
RF_control(index_control_filter)=NaN;
index_BMS_filter=find(freq_BMS_P1==0);
RF_BMS(index_BMS_filter)=NaN;

%% save data: reference sequence, frequency, fitness
save data/HCV_HQdata_RF


