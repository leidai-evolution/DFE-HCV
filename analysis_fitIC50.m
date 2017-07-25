%% downstream analsis: drug resistance
%updated: July 2017
%fit dose response curves of validated mutants (raw data from Qi et al 2014)
%infer IC50 from drug resistnace score W: based on fitness profiles at DCV=0 and 10/40/100pM
%normalize W to 48 hr

clear;
clc;

%% find IC50 of validated mutants: fit dose response data 
load ./data/dosereponse %data reference: analysis_doseresponse.m

%include WT
mut_valid={'T24G','T24S','F28M','F28C','L31I','T54L','C92S','Y93F','Y93W','R56T','WT'};
%all mutants with validated IC50
% mut_valid={'T24G','T24S','F28M','F28C','L31I','T54L','C92S','Y93F','Y93W','R56T'};
%resistant mutants with validated IC50
% mut_valid={'T24G','T24S','F28M','F28C','L31I','T54L','C92S','Y93F','Y93W'};

for i=1:length(mut_valid)
    temp=eval(mut_valid{i});
    dose{i}=repmat(temp(1,:),1,3);
    response{i}=reshape(temp(2:4,:)',1,[]); %transposed because reshape's default is along columns
end

%fit to Hill function
hillcoef=[];
IC50=[];
figure(1);
for i=1:length(mut_valid)
    subplot(ceil(length(mut_valid)/4),4,i);
%     [hillcoef(i),IC50(i)]=doseResponse(dose{i},response{i});
    [hillcoef(i),IC50(i)]=doseResponse_fixh(dose{i},response{i});
    title(mut_valid{i});
end

%%  fitness of validated mutants
%reference: analysis_RF3.m
clc
out_table=readtable('./data/fitness_singleaa.txt','Delimiter','\t');
mut_valid_resistant={'T24G','T24S','F28M','F28C','T54L','C92S','Y93F','Y93W'}; %all resistant mutants except L31I (missing)
% mut_valid={'T24G','T24S','F28M','F28C','L31I','T54L','C92S','Y93F','Y93W'}; %all resistant mutants

rf_mutant=[];
for i=1:length(mut_valid_resistant)
    ind_mut=find(ismember(out_table.mutation, mut_valid_resistant{i})==1);
    rf_mutant(i,:)=table2array(out_table(ind_mut,2:5));
end

%% plot W vs [drug]
%W is expected to increase with [drug]
W_mutant=[];
for j=1:3
    W_mutant(:,j)=rf_mutant(:,j+1)./rf_mutant(:,1);
end
%add WT
figure(2);
drug_plot=[10,40,100];
for i=1:size(W_mutant,1)
    plot(drug_plot,W_mutant(i,:),'-o');
    hold on;
end
xlabel('[DCV]');
ylabel('W');

%% plot W vs measured IC50 for validated mutants
% index_temp=[1,2,3,4,6,7,8,9];
% plot(log10(W_mutant(:,2)),log10(IC50(index_temp)),'o');

%% infer IC50 using W 
%estimate T using W of 10 validated mutants
h=1;
drug=[10,40,100];
IC50_wt=IC50(end);
W_pd=[]; %predicted W based on measured IC50 (48hr)
for i=1:length(mut_valid) %exclude WT
    for j=1:length(drug)
        W_pd(i,j)=(1/(1+(drug(j)/IC50(i))^h))/(1/(1+(drug(j)/IC50_wt)^h));
    end
end

%
for i=1:size(W_pd,1)
    plot(W_pd(i,:),W_mutant(i,:),'-o');
    hold on;
end

%% calibration: transform Wobserved (144hr) to W (48hr) used in PD model 
W_mutant_all=log(W_mutant(:)); %measured, log transformed
W_pd_all=log(W_pd(:)); %predicted by PD model. log transformed

[corr_temp,p_temp]=corr(W_mutant_all,W_pd_all,'type','pearson')

plot(W_mutant_all,W_pd_all,'o','markersize',8);
hold on;
P = polyfit(W_mutant_all,W_pd_all,1)
x = 0:0.1:6; 
yfit = P(1)*x+P(2);
plot(x,yfit,'r-');
% set(gca,'xlim',[0 3],'ylim',[0 1],'fontsize',12);

% pars=fit(W_mutant_all,W_pd_all,'poly1');
% plot(pars);
polyfit_str = ['y = ' num2str(P(1)) ' *x + ' num2str(P(2))]
% polyfit_str qill be : y = 4*x + 2
xpos=1;
ypos=2.3;
text(xpos,ypos,polyfit_str);

xlabel('logW (144 hours post infection)');
ylabel('logW_{predict} (48 hours post infection)');

%% calculate W (fold change in RF) for all mutants
inputfile2='./data/rf_aa_v4.mat'; %from analysis_RF2.m
load(inputfile2);

drug=[10,40,100];
%rf_aa{2,3,4}->W10,40,100
for k=1:length(drug)
    for i=1:86
        for j=1:20
            if rf_aa{1}(i,j)>0 %filter missing/lethal variants
                Wobserved{k}(i,j)=rf_aa{k+1}(i,j)/rf_aa{1}(i,j);
            else
                Wobserved{k}(i,j)=NaN; 
            end
        end
    end
end

%% predict IC50 from W observed
drug=[10,40,100];
IC50_wt=IC50(end);
h_wt=1;

for i=1:86
    for j=1:20
        for k=1:length(drug)
            if ~isnan(Wobserved{k}(i,j))
                Wpd{k}(i,j)=(Wobserved{k}(i,j)^P(1))*10^(P(2));
                %                     Wpd{k}(i,j)=(Wmeasure{k}(i,j)^P(1));
                %                     Wpd{k}(i,j)=(Wmeasure{k}(i,j)^(1/3));
                %                     Wpd{k}(i,j)=P(1)*Wmeasure{k}(i,j);
                IC50mut{k}(i,j)=fitIC50_simple(Wpd{k}(i,j),drug(k),IC50_wt,h_wt);
            else
                IC50mut{k}(i,j)=NaN;
            end
        end
    end
end

%% histrogram of IC50 
for k=1:3
    figure(k);
    
    %exclude WT
    temp=[];
    for i=1:86
        for j=1:20
            if strcmp(mutation(i,j),'WT')==0
                temp=[temp IC50mut{k}(i,j)];
            end
        end
    end
    IC50mut_filter{k}=temp;
    length(IC50mut_filter{k})
    
    %fraction of mutations with IC50 above upper limit
    temp_num=sum(isinf(IC50mut_filter{k})) 
%     temp_frac=temp_num/(86*19)
    IC50mut_filter{k}(isinf(IC50mut_filter{k}))=10^5; %set mutants with IC50 above upper limit to 10^5
    
    %fraction of resistant mutations: >IC50_wt
    resistance_num=length(find(IC50mut_filter{k}>IC50_wt))
    %total: not NaN
    total_num=sum(~isnan(IC50mut_filter{k}))
    resistance_frac=resistance_num/total_num
    
    histogram(log10(IC50mut_filter{k}/IC50_wt),15);
    xlabel('log_{10}(fold change in IC_{50})');
    ylabel('Count');
    % hist_custom(temp,0,'b',1);
    set(gca,'xlim',[-2 5]);
    title(num2str(drug_plot(k)));
    legend(strcat('fraction above upper bound =',num2str(temp_num),'/1634'));
    
end

%% plot W/IC50 vs RF_nodrug
rf_plot=[];
W_plot=[];
 for i=1:86
        for j=1:20
            if strcmp(mutation(i,j),'WT')==0 && rf_aa{1}(i,j)>0 && Wobserved{1}(i,j)>0
                rf_plot=[rf_plot; rf_aa{1}(i,j)];
                W_plot=[W_plot; Wobserved{1}(i,j)];
            end
        end
 end

 plot(log(W_plot),log(rf_plot),'+');
[corr_temp,p_temp]=corr(log(W_plot),log(rf_plot),'type','pearson')
xlabel('log(Resistance score)');
ylabel('log(Relative Fitness), without drug');

line([0 0],[-6 2],'color','k');
line([-6 4],[0 0],'color','k');
set(gca,'fontsize',12);

%count fraction
num_total=length(rf_plot);
frac_fitres=length(find(rf_plot>1 & W_plot>1))/num_total
frac_delres=length(find(rf_plot<1 & W_plot>1))/num_total
frac_fitsen=length(find(rf_plot>1 & W_plot<1))/num_total
frac_delsen=length(find(rf_plot<1 & W_plot<1))/num_total

%% histogram: W
% histogram(log(W_plot))
% resistance_frac=length(find(W_plot>1))/num_total
% set(gca,'xlim',[-6 4]);
% xlabel('log(Resistance score)');
% ylabel('Count');
% set(gca,'fontsize',12);




