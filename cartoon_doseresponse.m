%% simulate dose response curves
%updated: 08/11/2016

%WT: high fitness w/o drug, low IC50 
%mutant: low fitness w/o drug, high IC50

%x-axis: drug concentration, log scale
%first y-axis: fitness 
%second y-axis: plot relative fitness, mark the line where RF=1
%function: plotyy; or yyaxis (2016a)

clear;
clc

% h=1; %Hill coefficient
f0=10;
num_interval=100;

IC50=1*ones(1,num_interval);
IC50_mut=10*ones(1,num_interval);

x=logspace(-2,3,num_interval); %drug
y=f0*IC50./(IC50+x); %fitness
y_mut=0.33*f0*IC50_mut./(IC50_mut+x);
R_mut=y_mut./y;

% %% subplots
% subplot(2,1,1);
% plot(x,y);
% hold on;
% plot(x,y_mut);
% set(gca,'xscale','log');
% 
% subplot(2,1,2);
% plot(x,R_mut);
% set(gca,'xscale','log','yscale','log');

%double y-axis
[hAx,hLine1,hLine2] =plotyy([x' x'],[y' y_mut'],x,R_mut,'semilogx','semilogx');
% set(gca,'xscale','log');
xlabel('Drug concentration') % normalized by IC_{50} of WT virus
ylabel(hAx(1),'Absolute fitness') % left y-axis
hAx(1).YTick=[0 2 4 6 8 10];
ylabel(hAx(2),'Relative fitness to WT virus') % right y-axis
hAx(2).YTick=[0 1 2 3 4];
hLine2.LineStyle = '--';
% hAx(2).YColor='k';
% hLine2.Color='k';
title('Dose response curves')
box off;
legend('WT virus','Drug-resistant virus');
% hold on;
% line([10^-2 10^3],[3 3],'color','k');