function [IC50mut]=fitIC50_simple(Wmut,c,IC50wt,hwt)
%updated: 06/15/2016, LD
%assume hmut=hwt
%if Wmut exceeds the maximum possible value, output IC50mut=Inf

%parameters
%Wmut: RFmut(drug=c)/RFmut(drug=0). from data
%c: drug concentration used for selection
%IC50wt,hwt: IC50 and slope of WT dose-response curve. known parameters

%assumptions
%1)only use validated dose-reponse curve of WT
%2)hmut=hwt. this can be relaxed if Wmut are measured at multiple drug
%concentrations
%3)duration of selection experiment and duration of the dose-response
%experiment is supposed to be the same. this can be relaxed
Wmut_max=1+(c/IC50wt)^hwt;
if Wmut<Wmut_max
    IC50mut=c/(((1+(c/IC50wt)^(hwt))/Wmut-1)^(1/hwt));
else
    IC50mut=Inf;
end

end