clear all
close all
global nPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color nfig shape;
global fileNames outputNames;

densTissues = [1, 2, 5, 6, 8, 9, 11, 12, 19, 20, 21];
nonDensTissues = [3, 4, 7, 10, 13, 14, 15, 16, 17, 18];
vascTissues = [4, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20];
nonVascTissues = [1, 2, 3, 5, 6, 9, 15, 17, 19, 21];
allTissues = [densTissues, nonDensTissues];

%nPar = 24;
%b = {'$T_{tum}$', '$T_{heal}$', '$T_{end}$', '$\bar{v}$', '$\alpha_{heal}$'...
%   '$\alpha_{tumG1}$', '$\alpha_{tumS}$', '$\alpha_{tumG2}$'...
%    '$\alpha_{tumM}$', '$\alpha_{tumG0}$', '$\alpha_{preEnd}$'...
%   '$\alpha_{neoEnd}$', '$pO_2^{nec}$', '$d$', '$D$', '$V_{max}$', '$K_M$'...
%   '$D^{VEGF}$', '$V_{max}^{', '$K_M^{VEGF}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$'...
%   '$pO_2^{hyp}$', '$v^{hyp}$'};

%nPar = 33;
%b = {'$T_{tum}$', '$T_{heal}$', '$T_{end}$', '$\bar{v}$', '$\alpha_{heal}$'...
%     '$\alpha/\beta_{heal}$', '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$'...
%     '$\alpha_{tumS}$', '$\alpha/\beta_{tumS}$' '$\alpha_{tumG2}$'...
%     '$\alpha/\beta_{tumG2}$', '$\alpha_{tumM}$', '$\alpha/\beta_{tumM}$'...
%     '$\alpha_{tumG0}$', '$\alpha/\beta_{tumG0}$', '$\alpha_{preEnd}$'...
%     '$\alpha/\beta_{preEnd}$', '$\alpha_{neoEnd}$', '$\alpha/\beta_{neoEnd}$'...
%     '$T_{arrest}$', '$pO_2^{nec}$', '$d$', '$D$', '$V_{max}$', '$K_M$'...
%     '$D^{VEGF}$', '$V_{max}^{VEGF}$', '$K_M^{VEGF}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$'...
%     '$pO_2^{hyp}$', '$v^{hyp}$'};

nPar = 34;
b = {'$T_{tum}$', '$T_{heal}$', '$T_{end}$', '$\bar{v}$', '$\alpha_{heal}$'...
    '$\alpha/\beta_{heal}$', '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$'...
    '$\alpha_{tumS}$', '$\alpha/\beta_{tumS}$' '$\alpha_{tumG2}$'...
    '$\alpha/\beta_{tumG2}$', '$\alpha_{tumM}$', '$\alpha/\beta_{tumM}$'...
    '$\alpha_{tumG0}$', '$\alpha/\beta_{tumG0}$', '$\alpha_{preEnd}$'...
    '$\alpha/\beta_{preEnd}$', '$\alpha_{neoEnd}$', '$\alpha/\beta_{neoEnd}$'...
    '$d_{thres}$', '$T_{arrest}$', '$pO_2^{nec}$', '$d$', '$D^{O_2}$'...
    '$V_{max}^{O_2}$', '$K_M^{O_2}$', '$D^{VEGF}$', '$V_{max}^{VEGF}$'...
    '$K_M^{VEGF}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$', '$pO_2^{hyp}$', '$v^{hyp}$'};
color = linspace(0, 1, nPar);
% 
%  
% nPar = 20;
% b = {'$T_{tum}$', '$T_{heal}$', '$T_{end}$', '$\bar{v}$', '$\alpha$'...
%     '$\alpha/\beta$', '$d_{thres}$', '$T_{arrest}$', '$pO_2^{nec}$'...
%     '$d$', '$D^{O_2}$', '$V_{max}^{O_2}$', '$K_M^{O_2}$', '$D^{VEGF}$'...
%     '$V_{max}^{VEGF}$', '$K_M^{VEGF}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$'...
%     '$pO_2^{hyp}$', '$v^{hyp}$'};
% color = linspace(0, 1, nPar);
shape = ['o', 's', 'v', 'd'];

nfig = 0;
quit = 0;

% path = uigetdir('../../Carlos/Results');
path = '../../Carlos/Results/Morris100_34Par_Cluster';
% path = '../../Carlos/Results/Morris100_OneAlphaBeta_Tissue4';
% path = '../../Carlos/Results/Morris100_8OutputsTest';

fileNames = {'EndTreatTumDens', '3MonTumDens', 'RecTumDens',...
    'FinTumVol', 'IntTumDens', 'TimeTo95', 'TimeTo99'...
    'RecTime'};
outputNames = {'Tumour density at the end of treat.', 'Tumour density 3 months after the end of treat.'...
    'Tumour density at recurrence', 'Final tumour volume', 'Integral of tumour density'...
    'Time to kill 95% of tumour cells', 'Time to kill 99% of tumour cells', 'Recurrence time'};
nOut = 8;

while(~quit)
    selOut = input(['Select an output [endTreaTumDens (1), 3MonTumDens (2), recTumDens (3), tumVol (4),\n'...
    'intTumDens (5), timeTo95 (6), timeTo99 (7), recTime (8) or all of them (-1)] or quit (0): ']);

    if(selOut >= 1 && selOut <= 8)
            tissue = input('Select one tissue (from 1 to 21) or all of them (0) or a mean over them (-1): ');
            if(tissue >= 1 && tissue <= 21)
                plotOutput(path, tissue, selOut)
            elseif(tissue == 0)
                for i = 1 : 21
                    plotOutput(path, i, selOut)
                end
            elseif(tissue == -1)
                plotMeanOutput(path, selOut)
            end
            
    elseif(selOut == -1)
            tissue = input('Select one tissue (from 1 to 21) or all of them (0) or a mean over them (-1): ');
            if(tissue >= 1 && tissue <= 21)
                plotAllOutputs(path, tissue)
            elseif(tissue == 0)
                for i = 1 : 21
                    plotAllOutputs(path, i)
                end
            elseif(tissue == -1)
                plotMeanAllOutputs(path)
            end
            
    elseif(selOut == 0)
            quit = 1;
    end
end




