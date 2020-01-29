clear all
close all
global nPar nTotPar selPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color colorbar nfig shape;
global fileNames outputNames;

% densTissues = [1, 2, 5, 6, 8, 9, 11, 12, 19, 20, 21];
% nonDensTissues = [3, 4, 7, 10, 13, 14, 15, 16, 17, 18];
% vascTissues = [4, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20];
% nonVascTissues = [1, 2, 3, 5, 6, 9, 15, 17, 19, 21];
% allTissues = [densTissues, nonDensTissues];

densTissues = [1, 2, 5, 6, 8, 9, 11, 19, 20];
nonDensTissues = [3, 4, 7, 10, 14, 15, 16, 17];
vascTissues = [4, 7, 8, 10, 11, 14, 16];
nonVascTissues = [1, 2, 3, 5, 6, 9, 15, 17, 19];
allTissues = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21];

green = [153, 255, 102];
darkGreen = [127, 207, 127];
lightGreen = [205, 255, 105];
blue = [127, 185, 225];
red = [255, 124, 128];
darkRed = [232, 136, 136];
lightRed = [255, 101, 152];
orange = [255, 174, 101];
gray = [173, 185, 202];
darkPurple = [160, 120, 196];
lightPurple = [146, 141, 251];

% path = '../../Carlos/Results/Sobol1000_38Par';
% b = {'$N$', '$tum$', '$T_{tum}$', '$hea$', '$T_{heal}$', '$ang$'...
%     '$T_{end}$', '$D^{VEGF}$', '$V_{max}^{VEGF}$', '$K_M^{VEGF}$'...
%     '$\bar{v}$', '$v^{hyp}$', '$\alpha_{heal}$', '$\alpha/\beta_{heal}$'...
%     '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$', '$\alpha_{tumS}$'...
%     '$\alpha/\beta_{tumS}$', '$\alpha_{tumG2}$'...
%     '$\alpha/\beta_{tumG2}$', '$\alpha_{tumM}$', '$\alpha/\beta_{tumM}$'...
%     '$\alpha_{tumG0}$', '$\alpha/\beta_{tumG0}$', '$\alpha_{preEnd}$'...
%     '$\alpha/\beta_{preEnd}$', '$\alpha_{neoEnd}$'...
%     '$\alpha/\beta_{neoEnd}$', '$d$', '$d_{thres}$', '$T_{arrest}$'...
%     '$oxy$', '$pO_2^{nec}$', '$D^{O_2}$', '$V_{max}^{O_2}$'...
%     '$K_M^{O_2}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$', '$pO_2^{hyp}$'};
% nTotPar = length(b);
% selPar = ones(1, nTotPar);
% selPar(1) = 0;
% selPar(4) = 0;
% selPar(6) = 0;
% selPar(32) = 0;

% path = '../../Carlos/Results/Sobol1000_38Par_2Gy';
% b = {'$N$', '$tum$', '$T_{tum}$', '$hea$'};
% nTotPar = length(b);
% selPar = ones(1, nTotPar);

path = '../../Carlos/Results/Sobol1000_12Par_2Gy';
b = {'$T_{tum}$', '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$'...
    '$\alpha_{tumS}$', '$\alpha_{tumG2}$', '$\alpha_{tumM}$'...
    '$T_{arrest}$', '$pO_2^{nec}$', '$D^{O_2}$', '$V_{max}^{O_2}$'...
    '$K_M^{O_2}$', '$pO_2^{preEnd}$'};
colorbar = [blue; lightPurple; lightPurple;
    lightPurple; lightPurple; lightPurple; darkPurple; darkGreen;
    lightGreen; lightGreen; lightGreen; lightGreen];
nTotPar = length(b);
selPar = ones(1, nTotPar);

% path = uigetdir('../../Carlos/Results');

b = b(selPar == 1);
nPar = nnz(selPar);
color = linspace(0, 1, nPar);
colorbar = colorbar(selPar == 1, :)./255;
shape = ['o', 's', 'v', 'd'];

nfig = 0;
quit = 0;

fileNames = {'EndTreatTumDens', '3MonTumDens'...
    'FinTumVol', 'IntTumDens', 'Killed50'...
    'Killed80', 'Killed90', 'Killed95', 'Killed99'...
    'Killed999', 'TimeTo95', 'TimeTo99', 'Rec'...
    'RecTumDens', 'RecTime'};
outputNames = {'Tumor density at the end of treat.'...
    'Tumor density at t = 12 weeks'...
    'Final tumor volume', 'Integral of tumor density'...
    '50% of tumor cells killed', '80% of tumor cells killed'...
    '90% of tumor cells killed', '95% of tumor cells killed'...
    '99% of tumor cells killed', '99.9% of tumor cells killed'...
    'Time to kill 95% of tumor cells'...
    'Time to kill 99% of tumor cells', 'Recurrence'...
    'Tumor density at recurrence', 'Recurrence time'};
nOut = 15;

while(~quit)
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\n intTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n '...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\n recTumDens (14), recTime (15), '...
        'or all of them (-1)]: ']);
    
    if(selOut >= 1 && selOut <= nOut)
        tissue = input(['Select one tissue (from 1 to 21) or all of '...
            'them (0) or a mean over them (-1): ']);
        if(tissue >= 1 && tissue <= 21)
            plotSobolOutput(path, tissue, selOut)
        elseif(tissue == 0)
            for i = 1 : 21
                plotSobolOutput(path, i, selOut)
            end
        elseif(tissue == -1)
            plotSobolMeanOutput(path, selOut)
        end
        
    elseif(selOut == -1)
        tissue = input(['Select one tissue (from 1 to 21) or all of '...
            'them (0) or a mean over them (-1): ']);
        if(tissue >= 1 && tissue <= 21)
%             plotSobolAllOutputs(path, tissue)
        elseif(tissue == 0)
            for i = 1 : 21
%                 plotSobolAllOutputs(path, i)
            end
        elseif(tissue == -1)
%             plotSobolAllOutputs(path)
        end
        
    elseif(selOut == 0)
        quit = 1;
    end
end
