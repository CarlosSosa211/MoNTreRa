clear all
close all
global nPar
global varRange;
global b color nfig shape;
global fileNames outputNames;

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

nfig = 0;
quit = 0;

% path = uigetdir('../../Carlos/Results');
% path = '../../Carlos/Results/MorrisVarToyModel_Par2/';
path = '../../Carlos/Results/MorrisVarRangeCluster/';

fileNames = {'EndTreatTumDens', '3MonTumDens', 'RecTumDens',...
    'FinTumVol', 'IntTumDens', 'TimeTo95', 'TimeTo99'...
    'RecTime'};
outputNames = {'Tumour density at the end of treat.', 'Tumour density 3 months after the end of treat.'...
    'Tumour density at recurrence', 'Final tumour volume', 'Integral of tumour density'...
    'Time to kill 95\% of tumour cells', 'Time to kill 99\% of tumour cells', 'Recurrence time'};
nOut = 8;

varRange = load([path, 'morrisVarRange.res']);
par = varRange(1) + 1;
varRange = varRange(2 : length(varRange));

% nPar = 5;
% b = {'$x_{0}$', '$x_{1}$', '$x_{2}$', '$x_{3}$', '$x_{4}$'};
shape = ['o', 's', 'v', 'd'];

while(~quit)
    selOut = input(['Select an output [endTreaTumDens (1), 3MonTumDens (2), recTumDens (3), tumVol (4),\n'...
        'intTumDens (5), timeTo95 (6), timeTo99 (7) or recTime (8)] or quit (0): ']);
    
    if(selOut >= 1 && selOut <= 8)
        sel = input('Select the studied parameter (1) or all of them (2) or: ');
        if(sel == 1)
            color = linspace(0, 1, length(varRange));
            plotParOutput(path, par, 6, selOut)
        elseif(sel == 2)
            for i = 1 : length(varRange)
                color = linspace(0, 1, nPar);
                plotOutput(path, varRange(i), selOut)
            end
        end
        
    elseif (selOut == -2)
        if(sel == 1)
            color = linspace(0, 1, length(varRange));
            plotParY(path, par)
        elseif(sel == 2)
            for i = 1 : length(varRange)
                color = linspace(0, 1, nPar);
                plotY(path, i - 1)
            end
        end
        
    elseif(selOut == 0)
        quit = 1;
    end
end
