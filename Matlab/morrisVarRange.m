clear all
close all
global nPar
global varRange;
global b color nfig shape;

% nPar = 34;
% b = {'$T_{tum}$', '$T_{heal}$', '$T_{end}$', '$\bar{v}$', '$\alpha_{heal}$'...
%     '$\alpha/\beta_{heal}$', '$\alpha_{tumG1}$', '$\alpha/\beta_{tumG1}$'...
%     '$\alpha_{tumS}$', '$\alpha/\beta_{tumS}$' '$\alpha_{tumG2}$'...
%     '$\alpha/\beta_{tumG2}$', '$\alpha_{tumM}$', '$\alpha/\beta_{tumM}$'...
%     '$\alpha_{tumG0}$', '$\alpha/\beta_{tumG0}$', '$\alpha_{preEnd}$'...
%     '$\alpha/\beta_{preEnd}$', '$\alpha_{neoEnd}$', '$\alpha/\beta_{neoEnd}$'...
%     '$d_{thres}$', '$T_{arrest}$', '$pO_2^{nec}$', '$d$', '$D^{O_2}$'...
%     '$V_{max}^{O_2}$', '$K_M^{O_2}$', '$D^{VEGF}$', '$V_{max}^{VEGF}$'...
%     '$K_M^{VEGF}$', '$pO_2^{preEnd}$', '$pO_2^{neoEnd}$', '$pO_2^{hyp}$', '$v^{hyp}$'};

nfig = 0;
quit = 0;

% path = uigetdir('../../Carlos/Results');
path = '../../Carlos/Results/MorrisVarToyModel_Par2/';

varRange = load([path, 'morrisVarRange.res']);
par = varRange(1) + 1;
varRange = varRange(2 : length(varRange));

nPar = 5;
b = {'$x_{0}$', '$x_{1}$', '$x_{2}$', '$x_{3}$', '$x_{4}$'};
color = linspace(0, 1, length(varRange));
shape = ['o', 's', 'v', 'd'];

while(~quit)
    output = input('Select an output [timeTo95 (1), timeTo99 (2), tumDens (3) or intTumDens (4)] or quit (0): ');
    sel = 0;
    switch output
        case 1
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTimeTo95(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTimeTo95(path, par, 6)
            end
            
        case 2
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTimeTo99(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTimeTo99(path, par, 6)
            end
            
        case 3
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTumDens(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTumDens(path, par, 6)
            end
            
        case 4
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotIntTumDens(path, varRange(i))
                end
            elseif(sel == 1)
                plotParIntTumDens(path, par, 6)
            end
            
        case -1
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotY(path, i - 1)
                end
            elseif(sel == 1)
                plotParY(path, par)
            end
            
        case 0
            quit = 1;
    end
end