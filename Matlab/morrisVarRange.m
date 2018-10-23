clear all
close all
global nPar
global varRange;
global b color nfig shape;

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
shape = ['o', 's', 'v', 'd'];

nfig = 1;
quit = 0;

% path = uigetdir('../../Carlos/Results');
path = '../../Carlos/Results/Morris2_VarTtum/';

varRange = load([path, 'varRange.dat']);
par = varRange(1) + 1;
varRange = varRange(2 : length(varRange));

while(~quit)
    output = input('Select an output [timeTo95 (1), timeTo99 (2), tumDens (3) or intTumDens (4)] or quit (0): ');
    switch output
        case 1
            sel = 1;
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTumDens(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTumDens(path, par, 6)
            end
            
        case 2
            sel = 1;
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTumDens(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTumDens(path, par, 6)
            end
            
        case 3
            sel = 1;
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTumDens(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTumDens(path, par, 6)
            end
            
        case 4
            sel = 1;
            if(sel == 0)
                for i = 1 : length(varRange)
                    plotTumDens(path, varRange(i))
                end
            elseif(sel == 1)
                plotParTumDens(path, par, 6)
            end
                     
        case 0
            quit = 1;
    end
end
