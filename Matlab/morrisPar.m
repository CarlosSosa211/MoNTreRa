clear all
close all
global nPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color nfig shape;

densTissues = [1, 2, 5, 6, 8, 9, 11, 12, 19, 20, 21];
nonDensTissues = [3, 4, 7, 10, 13, 14, 15, 16, 17, 18];
vascTissues = [4, 7, 8, 10, 11, 12, 13, 14, 16, 18, 20];
nonVascTissues = [1, 2, 3, 5, 6, 9, 15, 17, 19, 21];
allTissues = [densTissues, nonDensTissues];

nDensTissues = length(densTissues);
nVascTissues = length(vascTissues);
nTissues = length(allTissues); 

nPar = 33;
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

path = uigetdir('../../Carlos/Results');
while(~quit)
    output = input(['Select an output [timeTo95 (1), timeTo99 (2), tumDens (3), intTumDens (4)] '...
        'or quit (0): ']);
    switch output
        case 1
            par = input(['Select one parameter (from 1 to ', num2str(nPar)...
                ') or all of them (0): ']);
            tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), ' ...
            'vascularized (4) or non-vascularized (5) tissues: ']);
            if(par >= 1 && par <= nPar)
                plotParTimeTo95(path, par, tissueSet)
            elseif(par == 0)
                for i = 1:nPar
                    plotParTimeTo95(path, i, tissueSet)
                end
            end
            
        case 2
            par = input(['Select one parameter (from 1 to ', num2str(nPar)...
                ') or all of them (0): ']);
            tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), ' ...
            'vascularized (4) or non-vascularized (5) tissues: ']);
            if(par >= 1 && par <= nPar)
                plotParTimeTo99(path, par, tissueSet)
            elseif(par == 0)
                for i = 1:nPar
                    plotParTimeTo99(path, i, tissueSet)
                end
            end
            
        case 3
            par = input(['Select one parameter (from 1 to ', num2str(nPar)...
                ') or all of them (0): ']);
            tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), ' ...
            'vascularized (4) or non-vascularized (5) tissues: ']);
            if(par >= 1 && par <= nPar)
                plotParTumDens(path, par, tissueSet)
            elseif(par == 0)
                for i = 1:nPar
                    plotParTumDens(path, i, tissueSet)
                end
            end
            
        otherwise
            quit = 1;
    end
end
