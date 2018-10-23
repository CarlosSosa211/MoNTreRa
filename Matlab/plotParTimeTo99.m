function plotParTimeTo99(path, par, tissueSet)
global nPar;
global allTissues;
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global varRange;
global b color nfig shape;

switch tissueSet
    case 1
        tTissues = allTissues;
    case 2
        tTissues = densTissues;
    case 3
        tTissues = nonDensTissues;
    case 4
        tTissues = vascTissues;
    case 5
        tTissues = nonVascTissues;
    case 6
        tTissues = varRange;
end

nTissues = length(tTissues);
timeTo99 = zeros(nPar, 2, nTissues);

for i = 1:length(tTissues)
    timeTo99(:, :, i) = load([path, '/morrisTimeTo99_', num2str(tTissues(i)), '.res']);
end

figure(nfig);
nfig = nfig + 1;

hold on
colormap(jet)
if(tissueSet == 1)
    densAndVasc = ismember(tTissues, densTissues) & ismember(tTissues, vascTissues);
    densAndNonVasc = ismember(tTissues, densTissues) & ismember(tTissues, nonVascTissues);
    nonDensAndVasc = ismember(tTissues, nonDensTissues) & ismember(tTissues, vascTissues);
    nonDensAndNonVasc = ismember(tTissues, nonDensTissues) & ismember(tTissues, nonVascTissues);
    
    scatter(timeTo99(par, 1, densAndVasc), timeTo99(par, 2, densAndVasc), 200, 'r', 'filled', 'v')
    scatter(timeTo99(par, 1, densAndNonVasc), timeTo99(par, 2, densAndNonVasc), 200, 'r', 'filled', 'o')
    scatter(timeTo99(par, 1, nonDensAndVasc), timeTo99(par, 2, nonDensAndVasc), 200, 'g', 'filled', 'v')
    scatter(timeTo99(par, 1, nonDensAndNonVasc), timeTo99(par, 2, nonDensAndNonVasc), 200, 'g', 'filled', 'o')
    
else
    for i = 1:nTissues
        scatter(timeTo99(par, 1, i), timeTo99(par, 2, i), 200, color(i), 'filled', shape(mod(i, length(shape)) + 1))
    end
end

maxVal = 1.1 * max([reshape(timeTo99(par, 1, :), 1, []), reshape(timeTo99(par, 2, :), 1, [])]);
plot([0, maxVal], [0, maxVal], '--k')
hold off

xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)
title99 = strcat(string(b(par)), ' - Time to kill 99\% of tumor cells');
title(title99, 'interpreter', 'latex', 'fontsize', 20)
axis([0, maxVal, 0, maxVal])
grid on
switch tissueSet
    case 1
        legend({'Dense vascuralized tissues', 'Dense non-vascularized tissues',...
            'Non-dense vascularized tissues', 'Non-dense non-vascularized tissues'},...
            'location', 'northwest', 'fontsize', 20)
    case 6
        for i = 1:length(varRange)
            leg(i) = string(varRange(i));
        end
        legend(leg,'location', 'bestoutside', 'fontsize', 20)
end
end