function plotParTumDens(path, par, tissueSet)
global nPar;
global allTissues;
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
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
end

nTissues = length(tTissues);
tumDens = zeros(nPar, 2, nTissues);

for i = 1:length(tTissues)
    tumDens(:, :, i) = load([path, '/morrisTumDens_', num2str(tTissues(i)), '.res']);
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
    
    scatter(tumDens(par, 1, densAndVasc), tumDens(par, 2, densAndVasc), 20, 'r', 'filled', 'v')
    scatter(tumDens(par, 1, densAndNonVasc), tumDens(par, 2, densAndNonVasc), 20, 'r', 'filled', 'o')
    scatter(tumDens(par, 1, nonDensAndVasc), tumDens(par, 2, nonDensAndVasc), 20, 'g', 'filled', 'v')
    scatter(tumDens(par, 1, nonDensAndNonVasc), tumDens(par, 2, nonDensAndNonVasc), 20, 'g', 'filled', 'o')
    
    legend('Dense vascuralized tissues', 'Dense non-vascularized tissues',...
        'Non-dense vascularized tissues', 'Non-dense non-vascularized tissues',...
        'Location', 'northwest')
else
    for i = 1:nTissues
        scatter(tumDens(par, 1, i), tumDens(par, 2, i), 20, color(i), 'filled', shape(mod(i, length(shape)) + 1))
    end
end
plot([0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])], [0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])], '--k')
xlabel('\mu*')
ylabel('\sigma')
titleDens = strcat(string(b(par)), ' - Final tumor density');
title(titleDens, 'Interpreter', 'Latex')
axis([0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)]), 0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])])
grid on
hold off