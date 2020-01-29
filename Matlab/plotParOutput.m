function plotParOutput(path, par, tissueSet, selOut)
global nPar;
global allTissues;
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global varRange;
global b color nfig shape;
global fileNames outputNames;

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
output = zeros(nPar, 3, nTissues);

initTumDens = load([path, '/initTumDens.dat']);
initVascDens = load([path, '/initVascDens.dat']);
initTumDens = initTumDens(tTissues);
initVascDens = initVascDens(tTissues);
[initTumDens, tissuesTumDens] = sort(initTumDens);
[initVascDens, tissuesVascDens] = sort(initVascDens);

for i = 1:length(tTissues)
    output(:, 1:2, i) = load([path, '/morris', char(fileNames(selOut))...
        '_', num2str(tTissues(i)), '.res']);
    output(:, 3, i) = sqrt(output(:, 1, i).^2 + output(:, 2, i).^2);
end

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));
hold on
colormap(jet)
if(tissueSet == 1)
    densAndVasc = ismember(tTissues, densTissues) &...
        ismember(tTissues, vascTissues);
    densAndNonVasc = ismember(tTissues, densTissues) &...
        ismember(tTissues, nonVascTissues);
    nonDensAndVasc = ismember(tTissues, nonDensTissues) &...
        ismember(tTissues, vascTissues);
    nonDensAndNonVasc = ismember(tTissues, nonDensTissues) &...
        ismember(tTissues, nonVascTissues);
    
    scatter(output(par, 1, densAndVasc), output(par, 2, densAndVasc),...
        400, 'r', 'filled', 'v')
    scatter(output(par, 1, densAndNonVasc),...
        output(par, 2, densAndNonVasc), 400, 'r', 'filled', 'o')
    scatter(output(par, 1, nonDensAndVasc),...
        output(par, 2, nonDensAndVasc), 400, 'g', 'filled', 'v')
    scatter(output(par, 1, nonDensAndNonVasc),...
        output(par, 2, nonDensAndNonVasc), 400, 'g', 'filled', 'o')
else
    for i = 1:nTissues
        scatter(output(par, 1, i), output(par, 2, i), 400, color(i),...
            'filled', shape(mod(i, length(shape)) + 1))
    end
end

maxValMu = 1.1 * max(reshape(output(par, 1, :), 1, []));
maxValSigma = 1.1 * max(reshape(output(par, 2, :), 1, []));
maxVal = max([maxValMu, maxValSigma]);
plot([0, maxVal], [0, maxVal], '--k')
hold off

xlabel('\mu*', 'fontsize', 42)
ylabel('\sigma', 'fontsize', 42)
ax = gca;
ax.FontSize = 42;
titleDens = strcat(string(b(par)), ' - ', char(outputNames(selOut)));
title(titleDens, 'interpreter', 'latex', 'fontsize', 42)
axis([0, maxValMu, 0, maxValSigma])
grid on
switch tissueSet
    case 1
        legend({'Dense well-vascuralized tissues'...
            'Dense poorly-vascularized tissues'...
            'Non-dense well-vascularized tissues'...
            'Non-dense poorly-vascularized tissues'},...
            'location', 'northwest', 'fontsize', 42)
    case 6
        for i = 1:length(varRange)
            leg(i) = string(varRange(i));
        end
        legend(leg,'location', 'bestoutside', 'fontsize', 42)
end
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, strcat(string(b(par)), char(fileNames(selOut)),...
    '.pdf'), '-dpdf', '-r0')

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));
plot(initVascDens',  permute(output(par, 3, tissuesVascDens),...
    [2, 3, 1]), '-o', 'Linewidth', 6, 'Markersize', 10)
ax = gca;
ax.FontSize = 42;
grid on
title(titleDens, 'interpreter', 'latex', 'Fontsize', 42)
xlabel('Initial vascular density (%)', 'FontSize', 42)
ylabel('Euclidean distance', 'FontSize', 42)
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, strcat(string(b(par)), char(fileNames(selOut)),...
    'VascDens.pdf'), '-dpdf', '-r0')

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));
plot(initTumDens',  permute(output(par, 3, tissuesTumDens),...
    [2, 3, 1]), '-o', 'Linewidth', 6, 'Markersize', 10)
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
ax.FontSize = 42;
grid on
title(titleDens, 'interpreter', 'latex', 'Fontsize', 42)
xlabel('Initial tumour density (%)', 'FontSize', 42)
ylabel('Euclidean distance', 'FontSize', 42)

% set(fig, 'units', 'inches')
% pos = get(fig, 'position');
% set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
%     'paperSize', [pos(3), pos(4)])
% print(fig, strcat(string(b(par)), char(fileNames(selOut)),...
%     'TumDens.pdf'), '-dpdf', '-r0')
end


