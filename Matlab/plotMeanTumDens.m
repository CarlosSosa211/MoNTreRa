function plotMeanTumDens(path)
global nPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color nfig shape;

tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), ' ...
    'vascularized (4) or non-vascularized (5) tissues]: ']);

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
tumDens = zeros(nPar, 4, nTissues);

for i = 1:length(tTissues)
    tumDens(:, 1:2, i) = load([path, '/morrisTumDens_', num2str(tTissues(i)), '.res']);
    tumDens(:, 4, i) = sqrt(tumDens(:, 1, i).^2 + tumDens(:, 2, i).^2);
    tumDens(:, 3, i) = tumDens(:, 1, i).^2 ./ tumDens(:, 4, i);
end

meanTumDens = mean(tumDens, 3);
stdTumDens = std(tumDens, [], 3);

cTumDens = [num2cell(meanTumDens), num2cell(stdTumDens), b'];
cTumDens = sortrows(cTumDens, 4);

figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(tumDens, 1)
    scatter(meanTumDens(i,1), meanTumDens(i,2), 500, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])], [0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])], '--k')
hold off

switch tissueSet
    case 1
        title('21 tissues - Final tumor density', 'fontsize', 20)
    case 2
        title('11 dense tissues - Final tumor density', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Final tumor density', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Final tumor density', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Final tumor density', 'fontsize', 20)
        
end
axis([0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)]), 0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])])
grid on
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex', 'fontsize', 16)
xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)

figure(nfig);
nfig = nfig + 1;

hold on
hBar = bar(cell2mat(cTumDens(:, 3:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTumDens(:, 9), 'fontsize', 20);
ax.YGrid = 'on';
switch tissueSet
    case 1
        title('21 tissues - Final tumor density', 'fontsize', 20)
    case 2
        title('11 dense tissues - Final tumor density', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Final tumor density', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Final tumor density', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Final tumor density', 'fontsize', 20)
end

ctr = zeros(nPar, 2);
ydt = zeros(nPar, 2);
for k = 1:2
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cTumDens(:, 7:8)), '.k')
hold off
ylim([0, inf])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'fontsize', 20, 'location', 'northwest', 'interpreter', 'latex')
xtickangle(45)
end
