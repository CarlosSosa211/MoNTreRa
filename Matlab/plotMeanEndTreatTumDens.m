function plotMeanEndTreatTumDens(path)
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
endTreatTumDens = zeros(nPar, 4, nTissues);

for i = 1:length(tTissues)
    endTreatTumDens(:, 1:2, i) = load([path, '/morrisEndTreatTumDens_', num2str(tTissues(i)), '.res']);
    endTreatTumDens(:, 4, i) = sqrt(endTreatTumDens(:, 1, i).^2 + endTreatTumDens(:, 2, i).^2);
    endTreatTumDens(:, 3, i) = endTreatTumDens(:, 1, i).^2 ./ endTreatTumDens(:, 4, i);
end

meanEndTreatTumDens = mean(endTreatTumDens, 3);
stdEndTreatTumDens = std(endTreatTumDens, [], 3);

cEndTreatTumDens = [num2cell(meanEndTreatTumDens), num2cell(stdEndTreatTumDens), b'];
cEndTreatTumDens = sortrows(cEndTreatTumDens, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : size(endTreatTumDens, 1)
    scatter(meanEndTreatTumDens(i,1), meanEndTreatTumDens(i,2), 500, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)])], [0, 1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)])], '--k')
hold off

switch tissueSet
    case 1
        title('21 tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 2
        title('11 dense tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Tumor density at the end of treat.', 'fontsize', 20)
        
end
axis([0, 1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)]), 0, 1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)])])
grid on
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex', 'fontsize', 16)
xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)

nfig = nfig + 1;
figure(nfig);
hold on
hBar = bar(cell2mat(cEndTreatTumDens(:, 3:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cEndTreatTumDens(:, 9), 'fontsize', 20);
ax.YGrid = 'on';
switch tissueSet
    case 1
        title('21 tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 2
        title('11 dense tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Tumor density at the end of treat.', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Tumor density at the end of treat.', 'fontsize', 20)
end

ctr = zeros(nPar, 2);
ydt = zeros(nPar, 2);
for k = 1:2
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cEndTreatTumDens(:, 7:8)), '.k')
hold off
ylim([0, inf])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'fontsize', 20, 'location', 'northwest', 'interpreter', 'latex')
xtickangle(45)
end
