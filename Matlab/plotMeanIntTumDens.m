function plotMeanIntTumDens(path)
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
intTumDens = zeros(nPar, 4, nTissues);

for i = 1:length(tTissues)
    intTumDens(:, 1:2, i) = load([path, '/morrisIntTumDens_', num2str(tTissues(i)), '.res']);
    intTumDens(:, 4, i) = sqrt(intTumDens(:, 1, i).^2 + intTumDens(:, 2, i).^2);
    intTumDens(:, 3, i) = intTumDens(:, 1, i).^2 ./ intTumDens(:, 4, i);
end

meanIntTumDens = mean(intTumDens, 3);
stdIntTumDens = std(intTumDens, [], 3);

cIntTumDens = [num2cell(meanIntTumDens), num2cell(stdIntTumDens), b'];
cIntTumDens = sortrows(cIntTumDens, 4);

figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(intTumDens, 1)
    scatter(meanIntTumDens(i,1), meanIntTumDens(i,2), 500, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)])], [0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)])], '--k')
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
axis([0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)]), 0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)])])
grid on
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex', 'fontsize', 18)
xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)

figure(nfig);
nfig = nfig + 1;

hold on
hBar = bar(cell2mat(cIntTumDens(:, 3:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cIntTumDens(:, 9), 'fontsize', 20);
ax.YGrid = 'on';
switch tissueSet
    case 1
        title('21 tissues - Integral of tumor density', 'fontsize', 20)
    case 2
        title('11 dense tissues - Integral of tumor density', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Integral of tumor density', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Integral of tumor density', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Integral of tumor density', 'fontsize', 20)
end
xtickangle(45)

ctr = zeros(nPar, 2);
ydt = zeros(nPar, 2);
for k = 1:2
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cIntTumDens(:, 7:8)), '.k')
hold off
ylim([0, inf])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'fontsize', 20, 'location', 'northwest', 'interpreter', 'latex')
end
