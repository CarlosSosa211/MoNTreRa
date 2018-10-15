function plotTumDens(path, nTissue)
global nPar
global b color nfig shape;

tumDens = load([path, '/morrisTumDens_', num2str(nTissue), '.res']);
tumDens(:, 4) = sqrt(tumDens(:, 1).^2 + tumDens(:, 2).^2);
tumDens(:, 3) = tumDens(:, 1).^2 ./ tumDens(:, 4);

figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(tumDens, 1)
    scatter(tumDens(i,1), tumDens(i,2), 20, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])], [0, 1.1 * max([tumDens(:, 1); tumDens(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'interpreter', 'latex')
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(nTissue), ' - Final tumor density'])
axis([0,  1.1 * max([tumDens(:, 1); tumDens(:, 2)]), 0,  1.1 * max([tumDens(:, 1); tumDens(:, 2)])])
grid on
hold off

cTumDens = [num2cell(tumDens), b'];
cTumDens = sortrows(cTumDens, 4);
figure(nfig);
nfig = nfig + 1;
bar(cell2mat(cTumDens(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTumDens(:, 5));
ax.YGrid = 'on';
title(['Tissue ', num2str(nTissue), ' - Final tumor density'])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'location', 'northwest', 'interpreter', 'latex')
xtickangle(45)
end