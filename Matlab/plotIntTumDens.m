function plotIntTumDens(path, nTissue)
global nPar
global b color nfig shape;

intTumDens = load([path, '/morrisIntTumDens_', num2str(nTissue), '.res']);
intTumDens(:, 4) = sqrt(intTumDens(:, 1).^2 + intTumDens(:, 2).^2);
intTumDens(:, 3) = intTumDens(:, 1).^2 ./ intTumDens(:, 4);

figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(intTumDens, 1)
    scatter(intTumDens(i,1), intTumDens(i,2), 20, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)])], [0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex')
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(nTissue), ' - Integral tumor density'])
axis([0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)]), 0, 1.1 * max([intTumDens(:, 1); intTumDens(:, 2)])])
grid on
hold off

cIntTumDens = [num2cell(intTumDens), b'];
cIntTumDens = sortrows(cIntTumDens, 4);
figure(nfig);
nfig = nfig + 1;
bar(cell2mat(cIntTumDens(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cIntTumDens(:, 5));
ax.YGrid = 'on';
title(['Tissue ', num2str(nTissue), ' - Integral of tumor density'])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'location', 'northwest', 'interpreter', 'latex')
xtickangle(45)
end