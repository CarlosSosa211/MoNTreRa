function plotEndTreatTumDens(path, n)
global nPar
global b color nfig shape;

endTreatTumDens = load([path, '/morrisEndTreatTumDens_', num2str(n), '.res']);
endTreatTumDens(:, 4) = sqrt(endTreatTumDens(:, 1).^2 + endTreatTumDens(:, 2).^2);
endTreatTumDens(:, 3) = endTreatTumDens(:, 1).^2 ./ endTreatTumDens(:, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : size(endTreatTumDens, 1)
    scatter(endTreatTumDens(i,1), endTreatTumDens(i,2), 200, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)])], [0, 1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'interpreter', 'latex', 'fontsize', 18)
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(n), ' - Tumor density at the end of treat.'], 'fontsize', 20)
axis([0,  1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)]), 0,  1.1 * max([endTreatTumDens(:, 1); endTreatTumDens(:, 2)])])
grid on
hold off

cEndTreatTumDens = [num2cell(endTreatTumDens), b'];
cEndTreatTumDens = sortrows(cEndTreatTumDens, 4);
nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cEndTreatTumDens(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cEndTreatTumDens(:, 5));
ax.YGrid = 'on';
title(['Tissue ', num2str(n), ' - Tumor density at the end of treat.'], 'fontsize', 20)
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'location', 'northwest', 'interpreter', 'latex', 'fontsize', 20)
xtickangle(45)
end