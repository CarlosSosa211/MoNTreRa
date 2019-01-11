function plotTimeTo99(path, nTissue)
global nPar
global b color nfig shape;

timeTo99 = load([path, '/morrisTimeTo99_', num2str(nTissue), '.res']);
timeTo99(:, 4) = sqrt(timeTo99(:, 1).^2 + timeTo99(:, 2).^2);
timeTo99(:, 3) = timeTo99(:, 1).^2 ./ timeTo99(:, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : nPar
    scatter(timeTo99(i,1), timeTo99(i,2), 200, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])], [0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'interpreter', 'latex', 'fontsize', 18)
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(nTissue), ' - Time to kill 99% of tumor cells'], 'fontsize', 20)
axis([0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)]), 0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])])
grid on
hold off

cTimeTo99 = [num2cell(timeTo99), b'];
cTimeTo99 = sortrows(cTimeTo99, 4);
nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cTimeTo99(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTimeTo99(:, 5));
ax.YGrid = 'on';
title(['Tissue ', num2str(nTissue), ' - Time to kill 99% of tumor cells'])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'location', 'northwest', 'interpreter', 'latex', 'fontsize', 20)
xtickangle(45)
end