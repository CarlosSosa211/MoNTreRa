function plotTimeTo95(path, nTissue)
global nPar
global b color nfig shape;

timeTo95 = load([path, '/morrisTimeTo95_', num2str(nTissue), '.res']);
timeTo95(:, 4) = sqrt(timeTo95(:, 1).^2 + timeTo95(:, 2).^2);
timeTo95(:, 3) = timeTo95(:, 1).^2 ./ timeTo95(:, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : size(timeTo95, 1)
    scatter(timeTo95(i,1), timeTo95(i,2), 200, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])], [0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])], '--k')
legend(b, 'location', 'bestoutside', 'interpreter', 'latex',  'fontsize', 18)
xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)
title(['Tissue ', num2str(nTissue), ' - Time to kill 95% of tumor cells'], 'fontsize', 20)
axis([0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)]), 0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])])
grid on
hold off

cTimeTo95 = [num2cell(timeTo95), b'];
cTimeTo95 = sortrows(cTimeTo95, 4);
nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cTimeTo95(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTimeTo95(:, 5));
ax.YGrid = 'on';
title(['Tissue ', num2str(nTissue), ' - Time to kill 95% of tumor cells'], 'fontsize', 20)
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'location', 'northwest', 'interpreter', 'latex', 'fontsize', 20)
xtickangle(45)
end
