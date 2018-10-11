function plotTimeTo95(path, nTissue)
global nPar
global b color nfig shape;

timeTo95 = load([path, '/morrisTimeTo95_', num2str(nTissue), '.res']);
figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(timeTo95, 1)
    scatter(timeTo95(i,1), timeTo95(i,2), 20, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])], [0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex')
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(nTissue), ' - Time to kill 95% of tumor cells'])
axis([0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)]), 0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])])
grid on
hold off

distTimeTo95 = sqrt(timeTo95(:, 1).^2 + timeTo95(:, 2).^2);
cTimeTo95 = [num2cell(timeTo95), num2cell(distTimeTo95), b'];
cTimeTo95 = sortrows(cTimeTo95, 3);
figure(nfig);
nfig = nfig + 1;
bar(cell2mat(cTimeTo95(:, 1:3)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTimeTo95(:, 4));
ax.YGrid = 'on';
title(['Tissue ', num2str(nTissue), ' - Time to kill 95% of tumor cells'])
legend('\mu*', '\sigma', 'dist', 'location', 'northwest') 
xtickangle(45)
end
