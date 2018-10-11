function plotTimeTo99(path, nTissue)
global nPar
global b color nfig shape;

timeTo99 = load([path, '/morrisTimeTo99_', num2str(nTissue), '.res']);
figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(timeTo99, 1)
    scatter(timeTo99(i,1), timeTo99(i,2), 20, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])], [0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex')
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(nTissue), ' - Time to kill 99% of tumor cells'])
axis([0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)]), 0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])])
grid on
hold off

distTimeTo99 = sqrt(timeTo99(:, 1).^2 + timeTo99(:, 2).^2);
cTimeTo99 = [num2cell(timeTo99), num2cell(distTimeTo99), b'];
cTimeTo99 = sortrows(cTimeTo99, 3);
figure(nfig);
nfig = nfig + 1;
bar(cell2mat(cTimeTo99(:, 1:3)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTimeTo99(:, 4));
ax.YGrid = 'on';
title(['Tissue ', num2str(nTissue), ' - Time to kill 99% of tumor cells'])
legend('\mu*', '\sigma', 'dist', 'location', 'northwest') 
xtickangle(45)
end