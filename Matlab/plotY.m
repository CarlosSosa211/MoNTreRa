function plotY(path, n)
global nPar
global b color nfig shape;

Y = load([path, '/morrisY_', num2str(n), '.res']);
Y(:, 4) = sqrt(Y(:, 1).^2 + Y(:, 2).^2);
Y(:, 3) = Y(:, 1).^2 ./ Y(:, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : size(Y, 1)
    scatter(Y(i,1), Y(i,2), 200, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([Y(:, 1); Y(:, 2)])], [0, 1.1 * max([Y(:, 1); Y(:, 2)])], '--k')
legend(b, 'location', 'bestoutside', 'interpreter', 'latex',  'fontsize', 18)
xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)
title(['Range ', num2str(n), ' - Y'], 'fontsize', 20)
axis([0, 1.1 * max([Y(:, 1); Y(:, 2)]), 0, 1.1 * max([Y(:, 1); Y(:, 2)])])
grid on
hold off

cY = [num2cell(Y), b'];
cY = sortrows(cY, 4);
nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cY(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cY(:, 5));
ax.YGrid = 'on';
title(['Range ', num2str(n), ' - Y'], 'fontsize', 20)
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'location', 'northwest', 'interpreter', 'latex', 'fontsize', 20)
xtickangle(45)
end
