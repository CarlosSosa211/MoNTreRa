function plotOutput(path, n, selOut)
global nPar
global b color nfig shape;
global fileNames outputNames;

output = load([path, '/morris', char(fileNames(selOut)), '_', num2str(n)...
    '.res']);
output(:, 4) = sqrt(output(:, 1).^2 + output(:, 2).^2);
output(:, 3) = output(:, 1).^2 ./ output(:, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : nPar
    scatter(output(i,1), output(i,2), 200, color(i), 'filled',...
        shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([output(:, 1); output(:, 2)])],...
    [0, 1.1 * max([output(:, 1); output(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'interpreter', 'latex',...
    'fontsize', 16)
xlabel('\mu*')
ylabel('\sigma')
title(['Tissue ', num2str(n), ' - ', char(outputNames(selOut))],...
    'fontsize', 20)
axis([0,  1.1 * max([output(:, 1); output(:, 2)]),...
    0,  1.1 * max([output(:, 1); output(:, 2)])])
grid on
hold off

cOutput = [num2cell(output), b'];
cOutput = sortrows(cOutput, 4);
nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cOutput(:, 3:4)))
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOutput(:, 5));
ax.YGrid = 'on';
title(['Tissue ', num2str(n), ' - ', char(outputNames(selOut))],...
    'fontsize', 20)
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$'...
    '$\sqrt{\mu*^2 + \sigma^2}$'}, 'location', 'northwest',...
    'interpreter', 'latex', 'fontsize', 20)
xtickangle(45)
end