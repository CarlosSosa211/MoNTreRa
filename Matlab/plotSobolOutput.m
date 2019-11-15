function plotSobolOutput(path, n, selOut)
global nPar selPar
global b nfig;
global fileNames outputNames;

output = load([path, '/sobol', char(fileNames(selOut)), '_', num2str(n)...
    '.res']);
output = output(selPar == 1, :);
output(:, 3) = (output(:, 2) - output(:, 1));

cOutput = [num2cell(output), b'];
cOutput = sortrows(cOutput, 2);

nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cOutput(:, [1, 3])), 'stacked')
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOutput(:, 4));
ax.YGrid = 'on';
ax.FontSize = 20;
ylim([0, 1])
title(['Tissue ', num2str(n), ' - ', char(outputNames(selOut))],...
    'fontsize', 22)
legend({'Main', 'Interactions'}, 'location', 'northwest',...
    'interpreter', 'latex', 'fontsize', 22)
xtickangle(45)
end