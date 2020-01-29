function plotSobolOutput(path, n, selOut)
global nPar selPar
global b colorbar nfig;
global fileNames outputNames;

output = load([path, '/sobol', char(fileNames(selOut)), '_', num2str(n)...
    '.res']);
output = output(2:end, :);
output = output(selPar == 1, :);
output(:, 3) = (output(:, 2) - output(:, 1));

cOutput = [num2cell(output), b', num2cell(colorbar)];
cOutput = sortrows(cOutput, 2, 'descend');

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));

hBar = bar(cell2mat(cOutput(:, [1, 3])), 'stacked', 'faceColor',...
    'flat');
hBar(1).CData = cell2mat(cOutput(:, 5:7));
hBar(2).CData = cell2mat(cOutput(:, 5:7));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOutput(:, 4));
ax.YGrid = 'on';
ax.FontSize = 42;
ylim([0, 1])
title(['Tissue ', num2str(n), ' - ', char(outputNames(selOut))],...
    'fontsize', 42)
legend({'Main', 'Interactions'}, 'location', 'northeast',...
    'interpreter', 'latex', 'fontsize', 42)
xtickangle(45)

outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, ['sobol', char(fileNames(selOut)), '.pdf'], '-dpdf', '-r0')
end