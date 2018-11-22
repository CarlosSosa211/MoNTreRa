function plotAllOutputs(path, nTissue)
global nPar
global b nfig;

timeTo95 = load([path, '/morrisTimeTo95_', num2str(nTissue), '.res']);
timeTo95(:, 3) = sqrt(timeTo95(:, 1).^2 + timeTo95(:, 2).^2);
timeTo99 = load([path, '/morrisTimeTo99_', num2str(nTissue), '.res']);
timeTo99(:, 3) = sqrt(timeTo99(:, 1).^2 + timeTo99(:, 2).^2);
tumDens = load([path, '/morrisTumDens_', num2str(nTissue), '.res']);
tumDens(:, 3) = sqrt(tumDens(:, 1).^2 + tumDens(:, 2).^2);
intTumDens = load([path, '/morrisIntTumDens_', num2str(nTissue), '.res']);
intTumDens(:, 3) = sqrt(intTumDens(:, 1).^2 + intTumDens(:, 2).^2);

cOut = [num2cell(timeTo95(:, 3)), num2cell(timeTo99(:, 3))...
    num2cell(tumDens(:, 3)), num2cell(intTumDens(:, 3)), b'];
cOut = sortrows(cOut, 1);

figure(nfig);
nfig = nfig + 1;
bar(cell2mat(cOut(:, 1:2)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOut(:, 5));
ax.YGrid = 'on';
ylim([0, inf])
title(['Tissue ', num2str(nTissue), ' - All outputs'], 'fontsize', 20)
legend('Time to kill 95%', 'Time to kill 99%', 'Final tumor density',...
    'Integral of tumor density', 'location', 'northwest', 'fontsize', 20)
xtickangle(45)