function plotSobolMeanOutput(path, selOut)
global nPar nTotPar selPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b colorbar colorbar2 nfig;
global fileNames outputNames;

tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), '...
    'vascularised (4) or non-vascularised (5) tissues]: ']);

switch tissueSet
    case 1
        tTissues = allTissues;
    case 2
        tTissues = densTissues;
    case 3
        tTissues = nonDensTissues;
    case 4
        tTissues = vascTissues;
    case 5
        tTissues = nonVascTissues;
end

nTissues = length(tTissues);
output = zeros(nTotPar, 3, nTissues);

for i = 1:length(tTissues)
    temp = load([path, '/sobol', char(fileNames(selOut)), '_'...
        num2str(tTissues(i)), '.res']);
    output(:, 1:2, i) = temp(2:end, :);
    output(:, 3, i) = (output(:, 2, i) - output(:, 1, i));
end

output = output(selPar == 1, :, :);

meanOutput = mean(output, 3);
stdOutput = std(output, [], 3);

cOutput = [num2cell(meanOutput), num2cell(stdOutput), b'...
    num2cell(colorbar), num2cell(colorbar2)];
cOutput = sortrows(cOutput, 2, 'descend');

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));

hold on
hBar = bar(cell2mat(cOutput(:, [1, 3])), 'stacked', 'faceColor', 'flat',...
    'Linewidth', 2);
% hBar = bar(cell2mat(cOutput(:, [1, 2])), 'faceColor', 'flat');
hBar(1).CData = cell2mat(cOutput(:, 8:10));
hBar(2).CData = cell2mat(cOutput(:, 11:13));

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOutput(:, 7));
ax.YGrid = 'on';
ax.FontSize = 42;

xpos = hBar(1).XData;
ypos = hBar(1).YData;
errorbar(xpos, ypos, cell2mat(cOutput(:, 4)), 'ok', 'Linewidth', 2)

xpos = hBar(1).XData;
ypos = hBar(1).YData + hBar(2).YData;
errorbar(xpos, ypos, cell2mat(cOutput(:, 5)), 'ok', 'Linewidth', 2)
hold off

ylim([0, 1])
title(['76 tissues - ', char(outputNames(selOut))], 'fontsize', 42)
legend({'Main', 'Interactions'}, 'location', 'northeast',...
    'interpreter', 'latex', 'fontsize', 42)
% legend({'SI', 'TSI'}, 'location', 'northeast',...
%     'interpreter', 'latex', 'fontsize', 42)
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