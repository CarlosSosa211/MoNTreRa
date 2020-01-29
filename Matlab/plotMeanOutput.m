function plotMeanOutput(path, selOut)
global nPar nImpPar nTotPar selPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color colorbar nfig shape;
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
output = zeros(nTotPar, 4, nTissues);

for i = 1:length(tTissues)
    output(:, 1:2, i) = load([path, '/morris', char(fileNames(selOut))...
        '_', num2str(tTissues(i)), '.res']);
    output(:, 4, i) = sqrt(output(:, 1, i).^2 + output(:, 2, i).^2);
    output(:, 3, i) = output(:, 1, i).^2 ./ output(:, 4, i);
end

output = output(selPar == 1, :, :);

meanOutput = mean(output, 3);
stdOutput = std(output, [], 3);

cOutput = [num2cell(meanOutput), num2cell(stdOutput), b'...
    num2cell(colorbar)];
cOutput = sortrows(cOutput, 4, 'descend');

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));

hold on
colormap(jet)
for i = nPar:-1:nPar - nImpPar + 1
    scatter(cell2mat(cOutput(i,1)), cell2mat(cOutput(i,2)), 500,...
        color(i), 'filled', shape(mod(i, length(shape)) + 1))
end

maxValMu = 1.1 * max(meanOutput(:, 1));
maxValSigma = 1.1 * max(meanOutput(:, 2));
maxVal = max([maxValMu, maxValSigma]);
plot([0, maxVal], [0, maxVal], '--k', 'Linewidth', 4)
hold off

ax = gca;
ax.FontSize = 42;
switch tissueSet
    case 1
        title(['21 tissues - ', char(outputNames(selOut))], 'fontsize', 42)
    case 2
        title(['11 dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 42)
    case 3
        title(['10 non-dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 42)
    case 4
        title(['11 vascularised tissues - ', char(outputNames(selOut))],...
            'fontsize', 42)
    case 5
        title(['10 non-vascularised tissues - '...
            char(outputNames(selOut))], 'fontsize', 42)
end
axis([0, maxValMu, 0, maxValSigma])
grid on
legend(flip(cOutput(nPar - nImpPar + 1:nPar, 9)), 'Location',...
    'bestoutside', 'Interpreter', 'Latex', 'fontsize', 18)
xlabel('\mu*', 'fontsize', 42)
ylabel('\sigma', 'fontsize', 42)

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
print(fig, [char(fileNames(selOut)), '.pdf'], '-dpdf', '-r0')

% nfig = nfig + 1;
% figure(nfig);
fig = figure('position', get(0, 'screensize'));

hold on
% hBar = bar(cell2mat(cOutput(end - nImpPar + 1:end, 3:4)), 'faceColor',...
%     'flat');
hBar = bar(cell2mat(cOutput(end - nImpPar + 1:end, 4)), 'faceColor',...
    'flat');
hBar.CData(:, :) = cell2mat(cOutput(end - nImpPar + 1:end, 10:12));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOutput(end - nImpPar + 1:end, 9), 'fontsize', 20);
ax.YGrid = 'on';
ax = gca;
ax.FontSize = 32;
ylabel('Euclidean distance to the origin')
switch tissueSet
    case 1
        title(['21 tissues - ', char(outputNames(selOut))], 'fontsize', 42)
    case 2
        title(['11 dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 42)
    case 3
        title(['10 non-dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 42)
    case 4
        title(['11 vascularised tissues - ', char(outputNames(selOut))],...
            'fontsize', 42)
    case 5
        title(['10 non-vascularised tissues - '...
            char(outputNames(selOut))], 'fontsize', 42)
end

% xpos = zeros(nImpPar, 2);
% ypos = zeros(nImpPar, 2);
% for k = 1:2
%     xpos(:, k) = hBar(1).XData + hBar(k).XOffset;
%     ypos(:, k) = hBar(k).YData;
% end

% errorbar(xpos, ypos, cell2mat(cOutput(end - nImpPar + 1:end, 7:8)),...
%     '.k', 'Linewidth', 2)

xpos = hBar.XData;
ypos = hBar.YData;

errorbar(xpos, ypos, cell2mat(cOutput(end - nImpPar + 1:end, 8)),...
    '.k', 'Linewidth', 2)

% dthresi =  find(strcmp(cOutput(end - nImpPar + 1:end,9), '$d_{thres}$'));
% if dthresi
%     dthrespos = 0.5 * (xpos(dthresi, 2) + xpos(dthresi + 1, 1));
%     plot([dthrespos, dthrespos],...
%         [0, max(ypos(:, 2) + cell2mat(cOutput(:,8)))], '--k',...
%         'Linewidth', 2)
% end

dthresi =  find(strcmp(cOutput(end - nImpPar + 1:end,9), '$d_{thres}$'));
if dthresi
    dthrespos = 0.5 * (xpos(dthresi - 1) + xpos(dthresi));
    plot([dthrespos, dthrespos],...
        [0, max(ypos' + cell2mat(cOutput(:,8)))], '--k',...
        'Linewidth', 2)
end
hold off

ylim([0, inf])
% legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$'...
%     '$\sqrt{\mu*^2 + \sigma^2}$'}, 'fontsize', 32, 'location',...
%     'northwest', 'interpreter', 'latex')
xtickangle(90)

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
print(fig, [char(fileNames(selOut)), 'Dist.pdf'], '-dpdf', '-r0')
end