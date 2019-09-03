function plotMeanOutput(path, selOut)
global nPar nTotPar selPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color nfig shape;
global fileNames outputNames;

tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), '...
    'vascularized (4) or non-vascularized (5) tissues]: ']);

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

cOutput = [num2cell(meanOutput), num2cell(stdOutput), b'];
cOutput = sortrows(cOutput, 4);

nfig = nfig + 1;
figure(nfig);
hold on
colormap(jet)
for i = 1 : size(output, 1)
    scatter(meanOutput(i,1), meanOutput(i,2), 500, color(i), 'filled',...
        shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([output(:, 1); output(:, 2)])],...
    [0, 1.1 * max([output(:, 1); output(:, 2)])], '--k')
hold off

ax = gca;
ax.FontSize = 22;
switch tissueSet
    case 1
        title(['21 tissues - ', char(outputNames(selOut))], 'fontsize', 22)
    case 2
        title(['11 dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 22)
    case 3
        title(['10 non-dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 22)
    case 4
        title(['11 vascularized tissues - ', char(outputNames(selOut))],...
            'fontsize', 22)
    case 5
        title(['10 non-vascularized tissues - '...
            char(outputNames(selOut))], 'fontsize', 22)
end
axis([0, 1.1 * max([output(:, 1); output(:, 2)]),...
    0, 1.1 * max([output(:, 1); output(:, 2)])])
grid on
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex',...
    'fontsize', 16)
xlabel('\mu*', 'fontsize', 22)
ylabel('\sigma', 'fontsize', 22)

nfig = nfig + 1;
figure(nfig);
hold on
hBar = bar(cell2mat(cOutput(:, 3:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOutput(:, 9), 'fontsize', 22);
ax.YGrid = 'on';
ax = gca;
ax.FontSize = 22;
switch tissueSet
    case 1
        title(['21 tissues - ', char(outputNames(selOut))], 'fontsize', 20)
    case 2
        title(['11 dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 22)
    case 3
        title(['10 non-dense tissues - ', char(outputNames(selOut))],...
            'fontsize', 22)
    case 4
        title(['11 vascularized tissues - ', char(outputNames(selOut))],...
            'fontsize', 22)
    case 5
        title(['10 non-vascularized tissues - '...
            char(outputNames(selOut))], 'fontsize', 22)
end

xpos = zeros(nPar, 2);
ypos = zeros(nPar, 2);
for k = 1:2
    xpos(:, k) = hBar(1).XData + hBar(k).XOffset;
    ypos(:, k) = hBar(k).YData;
end

errorbar(xpos, ypos, cell2mat(cOutput(:, 7:8)), '.k')
hold off
ylim([0, inf])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$'...
    '$\sqrt{\mu*^2 + \sigma^2}$'}, 'fontsize', 20, 'location',...
    'northwest', 'interpreter', 'latex')
xtickangle(45)
end
