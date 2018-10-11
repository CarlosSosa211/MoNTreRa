function plotMeanAllOutputs(path)
global nPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b nfig;

tTissues = allTissues;
nTissues = length(tTissues);
timeTo95 = zeros(nPar, 2, nTissues);
timeTo99 = zeros(nPar, 2, nTissues);
tumDens = zeros(nPar, 2, nTissues);
intTumDens = zeros(nPar, 2, nTissues);

for i = 1:length(tTissues)
    timeTo95(:, 1:2, i) = load([path, '/morrisTimeTo95_', num2str(tTissues(i)), '.res']);
    timeTo95(:, 3, i) = sqrt(timeTo95(:, 1, i).^2 + timeTo95(:, 2, i).^2);
    timeTo99(:, 1:2, i) = load([path, '/morrisTimeTo99_', num2str(tTissues(i)), '.res']);
    timeTo99(:, 3, i) = sqrt(timeTo99(:, 1, i).^2 + timeTo99(:, 2, i).^2);
    tumDens(:, 1:2, i) = load([path, '/morrisTumDens_', num2str(tTissues(i)), '.res']);
    tumDens(:, 3, i) = sqrt(tumDens(:, 1, i).^2 + tumDens(:, 2, i).^2);
    intTumDens(:, 1:2, i) = load([path, '/morrisIntTumDens_', num2str(tTissues(i)), '.res']);
    intTumDens(:, 3, i) = sqrt(intTumDens(:, 1, i).^2 + intTumDens(:, 2, i).^2);
end

meanTimeTo95 = mean(timeTo95, 3);
stdTimeTo95 = std(timeTo95, [], 3);

meanTimeTo99 = mean(timeTo99, 3);
stdTimeTo99 = std(timeTo99, [], 3);

meanTumDens = mean(tumDens, 3);
stdTumDens = std(tumDens, [], 3);

meanIntTumDens = mean(intTumDens, 3);
stdIntTumDens = std(intTumDens, [], 3);

cOut = [num2cell(meanTimeTo95(:, 3)), num2cell(meanTimeTo99(:, 3))...
    num2cell(meanTumDens(:, 3)), num2cell(meanIntTumDens(:, 3))...
    num2cell(stdTimeTo95(:, 3)), num2cell(stdTimeTo99(:, 3))...
    num2cell(stdTumDens(:, 3)), num2cell(stdIntTumDens(:, 3)), b'];
cOut = sortrows(cOut, 1);

figure(nfig);
nfig = nfig + 1;
hold on
hBar = bar(cell2mat(cOut(:, 1:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOut(:, 9));
ax.YGrid = 'on';
title('21 tissues - Morris Euclidean distance from origin')
legend('Time to kill 95%', 'Time to kill 99%', 'location', 'northwest')
xtickangle(45)

ctr = zeros(nPar, 2);
ydt = zeros(nPar, 2);
for k = 1:2
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cOut(:, 5:8)), '.k')
hold off