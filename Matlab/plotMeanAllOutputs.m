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
endTreatTumDens = zeros(nPar, 2, nTissues);
intTumDens = zeros(nPar, 2, nTissues);

for i = 1:length(tTissues)
    timeTo95(:, 1:2, i) = load([path, '/morrisTimeTo95_', num2str(tTissues(i)), '.res']);
    timeTo95(:, 3, i) = sqrt(timeTo95(:, 1, i).^2 + timeTo95(:, 2, i).^2);
    timeTo99(:, 1:2, i) = load([path, '/morrisTimeTo99_', num2str(tTissues(i)), '.res']);
    timeTo99(:, 3, i) = sqrt(timeTo99(:, 1, i).^2 + timeTo99(:, 2, i).^2);
    endTreatTumDens(:, 1:2, i) = load([path, '/morrisEndTreatTumDens_', num2str(tTissues(i)), '.res']);
    endTreatTumDens(:, 3, i) = sqrt(endTreatTumDens(:, 1, i).^2 + endTreatTumDens(:, 2, i).^2);
    intTumDens(:, 1:2, i) = load([path, '/morrisIntTumDens_', num2str(tTissues(i)), '.res']);
    intTumDens(:, 3, i) = sqrt(intTumDens(:, 1, i).^2 + intTumDens(:, 2, i).^2);
end

meanTimeTo95 = mean(timeTo95, 3);
maxMeanTimeTo95_ = 1./max(meanTimeTo95);
meanTimeTo95 = meanTimeTo95 .* maxMeanTimeTo95_;
stdTimeTo95 = std(timeTo95, [], 3);
stdTimeTo95 = stdTimeTo95 .* maxMeanTimeTo95_;

meanTimeTo99 = mean(timeTo99, 3);
maxMeanTimeTo99_ = 1./max(meanTimeTo99);
meanTimeTo99 = meanTimeTo99 .* maxMeanTimeTo99_;
stdTimeTo99 = std(timeTo99, [], 3);
stdTimeTo99 = stdTimeTo99 .* maxMeanTimeTo99_;

meanEndTreatTumDens = mean(endTreatTumDens, 3);
maxMeanEndTreatTumDens_ = 1./max(meanEndTreatTumDens);
meanEndTreatTumDens = meanEndTreatTumDens .* maxMeanEndTreatTumDens_;
stdEndTreatTumDens = std(endTreatTumDens, [], 3);
stdEndTreatTumDens = stdEndTreatTumDens .* maxMeanEndTreatTumDens_;

meanIntTumDens = mean(intTumDens, 3);
maxMeanIntTumDens_ = 1./max(meanIntTumDens);
meanIntTumDens = meanIntTumDens .* maxMeanIntTumDens_;
stdIntTumDens = std(intTumDens, [], 3);
stdIntTumDens = stdIntTumDens .* maxMeanIntTumDens_;

cOut = [num2cell(meanTimeTo95(:, 3)), num2cell(meanTimeTo99(:, 3))...
    num2cell(meanEndTreatTumDens(:, 3)), num2cell(meanIntTumDens(:, 3))...
    num2cell(stdTimeTo95(:, 3)), num2cell(stdTimeTo99(:, 3))...
    num2cell(stdEndTreatTumDens(:, 3)), num2cell(stdIntTumDens(:, 3)), b'];
cOut = sortrows(cOut, 1);


figure(nfig);
nfig = nfig + 1;

hold on
hBar = bar(cell2mat(cOut(:, 1:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOut(:, 9), 'fontsize', 20);
ax.YGrid = 'on';
title('21 tissues - Morris Euclidean distance from origin', 'fontsize', 20)
xtickangle(45)

ctr = zeros(nPar, 4);
ydt = zeros(nPar, 4);
for k = 1:4
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cOut(:, 5:8)), '.k')
hold off
ylim([0, inf])
legend({'Time to kill 95%', 'Time to kill 99%', 'Tumor density at the end of treat.', 'Integral of tumor density'},...
    'location', 'northwest', 'fontsize', 20)
end