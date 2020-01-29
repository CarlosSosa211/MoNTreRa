function plotMeanAllOutputs(path)
global nPar nImpPar nTotPar selPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b nfig;

tTissues = allTissues;
nTissues = length(tTissues);
endTreatTumDens = zeros(nTotPar, 3, nTissues);
threeMonTumDens = zeros(nTotPar, 3, nTissues);
killed99 = zeros(nTotPar, 3, nTissues);
intTumDens = zeros(nTotPar, 3, nTissues);

for i = 1:length(tTissues)
    endTreatTumDens(:, 1:2, i) = load([path, '/morrisEndTreatTumDens_'...
        num2str(tTissues(i)), '.res']);
    endTreatTumDens(:, 3, i) = sqrt(endTreatTumDens(:, 1, i).^2 +...
        endTreatTumDens(:, 2, i).^2);
    threeMonTumDens(:, 1:2, i) = load([path, '/morris3MonTumDens_'...
        num2str(tTissues(i)), '.res']);
    threeMonTumDens(:, 3, i) = sqrt(threeMonTumDens(:, 1, i).^2 +...
        threeMonTumDens(:, 2, i).^2);
    killed99(:, 1:2, i) = load([path, '/morrisKilled99_'...
        num2str(tTissues(i)), '.res']);
    killed99(:, 3, i) = sqrt(killed99(:, 1, i).^2 + killed99(:, 2, i).^2);
    intTumDens(:, 1:2, i) = load([path, '/morrisIntTumDens_'...
        num2str(tTissues(i)), '.res']);
    intTumDens(:, 3, i) = sqrt(intTumDens(:, 1, i).^2 +...
        intTumDens(:, 2, i).^2);
end

endTreatTumDens = endTreatTumDens(selPar == 1, :, :);
meanEndTreatTumDens = mean(endTreatTumDens, 3);
maxMeanEndTreatTumDens_ = 1./max(meanEndTreatTumDens);
meanEndTreatTumDens = meanEndTreatTumDens .* maxMeanEndTreatTumDens_;
stdEndTreatTumDens = std(endTreatTumDens, [], 3);
stdEndTreatTumDens = stdEndTreatTumDens .* maxMeanEndTreatTumDens_;

threeMonTumDens = threeMonTumDens(selPar == 1, :, :);
mean3MonTumDens = mean(threeMonTumDens, 3);
maxMean3MonTumDens_ = 1./max(mean3MonTumDens);
mean3MonTumDens = mean3MonTumDens .* maxMean3MonTumDens_;
std3MonTumDens = std(threeMonTumDens, [], 3);
std3MonTumDens = std3MonTumDens .* maxMean3MonTumDens_;

killed99 = killed99(selPar == 1, :, :);
meanKilled99 = mean(killed99, 3);
maxMeanKilled99_ = 1./max(meanKilled99);
meanKilled99 = meanKilled99 .* maxMeanKilled99_;
stdKilled99 = std(killed99, [], 3);
stdKilled99 = stdKilled99 .* maxMeanKilled99_;

intTumDens = intTumDens(selPar == 1, :, :);
meanIntTumDens = mean(intTumDens, 3);
maxMeanIntTumDens_ = 1./max(meanIntTumDens);
meanIntTumDens = meanIntTumDens .* maxMeanIntTumDens_;
stdIntTumDens = std(intTumDens, [], 3);
stdIntTumDens = stdIntTumDens .* maxMeanIntTumDens_;

cOut = [num2cell(meanEndTreatTumDens(:, 3))...
    num2cell(mean3MonTumDens(:, 3)), num2cell(meanKilled99(:, 3))...
    num2cell(meanIntTumDens(:, 3)), num2cell(stdEndTreatTumDens(:, 3))...
    num2cell(std3MonTumDens(:, 3)), num2cell(stdKilled99(:, 3))...
    num2cell(stdIntTumDens(:, 3)), b'];
cOut = sortrows(cOut, 1);


nfig = nfig + 1;
figure(nfig);
hold on
hBar = bar(cell2mat(cOut(end - nImpPar + 1:end, 1:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOut(end - nImpPar + 1:end, 9), 'fontsize', 26);
ax.YGrid = 'on';
ax.FontSize = 26;
title('21 tissues - Morris Euclidean distance from origin', 'fontsize', 26)
xtickangle(45)

xpos = zeros(nImpPar, 4);
ypos = zeros(nImpPar, 4);
for k = 1:4
    xpos(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ypos(:, k) = hBar(k).YData;
end

errorbar(xpos, ypos, cell2mat(cOut(end - nImpPar + 1:end, 5:8)), '.k')
hold off
ylim([0, inf])
legend({'Tumour density at the end of treat.'...
    'Tumour density 3 months after the end of treat.'...
    '99% of initial tumour cells killed', 'Integral of tumour density'},...
    'location', 'northwest', 'fontsize', 26)
end
