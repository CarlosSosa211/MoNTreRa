clear all
close all

nfig = 0;
nOut = 8;
nTissues = 21;

for i = 1:nTissues
    path = ['../../Carlos/Results/Diff_Ang_432Sim_AllTissues/Tissue'...
        num2str(i)];
    diffEndTreatTumDens(:, :, i) = load([path, '/errEndTreatTumDens.res']);
    diff3MonTumDens(:, :, i)  = load([path, '/err3MonTumDens.res']);
    diffRecTumDens(:, :, i)  = load([path, '/errRecTumDens.res']);
    diffTumVol(:, :, i)  = load([path, '/errFinTumVol.res']);
    diffIntTumDens(:, :, i)  = load([path, '/errIntTumDens.res']);
    diffTimeTo95(:, :, i)  = load([path, '/errTimeTo95.res']);
    diffTimeTo99(:, :, i)  = load([path, '/errTimeTo99.res']);
    diffRecTime(:, :, i)  = load([path, '/errRecTime.res']);
end

par = diffEndTreatTumDens(:, 3:5, 1);

colTTum = 1;
colDThres = 2;
% colTArrest = 3;
colDose = 3;

tTTum = unique(par(:, colTTum));
tDThres = unique(par(:, colDThres));
% tTArrest = unique(par(:, colTArrest));
tDose = unique(par(:, colDose));

meanDiffEndTreatTumDens = mean(diffEndTreatTumDens, 3);
meanDiff3MonTumDens  = mean(diff3MonTumDens, 3);
meanDiffRecTumDens  = mean(diffRecTumDens, 3);
meanDiffTumVol  = mean(diffTumVol, 3);
meanDiffIntTumDens  = mean(diffIntTumDens, 3);
meanDiffTimeTo95  = mean(diffTimeTo95, 3);
meanDiffTimeTo99  = mean(diffTimeTo99, 3);
meanDiffRecTime  = mean(diffRecTime, 3);

%%

b = {'$endTreatTumDens$', '$3MonTumDens$', '$recTumDens$', '$tumVol$'...
    '$intTumDens$', '$timeTo95$', '$timeTo99$', '$recTime$'};

meanDiffRel = [num2cell(mean(meanDiffEndTreatTumDens(:, 2)))...
    num2cell(mean(meanDiff3MonTumDens(:, 2)))...
    num2cell(mean(meanDiffRecTumDens(:, 2)))...
    num2cell(mean(meanDiffTumVol(:, 2)))...
    num2cell(mean(meanDiffIntTumDens(:, 2)))...
    num2cell(mean(meanDiffTimeTo95(:, 2)))...
    num2cell(mean(meanDiffTimeTo99(:, 2)))...
    num2cell(mean(meanDiffRecTime(:, 2)));...
    num2cell(std(meanDiffEndTreatTumDens(:, 2)))...
    num2cell(std(meanDiff3MonTumDens(:, 2)))...
    num2cell(std(meanDiffRecTumDens(:, 2)))...
    num2cell(std(meanDiffTumVol(:, 2)))...
    num2cell(std(meanDiffIntTumDens(:, 2)))...
    num2cell(std(meanDiffTimeTo95(:, 2)))...
    num2cell(std(meanDiffTimeTo99(:, 2)))...
    num2cell(std(meanDiffRecTime(:, 2))); b];

nfig = nfig + 1;
figure(nfig);
hold on
bar(cell2mat(meanDiffRel(1, :)));
errorbar(1:nOut, cell2mat(meanDiffRel(1, :)),...
    cell2mat(meanDiffRel(2, :)), '.k')
hold off

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nOut)
set(ax,'XTickLabel', meanDiffRel(3, :), 'fontsize', 20);
xtickangle(45)
ax.YGrid = 'on';
ylim([0, inf])
title('All tissues - Relative differences', 'fontsize', 20)
