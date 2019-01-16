clear all
close all

nfig = 0;

nTissue = 4;
path = ['../../Carlos/Results/Diff_Ang_432Sim_AllTissues/Tissue', num2str(nTissue)];

colTTum = 1;
colDThres = 2;
% colTArrest = 3;
colDose = 3;

par = load([path, '/combPar.res']);

tTTum = unique(par(:, colTTum));
tDThres = unique(par(:, colDThres));
% tTArrest = unique(par(:, colTArrest));
tDose = unique(par(:, colDose));

selOut = input(['Select an output [endTreaTumDens (1), 3MonTumDens (2), recTumDens (3), tumVol (4),\n'...
    'intTumDens (5), timeTo95 (6), timeTo99 (7), recTime (8) or all of them (-1)]: ']);

fileNames = {'/endTreatTumDens.res', '/3MonTumDens.res', '/recTumDens.res',...
    '/finTumVol.res', '/intTumDens.res', '/timeTo95.res', '/timeTo99.res'...
    '/recTime.res'};
outputNames = {'tumour density at the end of treat.', 'tumour density 3 months after the end of treat.'...
    'tumour density at recurrence', 'final tumour volume', 'integral of tumour density'...
    'time to kill 95% of tumour cells', 'time to kill 99% of tumour cells', 'recurrence time'};

if(selOut > 0)
    output = load([path, char(fileNames(selOut))]);
    
elseif(selOut == -1)
    for i = 1:length(fileNames)
        output(:, :, i) = load([path, char(fileNames(i))]);
    end
end

%%
outputDiff = output(:, 1)' - output(:, 2)';
outputAbsDiff = abs(outputDiff);
outputRelDiff = outputAbsDiff ./ max(output(:, 1)', output(:, 2)');

%%
diffMeanTTum = [];
diffStdTTum = [];
for i = 1:length(tTTum)
    diffMeanTTum = [diffMeanTTum...
        mean(outputAbsDiff(par(:, colTTum) == tTTum(i)))];
    diffStdTTum = [diffStdTTum...
        std(outputAbsDiff(par(:, colTTum) == tTTum(i)))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, diffMeanTTum, diffStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
grid on
xlabel('TTum (h)')
ylabel('Difference')

%%
diffMeanDThres = [];
diffStdDThres = [];
for i = 1:length(tDThres)
    diffMeanDThres = [diffMeanDThres...
        mean(outputAbsDiff(par(:, colDThres) == tDThres(i)))];
    diffStdDThres = [diffStdDThres...
        std(outputAbsDiff(par(:, colDThres) == tDThres(i)))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDThres, diffMeanDThres, diffStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
grid on
xlabel('Dose threshold(Gy)')
ylabel('Difference')

%%
diffMeanTArrest = [];
diffStdTArrest = [];
for i = 1:length(tTArrest)
    diffMeanTArrest = [diffMeanTArrest...
        mean(outputAbsDiff(par(:, colTArrest) == tTArrest(i)))];
    diffStdTArrest = [diffStdTArrest...
        std(outputAbsDiff(par(:, colTArrest) == tTArrest(i)))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTArrest, diffMeanTArrest, diffStdTArrest,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
grid on
xlabel('TArrest (h)')
ylabel('Difference')

%%
diffMeanDose = [];
diffStdDose = [];
for i = 1:length(tDose)
    diffMeanDose = [diffMeanDose...
        mean(outputAbsDiff(par(:, colDose) == tDose(i)))];
    diffStdDose = [diffStdDose...
        std(outputAbsDiff(par(:, colDose) == tDose(i)))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, diffMeanDose, diffStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
grid on
xlabel('Dose (Gy)')
ylabel('Difference')

%%
diffTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        diffTTumDose(i, j) =...
            mean(outputAbsDiff(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        diffDThresDose(i, j) =...
            mean(outputAbsDiff(par(:, colDThres) == tDThres(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffDThresDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
xlabel('Dose (Gy)')
ylabel('Dthres (Gy)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tDThres))
yticklabels(tDThres)

%%
diffTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        diffTTumDThres(i, j) =...
            mean(outputAbsDiff(par(:, colTTum) == tTTum(i) &...
            par(:, colDThres) == tDThres(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffTTumDThres, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
xlabel('DThres (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        diffTTumDThres(i, j) =...
            mean(outputAbsDiff(par(:, colTTum) == tTTum(i) &...
            par(:, colDThres) == tDThres(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffTTumDThres, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', char(outputNames(selOut))])
xlabel('DThres (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
nfig = nfig + 1;
figure(nfig)
hold on
histogram(output(:, 1), 25);
histogram(output(:, 2), 25);
hold off
title(['Tissue ', num2str(nTissue), ' - Histogram of ', char(outputNames(selOut))]);
legend('No angiogenesis', 'Angiogenesis')

% binSize = 1/16;
% bound = 2000;
% lBound = - bound - 0.5 * binSize;
% rBound = bound + 0.5 * binSize;
% edges = [min([outputDiff, lBound]), lBound:binSize:rBound, max([outputDiff, rBound])];
nfig = nfig + 1;
figure(nfig)
% histogram(outputDiff, edges);
histogram(outputDiff, 25);
title(['Tissue ', num2str(nTissue), ' - Histogram of difference in ', char(outputNames(selOut))])


%%

%%TO-DO verificar que las dimensiones de outputDiff 


b = {'$endTreatTumDens$'; '$3MonTumDens$'; '$recTumDens$'; '$tumVol$';...
    '$intTumDens$'; '$timeTo95$'; '$timeTo99$'; '$recTime$'};
nOut = 8;

outputDiff = output(:, 1, :) - output(:, 2, :);
outputDiff = permute(outputDiff, [3, 1, 2]);
maxOutput = max(output(:, 1, :), output(:, 2, :));
maxOutput = permute(maxOutput, [3, 1, 2]);
outputAbsDiff = abs(outputDiff);
outputRelDiff = outputAbsDiff ./ maxOutput;
outputRelDiff(isnan(outputRelDiff)) = 0;
meanOutputRelDiff = mean(outputRelDiff, 2);
stdOutputRelDiff = std(outputRelDiff, 0, 2);

diffRel = [num2cell(meanOutputRelDiff), num2cell(stdOutputRelDiff), b];

nfig = nfig + 1;
figure(nfig);
hold on
bar(cell2mat(diffRel(:, 1)));
errorbar(1:nOut, cell2mat(diffRel(:, 1)), cell2mat(diffRel(:, 2)), '.k')
hold off

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nOut)
set(ax,'XTickLabel', diffRel(:, 3), 'fontsize', 20);
xtickangle(45)
ax.YGrid = 'on';
ylim([0, inf])
title(['Tissue ', num2str(nTissue), ' - Relative differences'], 'fontsize', 20)
