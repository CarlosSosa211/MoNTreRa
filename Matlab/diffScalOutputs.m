clear all
close all

nfig = 0;

% path = '../../Carlos/Results/Err_Ang_250Sim_Tissue8';
% path = '../../Carlos/Results/Err_AlphaBeta_432Sim';
% path = '../../Carlos/Results/Err_AlphaBeta_2592Sim';
% path = '../../Carlos/Results/Err_AngAlphaBeta_2592Sim';
% path = '../../Carlos/Results/Err_Ang_432Sim';
nTissue = 4;
path = ['../../Carlos/Results/Err_Ang_432Sim_AllTissues/Tissue', num2str(nTissue)];

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
    'intTumDens (5), timeTo95 (6), timeTo99 (7) or recTime (8)] or quit (0): ']);

switch selOut
    case 1
        path = [path, '/endTreatTumDens.res'];
        outputName = 'tumour density at the end of treat.';
    case 2
        path = [path, '/err3MonTumDens.res'];
        outputName = 'tumour density 3 months after the end of treat.';
    case 3
        path = [path, '/recTumDens.res'];
        outputName = 'tumour density at recurrence';
    case 4
        path = [path, '/finTumVol.res'];
        outputName = 'final tumour volume';
    case 5
        path = [path, '/intTumDens.res'];
        outputName = 'integral of tumour density';
    case 6
        path = [path, '/timeTo95.res'];
        outputName = 'time to kill 95% of tumour cells';
    case 7
        path = [path, '/timeTo99.res'];
        outputName = 'time to kill 99% of tumour cells';
    case 8
        path = [path, '/recTime.res'];
        outputName = 'recurrence time';
end

output = load(path);

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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in ', outputName])
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
title(['Tissue ', num2str(nTissue), ' - Histogram of ', outputName]);
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
title(['Tissue ', num2str(nTissue), ' - Histogram of difference in ', outputName])
