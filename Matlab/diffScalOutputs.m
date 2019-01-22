clear all
close all

nfig = 0;
nTissue = 4;
% path = ['../../Carlos/Results/Diff_Ang_432Sim_AllTissues/Tissue'...
%     num2str(nTissue)];
% path = '../../Carlos/Results/Diff_Ang_432Sim_Tissue4_15Out';
path = '../../Carlos/Results/Diff_Ang_432x5Sim_Tissue4';

colTTum = 1;
colDThres = 2;
% colTArrest = 3;
colDose = 3;

par = load([path, '/combPar.res']);

tTTum = unique(par(:, colTTum));
tDThres = unique(par(:, colDThres));
% tTArrest = unique(par(:, colTArrest));
tDose = unique(par(:, colDose));

fileNames = {'/endTreatTumDens.res', '/3MonTumDens.res'...
    '/finTumVol.res', '/intTumDens.res', '/killed50.res'...
    '/killed80.res''/killed90.res', '/killed95.res', '/killed99.res'...
    '/killed999.res', '/timeTo95.res', '/timeTo99.res', '/rec.res'...
    '/recTumDens.res', '/recTime.res'};
outputNames = {'tumour density at the end of treat.'...
    'tumour density 3 months after the end of treat.'...
    'final tumour volume', 'integral of tumour density'...
    '50% of tumour cells killed', '80% of tumour cells killed'...
    '90% of tumour cells killed', '95% of tumour cells killed'...
    '99% of tumour cells killed', '99.9% of tumour cells killed'...
    'time to kill 95% of tumour cells'...
    'time to kill 99% of tumour cells', 'recurrence'...
    'tumour density at recurrence', 'recurrence time'};

selOut = input(['Select an output [endTreaTumDens (1), 3MonTumDens (2) '...
    'finTumVol (3),\nintTumDens (4), 50%killed (5), 80%killed (6) '...
    '90%killed (7), 95%killed (8),\n99%killed (9), 99.9%killed (10) '...
    'timeTo95 (11), timeTo99 (12), rec(13),\nrecTumDens (14), '...
    'recTime (15) or all of them (-1)]: ']);

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

%%
noAngMeanTTum = [];
noAngStdTTum = [];
angMeanTTum = [];
angStdTTum = [];
for i = 1:length(tTTum)
    noAngMeanTTum = [noAngMeanTTum...
        mean(output(par(:, colTTum) == tTTum(i), 1))];
    noAngStdTTum = [noAngStdTTum...
        mean(output(par(:, colTTum) == tTTum(i), 3))];
    angMeanTTum = [angMeanTTum...
        mean(output(par(:, colTTum) == tTTum(i), 2))];
    angStdTTum = [angStdTTum...
        mean(output(par(:, colTTum) == tTTum(i), 4))];
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(tTTum, noAngMeanTTum, noAngStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
errorbar(tTTum, angMeanTTum, angStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold off
title(['Tissue ', num2str(nTissue), ' - '...
    char(outputNames(selOut))])
legend('No angiogenesis', 'Angiogenesis', 'location', 'northwest')
grid on
xlabel('TTum (h)')
ylabel(outputNames(selOut))

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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
grid on
xlabel('TTum (h)')
ylabel('Difference')

%%
noAngMeanDThres = [];
noAngStdDThres = [];
angMeanDThres = [];
angStdDThres = [];
for i = 1:length(tDThres)
    noAngMeanDThres = [noAngMeanDThres...
        mean(output(par(:, colDThres) == tDThres(i), 1))];
    noAngStdDThres = [noAngStdDThres...
        mean(output(par(:, colDThres) == tDThres(i), 3))];
    angMeanDThres = [angMeanDThres...
        mean(output(par(:, colDThres) == tDThres(i), 2))];
    angStdDThres = [angStdDThres...
        mean(output(par(:, colDThres) == tDThres(i), 4))];
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(tDThres, noAngMeanDThres, noAngStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
errorbar(tDThres, angMeanDThres, angStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold off
title(['Tissue ', num2str(nTissue), ' - '...
    char(outputNames(selOut))])
legend('No angiogenesis', 'Angiogenesis', 'location', 'northwest')
grid on
xlabel('Dose threshold (Gy)')
ylabel(outputNames(selOut))

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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
grid on
xlabel('TArrest (h)')
ylabel('Difference')

%%
noAngMeanDose = [];
noAngStdDose = [];
angMeanDose = [];
angStdDose = [];
for i = 1:length(tDose)
    noAngMeanDose = [noAngMeanDose...
        mean(output(par(:, colDose) == tDose(i), 1))];
    noAngStdDose = [noAngStdDose...
        mean(output(par(:, colDose) == tDose(i), 3))];
    angMeanDose = [angMeanDose...
        mean(output(par(:, colDose) == tDose(i), 2))];
    angStdDose = [angStdDose...
        mean(output(par(:, colDose) == tDose(i), 4))];
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(tDose, noAngMeanDose, noAngStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
errorbar(tDose, angMeanDose, angStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold off
title(['Tissue ', num2str(nTissue), ' - '...
    char(outputNames(selOut))])
legend('No angiogenesis', 'Angiogenesis', 'location', 'northwest')
grid on
xlabel('Dose (Gy)')
ylabel(outputNames(selOut))

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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
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
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
    char(outputNames(selOut))])
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
title(['Tissue ', num2str(nTissue), ' - Histogram of '...
    char(outputNames(selOut))]);
legend('No angiogenesis', 'Angiogenesis')

nfig = nfig + 1;
figure(nfig)
hDiff = histogram(outputDiff, 'binmethod', 'scott');
hDiff.BinEdges = hDiff.BinEdges + hDiff.BinWidth/2;
title(['Tissue ', num2str(nTissue), ' - Histogram of difference in '...
    char(outputNames(selOut))])


%%
b = {'$endTreatTumDens$'; '$3MonTumDens$'; '$tumVol$'; '$intTumDens$';...
    '$50\%killed$'; '$80\%killed$'; '$90\%killed$'; '$95\%killed$';...
    '$99\%killed$'; '$99.9%killed$'; '$timeTo95$'; '$timeTo99$';...
    '$rec$', '$recTime$'; '$recTumDens$';};
nOut = 15;

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
title(['Tissue ', num2str(nTissue), ' - Relative differences'],...
    'fontsize', 20)
