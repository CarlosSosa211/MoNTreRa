clear all
close all

nfig = 0;
% path = ['../../Carlos/Results/Diff_Ang_432Sim_AllTissues/Tissue'...
%     num2str(nTissue)];
% path = '../../Carlos/Results/Diff_Ang_432Sim_Tissue4_15Out';
% path = ['../../Carlos/Results/Diff_Ang_432x5Sim_Tissue', num2str(nTissue)];
% path = '../../Carlos/Results/Diff_Ang_10x3Sim_AllTissues/Tissue';
% path = '../../Carlos/Results/Diff_Ang_Dose_5Values_5Rep_AllTissues/Tissue';
path = '../../Carlos/Results/Diff_Res_Dose_5Values_5Rep_AllTissues/Tissue';
% path = '../../Carlos/Results/Diff_Ang_10x3Sim_AllTissues_noHypNec/Tissue';
withoutN = 'No healthy cell division';
withN = 'Healthy cell division';
nTissues = 21;
% withoutN = 'No angiogenesis';
% withN = 'Angiogenesis';

fileNames = {'/endTreatTumDens.res', '/3MonTumDens.res'...
    '/finTumVol.res', '/intTumDens.res', '/killed50.res'...
    '/killed80.res''/killed90.res', '/killed95.res', '/killed99.res'...
    '/killed999.res', '/timeTo95.res', '/timeTo99.res', '/rec.res'...
    '/recTumDens.res', '/recTime.res'};

outputNames = {'tumour density at the end of treat.'...
    'tumour density 3 months after the beginning of treat.'...
    'final tumour volume', 'integral of tumour density'...
    '50% of tumour cells killed', '80% of tumour cells killed'...
    '90% of tumour cells killed', '95% of tumour cells killed'...
    '99% of tumour cells killed', '99.9% of tumour cells killed'...
    'time to kill 95% of tumour cells'...
    'time to kill 99% of tumour cells', 'recurrence'...
    'tumour density at recurrence', 'recurrence time'};

%%
nTissue = input(['Select one tissue (from 1 to ', num2str(nTissues)...
    ') or all of them (0) or a mean over them (-1): ']);

if(nTissue >= 1 && nTissue <= nTissues)
    pathTissue = [path, num2str(nTissue)];
    par = load([path, '/combPar.res']);
    
    % colTTum = 1;
    % colDThres = 2;
    % colTArrest = 3;
    colDose = 1;
    tDose = unique(par(:, colDose));
    % tTTum = unique(par(:, colTTum));
    % tDThres = unique(par(:, colDThres));
    % tTArrest = unique(par(:, colTArrest));
    
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\nintTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n'...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\nrecTumDens (14), recTime (15) '...
        'or all of them (-1)]: ']);
    
    if(selOut >= 1 && selOut <= 15)
        output = load([path, char(fileNames(selOut))]);
        
    elseif(selOut == -1)
        for i = 1:length(fileNames)
            output(:, :, i) = load([path, char(fileNames(i))]);
        end
    end
    
elseif(nTissue == 0)
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\nintTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n'...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\nrecTumDens (14), recTime (15) '...
        'or all of them (-1)]: ']);
    
    for i = 1:nTissues
        pathTissue = [path, num2str(i)];
        par(:, :, i) = load([pathTissue, '/combPar.res']);
        output(:, :, i) = load([pathTissue, char(fileNames(selOut))]);
    end
    colDose = 1;
    tDose = unique(par(:, colDose, 1));
    
elseif(nTissue == -1)
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\nintTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n'...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\nrecTumDens (14), recTime (15) '...
        'or all of them (-1)]: ']);
    
    for i = 1:nTissues
        pathTissue = [path, num2str(i)];
        par(:, :, i) = load([pathTissue, '/combPar.res']);
        output(:, :, i) = load([pathTissue, char(fileNames(selOut))]);
    end
    % colTTum = 1;
    % colDThres = 2;
    % colTArrest = 3;
    colDose = 1;
    tDose = unique(par(:, colDose, 1));
    % tTTum = unique(par(:, colTTum));
    % tDThres = unique(par(:, colDThres));
    % tTArrest = unique(par(:, colTArrest));
    output = mean(output, 3);
end

%%
outputDiff = output(:, 1)' - output(:, 2)';
outputAbsDiff = abs(outputDiff);

%%
withoutMeanTTum = [];
withoutStdTTum = [];
withMeanTTum = [];
withStdTTum = [];
for i = 1:length(tTTum)
    withoutMeanTTum = [withoutMeanTTum...
        mean(output(par(:, colTTum) == tTTum(i), 1))];
    withoutStdTTum = [withoutStdTTum...
        mean(output(par(:, colTTum) == tTTum(i), 3))];
    withMeanTTum = [withMeanTTum...
        mean(output(par(:, colTTum) == tTTum(i), 2))];
    withStdTTum = [withStdTTum...
        mean(output(par(:, colTTum) == tTTum(i), 4))];
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(tTTum, withoutMeanTTum, withoutStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
errorbar(tTTum, withMeanTTum, withStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold off
if(nTissue == - 1)
    title(['All tissues - ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - '...
        char(outputNames(selOut))])
end
legend(withoutN, withN, 'location', 'northwest')
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
grid on
xlabel('TTum (h)')
ylabel('Difference')

%%
withoutMeanDThres = [];
withoutStdDThres = [];
withMeanDThres = [];
withStdDThres = [];
for i = 1:length(tDThres)
    withoutMeanDThres = [withoutMeanDThres...
        mean(output(par(:, colDThres) == tDThres(i), 1))];
    withoutStdDThres = [withoutStdDThres...
        mean(output(par(:, colDThres) == tDThres(i), 3))];
    withMeanDThres = [withMeanDThres...
        mean(output(par(:, colDThres) == tDThres(i), 2))];
    withStdDThres = [withStdDThres...
        mean(output(par(:, colDThres) == tDThres(i), 4))];
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(tDThres, withoutMeanDThres, withoutStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
errorbar(tDThres, withMeanDThres, withStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
hold off
if(nTissue == -1)
    title(['All tissues - ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - '...
        char(outputNames(selOut))])
end
legend(withoutN, withN, 'location', 'northwest')
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
grid on
xlabel('TArrest (h)')
ylabel('Difference')

%%
withoutMeanDose = [];
withoutStdDose = [];
withMeanDose = [];
withStdDose = [];
for i = 1:length(tDose)
    withoutMeanDose = [withoutMeanDose...
        mean(output(par(:, colDose) == tDose(i), 1))];
    withoutStdDose = [withoutStdDose...
        mean(output(par(:, colDose) == tDose(i), 3))];
    withMeanDose = [withMeanDose...
        mean(output(par(:, colDose) == tDose(i), 2))];
    withStdDose = [withStdDose...
        mean(output(par(:, colDose) == tDose(i), 4))];
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(tDose, withoutMeanDose, withoutStdDose,...
    'sb', 'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
errorbar(tDose, withMeanDose, withStdDose,...
    'sr', 'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
hold off
if(nTissue == -1)
    title(['All tissues  - ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - '...
        char(outputNames(selOut))])
end
legend(withoutN, withN, 'location', 'northwest')
grid on
ylim([0, inf])
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
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
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
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
if(nTissue == -1)
    title(['All tissues  - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
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
if(nTissue == -1)
    title(['All tissues - Histogram of ', char(outputNames(selOut))]);
else
    title(['Tissue ', num2str(nTissue), ' - Histogram of '...
        char(outputNames(selOut))]);
end
legend(withoutN, withN, 'location', 'northwest')

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
xtickwithle(45)
ax.YGrid = 'on';
ylim([0, inf])
if(nTissue == -1)
    title('All tissues - Relative differences', 'fontsize', 20)
else
    title(['Tissue ', num2str(nTissue), ' - Relative differences'],...
        'fontsize', 20)
end

%%
withoutMeanDose = zeros(length(tDose), nTissues);
withoutStdDose = zeros(length(tDose), nTissues);
withMeanDose = zeros(length(tDose), nTissues);
withStdDose = zeros(length(tDose), nTissues);

for i = 1:length(tDose)
    for j = 1:nTissues
        withoutMeanDose(i, j) = output(par(:, colDose) == tDose(i), 1, j);
        withoutStdDose(i, j) = output(par(:, colDose) == tDose(i), 3, j);
        withMeanDose(i, j) = output(par(:, colDose) == tDose(i), 2, j);
        withStdDose(i, j) = output(par(:, colDose) == tDose(i), 4, j);
    end
end

for i = 1:length(tDose)
    nfig = nfig + 1;
    figure(nfig)
    hold on
    errorbar(1:nTissues, withoutMeanDose(i, :), withoutStdDose(i, :),...
        'sb', 'MarkerSize', 10, 'MarkerEdgeColor', 'b',...
        'MarkerFaceColor', 'b')
    errorbar(1:nTissues, withMeanDose(i, :), withStdDose(i, :),...
        'sr', 'MarkerSize', 10, 'MarkerEdgeColor', 'r',...
        'MarkerFaceColor', 'r')
    hold off
    title(['All tissues  - ', char(outputNames(selOut)), '- '...
        num2str(tDose(i)), ' Gy'])
    legend(withoutN, withN, 'location', 'northwest')
    grid on
    ylim([0, inf])
    xticks(1:nTissues)
    xlabel('Tissue')
    ylabel(outputNames(selOut))
end
