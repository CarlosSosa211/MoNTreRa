clear all
close all

%%
% This block initialises nfig and defines
% - path: the path of the ouput directory,
% - nTissues: the number of tissues,
% - withoutN: the name of the X_0_X simulations
% - withN: the name of the X_1_X simulations
% - fileNames: the names of the ouput files,
% - ouputNames: the names of the outputs,
% - outputCol: the column number of the output in the corresponding files.

nfig = 0;
path = '../../Carlos/Results/Diff/Ang_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/Res_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/Arrest_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/HypNec_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/Oxy_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/AngRes_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/ArrestNoAngNoRes_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/OxyNoHypNec_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/OxyNoAngNoHypNec_Dose_5Val_5Rep/Tissue';

initVascDens = load('../HistSpec/initVascDens.dat');

nTissues = 21;
[initVascDens, tissuesVascDens] = sort(initVascDens);
nOut = 15;

withoutN = 'No angiogenesis';
withN = 'Angiogenesis';
% withoutN = 'No healthy cell division';
% withN = 'Healthy cell division';
% withoutN = 'No arrest';
% withN = 'Arrest';
% withoutN = 'No hypoxic necrosis';
% withN = 'Hypoxic necrosis';
% withoutN = 'No oxyegenation (no hypoxic necrosis)';
% withN = 'Oxygenation';
% withoutN = 'No angiogenesis and no healthy cell division';
% withN = 'Angiogenesis and healthy cell division';
% withoutN = 'No arrest, no angiogenesis, no healthy cell division';
% withN = 'Arrest, no angiogenesis, no healthy cell division';
% withoutN = 'No oxyegenation (no hypoxic necrosis)';
% withN = 'Oxygenation (no hypoxic necrosis)';
% withoutN = 'No oxyegenation (no hypoxic necrosis)';
% withN = 'Oxygenation (no angiogenesis, no hypoxic necrosis)';

fileNames = {'/endTreatTumDens.res', '/3MonTumDens.res'...
    '/finTumVol.res', '/intTumDens.res', '/killed50.res'...
    '/killed80.res', '/killed90.res', '/killed95.res', '/killed99.res'...
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
% This block asks the user to select the tissues and the output to be
% studied. Then, it reads the corresponding files.
% The following intermidiate variables are considered:
% - par: a matrix containing the combinations of parameters simulated. Each
% row corresponds to a simulation and each column, to a parameter,
% - output: a matrix containing the values of the scalar output(s) in
% question. For the single tissue - single output case
% (1 <= nTissue <= nTissues and 1 <= selOut <= nOut), each row
% corresponds to a combination of parameters. The first column corresponds
% to the mean values for nRep repetitions of X_0_X simulations; the second
% one, to the mean values of X_1_X simulations; the third one, two the std
% values of X_0_X simulations and the fourth one, to the std values of
% X_1_X simulations. For the single tissue - multiple outputs case
% (1 <= nTissue <= nTissues and selOut = 0), each layer correponds to an
% output. For the multiple tissues - single output case, (nTissue = 0 and
% 1 <= selOut <= nOut) each layer corresponds to a tissue. Finally, for the
% mean over the tissues - single output case (nTissue = -1 and
% 1 <= selOut <= nOut), the matrix has a single layer with the mean values.

nTissue = input(['Select one tissue (from 1 to ', num2str(nTissues)...
    ') or all of them (0) or a mean over them (-1): ']);

if(nTissue >= 1 && nTissue <= nTissues)
    pathTissue = [path, num2str(nTissue)];
    par = load([pathTissue, '/combPar.res']);
    
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
    
    if(selOut >= 1 && selOut <= nOut)
        output = load([pathTissue, char(fileNames(selOut))]);
        
    elseif(selOut == -1)
        for i = 1:length(fileNames)
            output(:, :, i) = load([pathTissue, char(fileNames(i))]);
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
% This block calculates the difference and the absolute difference of the
% mean values of X_0_X and X_1_X simulations.
outputDiff = output(:, 1)' - output(:, 2)';
outputAbsDiff = abs(outputDiff);

%%
selPar = 4;
switch selPar
    case 1
        parName = 'TTum (h)';
        parT = tTum;
        parCol = colTTum;
    case 2
        parName = 'Dose threshold (Gy)';
        parT = tDose;
        parCol = colDose;
    case 3
        parName = 'TArrest (h)';
        parT = tTArrest;
        parCol = colTArrest;
    case 4
        parName = 'Dose (Gy)';
        parT = tDose;
        parCol = colDose;
end

%%
% This block plots the mean and std of the selected output for simulations
% X_0_X and X_1_1 as a function of the studied values of the parameter par

withoutMean = zeros(1, length(parT));
withoutStd = zeros(1, length(parT));
withMean = zeros(1, length(parT));
withStd = zeros(1, length(parT));

for i = 1:length(parT)
    withoutMean(i) = mean(output(par(:, parCol) == parT(i), 1));
    withoutStd(i) = mean(output(par(:, parCol) == parT(i), 3));
    withMean(i) = mean(output(par(:, parCol) == parT(i), 2));
    withStd(i) = mean(output(par(:, parCol) == parT(i), 4));
end

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(parT, withoutMean, withoutStd, 'sb', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
errorbar(parT, withMean, withStd, 'sr', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
hold off
if(nTissue == - 1)
    title(['All tissues - ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - '...
        char(outputNames(selOut))])
end
legend(withoutN, withN, 'location', 'northeast')
grid on
xlabel(parName)
ylabel(outputNames(selOut))

%%
% This block plots the mean and std of the abs. difference of the selected
% output for simulations X_0_X and X_1_1 as a function of the studied
% values of the parameter par

diffMean = zeros(1, length(parT));
diffStd = zeros(1, length(parT));

for i = 1:length(parT)
    diffMean(i) = mean(outputAbsDiff(par(:, parCol) == parT(i)));
    diffStd(i) = std(outputAbsDiff(par(:, parCol) == parT(i)));
end

nfig = nfig + 1;
figure(nfig)
errorbar(parT, diffMean, diffStd, '-s', 'MarkerSize', 10,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
grid on
xlabel(parName)
ylabel('Difference')

%%
selPar1 = 4;
switch selPar1
    case 1
        parName1 = 'TTum (h)';
        parT1 = tTum;
        parCol1 = colTTum;
    case 2
        parName1 = 'Dose threshold (Gy)';
        parT1 = tDose;
        parCol1 = colDose;
    case 3
        parName1 = 'TArrest (h)';
        parT1 = tTArrest;
        parCol1 = colTArrest;
    case 4
        parName1 = 'Dose (Gy)';
        parT1 = tDose;
        parCol1 = colDose;
end

selPar2 = 4;
switch selPar2
    case 1
        parName2 = 'TTum (h)';
        parT2 = tTum;
        parCol2 = colTTum;
    case 2
        parName2 = 'Dose threshold (Gy)';
        parT2 = tDose;
        parCol2 = colDose;
    case 3
        parName2 = 'TArrest (h)';
        parT2 = tTArrest;
        parCol2 = colTArrest;
    case 4
        parName2 = 'Dose (Gy)';
        parT2 = tDose;
        parCol2 = colDose;
end

%%
% This block plots the mean of the abs. difference of the selected output
% for simulations X_0_X and X_1_1 as a function of the studied values of
% parameters par1 and par2

diff = zeros(length(parT1), length(parT2));
for i = 1:length(parT1)
    for j = 1:length(parT2)
        diff(i, j) =...
            mean(outputAbsDiff(par(:, parCol1) == parT1(i) &...
            par(:, parCol2) == parT2(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diff, 'CDataMapping','scaled')
colorbar
if(nTissue == -1)
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))])
else
    title(['Tissue ', num2str(nTissue), ' - Abs. diff. in '...
        char(outputNames(selOut))])
end
xlabel(parName2)
ylabel(parName1)
xticks(1:length(parT2))
xticklabels(parT2)
yticks(1:length(parT1))
yticklabels(parT1)

%%
% This block plots the histogram of the values of the selected
% output for simulations X_0_X and X_1_1 and the histogram of the
% differences between them.

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
% This block plots the mean and std values of the abs. differences between
% simulations X_0_X and X_1_1 for all the outputs.

b = {'$endTreatTumDens$'; '$3MonTumDens$'; '$tumVol$'; '$intTumDens$';...
    '$50\%killed$'; '$80\%killed$'; '$90\%killed$'; '$95\%killed$';...
    '$99\%killed$'; '$99.9%killed$'; '$timeTo95$'; '$timeTo99$';...
    '$rec$'; '$recTime$'; '$recTumDens$';};

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
% This block plots for, all the values of dose, the mean and std of the
% selected output for simulations X_0_X and X_1_1 as a function of the
% tissue

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
%     errorbar(1:nTissues, withoutMeanDose(i, tissuesVascDens),...
%         withoutStdDose(i, tissuesVascDens), 'sb', 'MarkerSize', 10,...
%         'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
%     errorbar(1:nTissues, withMeanDose(i, tissuesVascDens),...
%         withStdDose(i, tissuesVascDens), 'sr', 'MarkerSize', 10,...
%         'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
    errorbar(initVascDens, withoutMeanDose(i, tissuesVascDens),...
        withoutStdDose(i, tissuesVascDens), 'sb', 'MarkerSize', 10,...
        'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b')
    errorbar(initVascDens, withMeanDose(i, tissuesVascDens),...
        withStdDose(i, tissuesVascDens), 'sr', 'MarkerSize', 10,...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
    hold off
    title(['All tissues - ', char(outputNames(selOut)), ' - '...
        num2str(tDose(i)), ' Gy'])
    legend(withoutN, withN, 'location', 'northeast')
    grid on
    ylim([0, inf])
    %     xticks(1:nTissues)
    %     xticklabels(num2cell(tissuesVascDens))
    %     xticklabels(num2cell(initVascDens))
    %     xlabel('Tissue')
    xlabel('Initial vascular density')
    ylabel(outputNames(selOut))
end

%%
% This block plots, for all the values of dose, the mean and std of the
% abs. differences of the selected output for simulations X_0_X and X_1_1
% as a function of the tissue

absDiffDose = abs(output(:, 1, :) - output(:, 2, :));
absDiffDose = permute(absDiffDose, [1, 3, 2]);

for i = 1:length(tDose)
    nfig = nfig + 1;
    figure(nfig)
    hold on
    plot(initVascDens, absDiffDose(i, tissuesVascDens), 'sk', 'MarkerSize',...
        10, 'MarkerFaceColor', 'k')
    title(['All tissues - Abs. diff. in ', char(outputNames(selOut))...
        ' - ', num2str(tDose(i)), ' Gy'])
    grid on
    ylim([0, inf])
    %     xticks(1:nTissues)
    %     xticklabels(num2cell(tissuesVascDens))
    %     xticklabels(num2cell(initVascDens))
    %     xlabel('Tissue')
    xlabel('Initial vascular density')
    ylabel(outputNames(selOut))
end

%%
% This block plots, for all the values of dose, the mean and std of the
% rel. differences of the selected output for simulations X_0_X and X_1_1
% as a function of the tissue
relDiffDose = (abs(output(:, 1, :) - output(:, 2, :))) ./ output(:, 2, :);
relDiffDose = permute(relDiffDose, [1, 3, 2]);

for i = 1:length(tDose)
    nfig = nfig + 1;
    figure(nfig)
    hold on
    plot(initVascDens, relDiffDose(i, tissuesVascDens), 'sk', 'MarkerSize',...
        10, 'MarkerFaceColor', 'k')
    title(['All tissues - Rel. diff. in ', char(outputNames(selOut))...
        ' - ', num2str(tDose(i)), ' Gy'])
    grid on
    ylim([0, inf])
    %     xticks(1:nTissues)
    %     xticklabels(num2cell(tissuesVascDens))
    %     xticklabels(num2cell(initVascDens))
    %     xlabel('Tissue')
    xlabel('Initial vascular density')
    ylabel(outputNames(selOut))
end