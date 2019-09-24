clear all
close all

%%
% This block initialises nfig and defines
% - path: the path of the diff ouput directory,
% - pathSim: the paths to the considered output directories,
% - nTissues: the number of tissues,
% - withoutN: the name of the X_0_X simulations
% - withN: the name of the X_1_X simulations
% - fileNames: the names of the ouput files,
% - ouputNames: the names of the outputs,
% - outputCol: the column number of the output in the corresponding files.

nfig = 0;
path = '../../Carlos/Results/Diff/';
% pathSim = {'Ang_Dose_5Val_5Rep/', 'Ang_Dose_5Val_5Rep/'...
%     'Res_Dose_5Val_5Rep/', 'Arrest_Dose_5Val_5Rep/'...
%     'HypNec_Dose_5Val_5Rep/'};
% tSim = [1, 0, 0, 0, 0];
% pathSim = {'AngRes_Dose_5Val_5Rep/', 'AngRes_Dose_5Val_5Rep/'...
%     'AngArrest_Dose_5Val_5Rep/', 'ResArrest_Dose_5Val_5Rep/'};
% tSim = [1, 0, 0, 0];
pathSim = {'Ang_Dose_5Val_5Rep/', 'ArrestNoAngNoRes_Dose_5Val_5Rep/'};
tSim = [1, 0];
% pathSim = {'Ang_Dose_5Val_5Rep/', 'OxyNoAngNoHypNec_Dose_5Val_5Rep/'};
% tSim = [1, 1];
% pathSim = {'TimeConstantOxy_Dose_5Val_5Rep/'...
%     'TimeConstantOxy_Dose_5Val_5Rep/'};
% tSim = [0, 1];
nSim = length(tSim);

initVascDens = load('../HistSpec/initVascDens.dat');

nTissues = 21;
[initVascDens, tissuesVascDens] = sort(initVascDens);
nOut = 15;

% simN = {'All processes', 'No angiogenesis', 'No healthy cell division'...
%     'No cycle arrest', 'No hypoxic necrosis'};
% simN = {'All processes', 'No angiogenesis, no healthy cell division'...
%     'No angiogenesis, no cycle arrest'...
%     'No healthy cell division, no cycle arrest'};
simN = {'All processes'...
    'No angiogenesis, no healthy cell division, no cycle arrest'};
% simN = {'All processes', 'No angiogenesis, no hypoxic necrosis'};
% simN = {'All processes', 'Time constant oxygenation'};

fileNames = {'/endTreatTumDens.res', '/3MonTumDens.res'...
    '/finTumVol.res', '/intTumDens.res', '/killed50.res'...
    '/killed80.res', '/killed90.res', '/killed95.res', '/killed99.res'...
    '/killed999.res', '/timeTo95.res', '/timeTo99.res', '/rec.res'...
    '/recTumDens.res', '/recTime.res'};

outputNames = {'Tumour density at the end of treat.'...
    'Tumour density 3 months after the beginning of treat.'...
    'Final tumour volume', 'Integral of tumour density'...
    '50% of tumour cells killed', '80% of tumour cells killed'...
    '90% of tumour cells killed', '95% of tumour cells killed'...
    '99% of tumour cells killed', '99.9% of tumour cells killed'...
    'Time to kill 95% of tumour cells'...
    'Time to kill 99% of tumour cells', 'Recurrence'...
    'Tumour density at recurrence', 'Recurrence time'};

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
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\nintTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n'...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\nrecTumDens (14), recTime (15) '...
        'or all of them (-1)]: ']);
    
    if(selOut >= 1 && selOut <= nOut)
        for j = 1:nSim
            pathTissue = [path, char(pathSim(j)), '/Tissue'...
                num2str(nTissue)];
            par = load([pathTissue, '/combPar.res']);
            colDose = 1;
            output(:, :, j) = load([pathTissue, char(fileNames(selOut))]);
        end
        
    elseif(selOut == -1)
        for j = 1:nSim
            pathTissue = [path, char(pathSim(j)), '/Tissue'...
                num2str(nTissue)];
            par = load([pathTissue, '/combPar.res']);
            colDose = 1;
            for i = 1:length(fileNames)
                output(:, :, i, j) = load([pathTissue,...
                    char(fileNames(i))]);
            end
        end
    end
    
elseif(nTissue == 0)
    selOut = input(['Select an output [endTreaTumDens (1), '...
        '3MonTumDens (2), finTumVol (3),\nintTumDens (4), '...
        '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n'...
        '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
        'timeTo99 (12), rec(13),\nrecTumDens (14), recTime (15) '...
        'or all of them (-1)]: ']);
    
    for j = 1:nSim
        for i = 1:nTissues
            pathTissue = [path, char(pathSim(j)), '/Tissue' num2str(i)];
            par(:, :, i, j) = load([pathTissue, '/combPar.res']);
            output(:, :, i, j) = load([pathTissue,...
                char(fileNames(selOut))]);
        end
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
    
    for j = 1:nSim
        for i = 1:nTissues
            pathTissue = [path, char(pathSim(j)), '/Tissue',...
                num2str(i)];
            par(:, :, i, j) = load([pathTissue, '/combPar.res']);
            output(:, :, i, j) = load([pathTissue,...
                char(fileNames(selOut))]);
        end
    end
    colDose = 1;
    tDose = unique(par(:, colDose, 1));
    output = mean(output, 3);
    output = permute(output, [1, 2, 4, 3]);
end

%%
% This block plots, for all the values of dose, the selected output as a
% function of the tissue, sorted by initial vascular density

meanDose = zeros(length(tDose), nTissues, nSim);
stdDose = zeros(length(tDose), nTissues, nSim);

for i = 1:length(tDose)
    for j = 1:nTissues
        for k = 1:nSim
            meanDose(i, j, k) = output(par(:, colDose) == tDose(i),...
                1 + tSim(k), j, k);
            stdDose(i, j, k) = output(par(:, colDose) == tDose(i),...
                3 + tSim(k), j, k);
        end
        
    end
end

for i = 1:length(tDose)
    nfig = nfig + 1;
    figure(nfig)
    
    hold on
    hBar = bar(permute(meanDose(i, tissuesVascDens, :), [2, 3, 1]));
    xticks(1:nTissues)
    xticklabels(num2cell(tissuesVascDens))
    ax = gca;
    ax.FontSize = 22;
    title(['All tissues - ', char(outputNames(selOut)), ' - '...
        num2str(tDose(i)), ' Gy'], 'FontSize', 22)
    
    xpos = zeros(nTissues, nSim);
    ypos = zeros(nTissues, nSim);
    for k = 1:nSim
        xpos(:, k) = hBar(1).XData + hBar(k).XOffset;
        ypos(:, k) = hBar(k).YData;
    end
    errorbar(xpos, ypos, permute(stdDose(i, tissuesVascDens, :),...
        [2, 3, 1]), '.k')
    hold off
    
    ylim([0, 1.1 * max(max(ypos))])
    grid on
    legend(simN, 'location', 'northeast', 'FontSize', 22)
    xlabel('Tissue', 'FontSize', 22)
    ylabel(char(outputNames(selOut)), 'FontSize', 22)
end

%%
% This block plots the selected output as a function of the dose

nfig = nfig + 1;
figure(nfig)

meanOut = zeros(length(tDose), nSim);
stdOut = zeros(length(tDose), nSim);

for i = 1:length(tDose)
    for k = 1:nSim
        meanOut(i, k) = output(par(:, colDose) == tDose(i),...
            1 + tSim(k), k);
        stdOut(i, k) = output(par(:, colDose) == tDose(i),...
            3 + tSim(k), k);
    end
end

hold on
hBar = bar(meanOut);
xticks(tDose')
ax = gca;
ax.FontSize = 22;
title(['All tissues - ', char(outputNames(selOut))], 'FontSize', 22)
xpos = zeros(length(tDose), nSim);
ypos = zeros(length(tDose), nSim);
for k = 1:nSim
    xpos(:, k) = hBar(1).XData + hBar(k).XOffset;
    ypos(:, k) = hBar(k).YData;
end

errorbar(xpos, ypos, stdOut, '.k')
hold off

ylim([0, inf])
grid on
legend(simN, 'location', 'northeast', 'FontSize', 22)
xlabel('Dose (Gy)', 'FontSize', 22)
ylabel(char(outputNames(selOut)), 'FontSize', 22)

%%
% This block plots, for all the values of dose, the selected output as a
% function of the initial vascular density

meanDose = zeros(length(tDose), nTissues, nSim);
stdDose = zeros(length(tDose), nTissues, nSim);

for i = 1:length(tDose)
    for j = 1:nTissues
        for k = 1:nSim
            meanDose(i, j, k) = output(par(:, colDose) == tDose(i),...
                1 + tSim(k), j, k);
            stdDose(i, j, k) = output(par(:, colDose) == tDose(i),...
                3 + tSim(k), j, k);
        end
        
    end
end

for i = 1:length(tDose)
    nfig = nfig + 1;
    figure(nfig)
    ax = gca;
    ax.FontSize = 22;
    ax = gca;
    hold on
    for j = 1:nSim
        plot(initVascDens', permute(meanDose(i, tissuesVascDens, j),...
            [2, 3, 1]), '-o', 'Linewidth', 4, 'Markersize', 6)
    end
    for j = 1:nSim
        errorbar(initVascDens', permute(meanDose(i, tissuesVascDens, j),...
            [2, 3, 1]), permute(stdDose(i, tissuesVascDens, j),...
            [2, 3, 1]), '.k', 'Markersize', 10)
    end
    hold off
    ylim([0, inf])
    grid on
    title(['All tissues - ', char(outputNames(selOut)), ' - '...
        num2str(tDose(i)), ' Gy'], 'FontSize', 22)
    xlabel('Initial vascular density (%)', 'FontSize', 22)
    ylabel(char(outputNames(selOut)), 'FontSize', 22)
    legend(simN, 'location', 'southeast', 'FontSize', 22)
end

%%
% This block plots, for all the values of dose, the relative difference of
% the selected output as a function of the initial vascular density

meanDose = zeros(length(tDose), nTissues, nSim);
stdDose = zeros(length(tDose), nTissues, nSim);

for i = 1:length(tDose)
    for j = 1:nTissues
        for k = 1:nSim
            meanDose(i, j, k) = output(par(:, colDose) == tDose(i),...
                1 + tSim(k), j, k);
            stdDose(i, j, k) = output(par(:, colDose) == tDose(i),...
                3 + tSim(k), j, k);
        end
        
    end
end

for i = 1:length(tDose)
    nfig = nfig + 1;
    figure(nfig)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 2;
    ax.FontSize = 22;
    for j = 2:nSim
        diff = (meanDose(i, tissuesVascDens, j) -...
            meanDose(i, tissuesVascDens, 1)) ./...
            meanDose(i, tissuesVascDens, 1) ;
        diff = permute(diff, [2, 3, 1]);
        plot(initVascDens', diff, '-o', 'Linewidth', 2, 'Markersize', 10)
    end
    hold off
    grid on
    title(['All tissues - ', 'Relative difference in ',...
        char(outputNames(selOut)), ' - ', num2str(tDose(i)), ' Gy'],...
        'FontSize', 22)
    xlabel('Initial vascular density (%)', 'FontSize', 22)
    ylabel(['Relative difference in ', char(outputNames(selOut)),],...
        'FontSize', 22)
    legend(simN(2:end), 'location', 'northeast', 'FontSize', 22)
end