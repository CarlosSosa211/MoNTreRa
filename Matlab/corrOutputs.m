clear all
close all

%%
% This block initialises nfig and defines
% - path: the path of the ouput directory,
% - nTissues: the number of tissues,
% - fileNames: the names of the ouput files,
% - ouputNames: the names of the outputs,
% - outputCol: the column number of the output in the corresponding files.

nfig = 0;
path = '../../Carlos/Results/Corr/Dose_5Val_5Rep/Tissue';
% Simulations were performed using the input files *Res.dat. To study the
% results of simulation considering all biological processes, output files
% X_1_X.res are used.
nTissues = 21;
P = 5;

fileNames = {'/tumDens','/vascDens', '/vascDens', '/vascDens'...
    '/deadDens', '/pO2Stat', '/pO2Stat'};

outputNames = {'Tumour density', 'Vascular density'...
    'Pre-existing vascular density', 'Neo-created vascular density'...
    'Dead cells density', 'Median pO2', 'Mean pO2'};

outputCol = [2, 2, 3, 4, 2, 2, 3];

%%
% This block loads the initial vascular densities values for the 21
% tissues.

initVascDens = load('../HistSpec/initVascDens.dat');
[initVascDens, tissuesVascDens] = sort(initVascDens);

%%
% This block studies the correlation between two time-dependent outputs
% selOut1 and selOut2. For every combination of parameters (values of
% dose), the following indicators are calculated:
% - sim: <u, v> / (||u|| * ||v||),
% - lsd: Sum of ui^2 - vi^2,
% - R: the correlation coefficients,
% - p: the polynomial of degree 1 fitting the curve u = v.

% For the case of nTissues = 0, these indicators are calculated
% independently for each tissue.

% The following intermidiate variables are considered:
% - par: a matrix containing the combinations of parameters simulated. Each
% row corresponds to a simulation and each column, to a parameter,
% - meanOutput1i: a matrix containing the mean values for P repetitions of
% output 1 for the combination of parameters i,
% - meanOutput2i: a matrix containing the mean values for P repetitions of
% output 2 for the combination of parameters i,
% - meanOutput1: a cell array containing the mean values for P repetitions
% of output 1 for all the combination of parameters,
% - meanOutput2: a cell array containing the mean values for P repetitions
% of output 2 for all the combination of parameters,
% - stdOutput1: a cell array containing the std values for P repetitions of
% output 1 for all the combination of parameters,
% - stdOutput2: a cell array containing the std values for P repetitions of
% output 2 for all the combination of parameters.

nTissue = input(['Select one tissue (from 1 to ', num2str(nTissues)...
    ') or all of them (0) or a mean over them (-1): ']);

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7)] or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7)] or quit (0): ']);


if(nTissue >= 1 && nTissue <= nTissues)
    pathTissue = [path, num2str(nTissue)];
    par = load([pathTissue, '/combPar.res']);
    nCombPar = size(par, 1);
    colDose = 1;
    tDose = unique(par(:, colDose));
    
    meanOutput1 = cell(1, nCombPar);
    meanOutput2 = cell(1, nCombPar);
    stdOutput1 = cell(1, nCombPar);
    stdOutput2 = cell(1, nCombPar);
    
    sim = zeros(nCombPar, 1);
    lsd = zeros(nCombPar, 1);
    R = zeros(2, 2, nCombPar);
    p = zeros(nCombPar, 2);
    
    for i = 1:nCombPar
        clear output1 output2
        for j = 1:P
            temp1 = load([pathTissue, char(fileNames(selOut1)), '_'...
                num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            temp2 = load([pathTissue, char(fileNames(selOut2)), '_'...
                num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            output1(:, :, j) = temp1(:, [1, outputCol(selOut1)]);
            output2(:, :, j) = temp2(:, [1, outputCol(selOut2)]);
        end
        
        meanOutput1i = mean(output1, 3);
        meanOutput2i = mean(output2, 3);
        meanOutput1(i) = {meanOutput1i};
        meanOutput2(i) = {meanOutput2i};
        
        stdOutput1(i) = {std(output1, 0, 3)};
        stdOutput2(i) = {std(output2, 0, 3)};
        
        sim(i) = meanOutput1i(2:end, 2)' * meanOutput2i(2:end, 2) /...
            (norm(meanOutput1i(2:end, 2)) * norm(meanOutput2i(2:end, 2)));
        
        lsd(i) = sum((meanOutput1i(2:end, 2) - meanOutput2i(2:end, 2)).^2);
        
        R(:, :, i) = corrcoef(meanOutput1i(2:end, 2),...
            meanOutput2i(2:end, 2));
        
        [p(i, :), S] = polyfit(meanOutput1i(2:end, 2)',...
            meanOutput2i(2:end, 2)', 1);
    end
end

if(nTissue == 0)
    par = load([path, '1/combPar.res']);
    nCombPar = size(par, 1);
    meanOutput1 = cell(1, nCombPar);
    meanOutput2 = cell(1, nCombPar);
    stdOutput1 = cell(1, nCombPar);
    stdOutput2 = cell(1, nCombPar);
    colDose = 1;
    tDose = unique(par(:, colDose));
    sim = zeros(nCombPar, nTissues);
    lsd = zeros(nCombPar, nTissues);
    R = zeros(2, 2, nCombPar, nTissues);
    p = zeros(nCombPar, 2, nTissues);
    
    for k = 1:nTissues
        pathTissue = [path, num2str(k)];
        for i = 1:nCombPar
            clear output1 output2
            for j = 1:P
                temp1 = load([pathTissue, char(fileNames(selOut1)), '_'...
                    num2str(i - 1), '_1_', num2str(j - 1), '.res']);
                temp2 = load([pathTissue, char(fileNames(selOut2)), '_'...
                    num2str(i - 1), '_1_', num2str(j - 1), '.res']);
                output1(:, :, j) = temp1(:, [1, outputCol(selOut1)]);
                output2(:, :, j) = temp2(:, [1, outputCol(selOut2)]);
            end
            
            meanOutput1i = mean(output1, 3);
            meanOutput2i = mean(output2, 3);
            meanOutput1(i) = {meanOutput1i};
            meanOutput2(i) = {meanOutput2i};
            
            stdOutput1(i) = {std(output1, 0, 3)};
            stdOutput2(i) = {std(output2, 0, 3)};
            
            sim(i, k) = meanOutput1i(2:end, 2)' *...
                meanOutput2i(2:end, 2) / (norm(meanOutput1i(2:end, 2)) *...
                norm(meanOutput2i(2:end, 2)));
            
            lsd(i, k) = sum((meanOutput1i(2:end, 2) -...
                meanOutput2i(2:end, 2)).^2);
            
            R(:, :, i, k) = corrcoef(meanOutput1i(2:end, 2),...
                meanOutput2i(2:end, 2));
            
            [p(i, :, k), S] = polyfit(meanOutput1i(2:end, 2)',...
                meanOutput2i(2:end, 2)', 1);
        end
    end
end
%%
pAllDoses = mean(p, 1);
pAllDosesStd = std(p, 0, 1);
pAllDoses = permute(pAllDoses, [3, 2, 1]);
pAllDosesStd = permute(pAllDosesStd, [3, 2, 1]);
pAllTissues = pAllDoses ./ [initVascDens, initVascDens];

%%
% This block plots, for given combination of parameters (dose value), the
% mean values of time-dependent output 1 as a function of the mean values
% of time-dependent output 2.

dose = 1;
meanOutput1Dose = cell2mat(meanOutput1(dose));
meanOutput2Dose = cell2mat(meanOutput2(dose));

nfig = nfig + 1;
figure(nfig)
plot(meanOutput1Dose(:, 2), meanOutput2Dose(:, 2), '-o')
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

%%
% This block fits with a polynomial of degree 1 the initial values of
% time-dependent output 1 and output 2 for all the tissues. The result
% is then plotted.

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

output1 = zeros(1, nTissues);
output2 = zeros(1, nTissues);
for i = 1:nTissues
    pathTissue = [path, num2str(i)];
    temp = load([pathTissue, char(fileNames(selOut1)),'_0_1_0.res']);
    output1(i) = temp(2, outputCol(selOut1));
    temp = load([pathTissue, char(fileNames(selOut2)),'_0_1_0.res']);
    output2(i) = temp(2, outputCol(selOut2));
end

[sortOutput1, ind] = sort(output1);
sortOutput2 = output2(ind);
p = polyfit(sortOutput1, sortOutput2, 2);

nfig = nfig + 1;
figure(nfig)
hold on
plot(sortOutput1, sortOutput2, '-o')
alpha = 0.2;
beta = 0.6;
plot(sortOutput1, sortOutput1.^(5/3))
% plot(sortOutput1, p(1) * sortOutput1 + p(2))
hold off
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

%%
% This block fits with a polynomial of degree 1 the time-dependent output 1
% and output 2. The result is then plotted.

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

meanP = zeros(nTissues, 2);
stdP = zeros(nTissues, 2);

for nTissue = 1:nTissues
    pathTissue = [path, num2str(nTissue)];
    par = load([pathTissue, '/combPar.res']);
    nCombPar = size(par, 1);
    
    meanOutput1 = cell(1, nCombPar);
    meanOutput2 = cell(1, nCombPar);
    stdOutput1 = cell(1, nCombPar);
    stdOutput2 = cell(1, nCombPar);
    
    p = zeros(nCombPar, 2);
    
    for i = 1:nCombPar
        for j = 1:P
            temp1 = load([pathTissue, char(fileNames(selOut1))...
                '_', num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            temp2 = load([pathTissue, char(fileNames(selOut2))...
                '_', num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            output1(:, :, j) = temp1(:, [1, outputCol(selOut1)]);
            output2(:, :, j) = temp2(:, [1, outputCol(selOut2)]);
        end
        output1 = mean(output1, 3);
        output2 = mean(output2, 3);
        meanOutput1(i) = {output1};
        meanOutput2(i) = {output2};
        
        p(i, :) = polyfit(output1(2:end, 2), output2(2:end, 2), 1);
    end
    
    meanP(nTissue, :) = mean(p, 1);
    stdP(nTissue, :) = std(p, 0, 1);
end

for i = 1:nCombPar
    nfig = nfig + 1;
    figure(nfig)
    hold on
    ax = gca;
    ax.FontSize = 22;
    output1 = cell2mat(meanOutput1(i));
    output2 = cell2mat(meanOutput2(i));
    plot(output1(2:end, 2), output2(2:end, 2))
    plot(output1(2:end, 2), meanP(nTissue, 1) * output1(2:end, 2) +...
        meanP(nTissue, 2))
    hold off
    grid on
    xlabel(char(outputNames(selOut1)), 'FontSize', 22)
    ylabel(char(outputNames(selOut2)), 'FontSize', 22)
    title(['Tissue ', num2str(nTissue), ' - ', num2str(par(i))...
        ' Gy'], 'FontSize', 22)
    
    nfig = nfig + 1;
    figure(nfig)
    hold on
    ax = gca;
    ax.FontSize = 22;
    plot(output1(2:end, 1), output1(2:end, 2))
    plot(output2(2:end, 1), output2(2:end, 2))
    plot(output2(2:end, 1), meanP(nTissue, 1) * output1(2:end, 2) +...
        meanP(nTissue, 2))
    hold off
    grid on
    xlabel('t (h)', 'FontSize', 22)
    ylabel(char(outputNames(selOut2)), 'FontSize', 22)
    title(['Tissue ', num2str(nTissue), ' - ', num2str(par(i))...
        ' Gy'], 'FontSize', 22)
    legend({char(outputNames(selOut1)), char(outputNames(selOut2)),...
        ['Fitted ', char(outputNames(selOut2))]}, 'FontSize', 22)
    
end

%%
% This block fits with a polynomial of degree 1 the time-dependent output 1
% and output 2. The result is then plotted in a same figure for all the
% tissues.

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

pathTissue = [path, '1'];
par = load([pathTissue, '/combPar.res']);
nCombPar = size(par, 1);
p = zeros(nCombPar, 2);

meanOutput1 = cell(nTissues, nCombPar);
meanOutput2 = cell(nTissues, nCombPar);
stdOutput1 = cell(nTissues, nCombPar);
stdOutput2 = cell(nTissues, nCombPar);

meanP = zeros(nTissues, 2);
stdP = zeros(nTissues, 2);

nrowPlot = 3;
ncolPlot = ceil(nTissues / nrowPlot);

for nTissue = 1:nTissues
    pathTissue = [path, num2str(nTissue)];
    p = zeros(nCombPar, 2);
    
    for i = 1:nCombPar
        for j = 1:P
            temp1 = load([pathTissue, char(fileNames(selOut1))...
                '_', num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            temp2 = load([pathTissue, char(fileNames(selOut2))...
                '_', num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            output1(:, :, j) = temp1(:, [1, outputCol(selOut1)]);
            output2(:, :, j) = temp2(:, [1, outputCol(selOut2)]);
        end
        output1 = mean(output1, 3);
        output2 = mean(output2, 3);
        meanOutput1(nTissue, i) = {output1};
        meanOutput2(nTissue, i) = {output2};
        
        p(i, :) = polyfit(output1(2:end, 2), output2(2:end, 2), 1);
    end
    
    meanP(nTissue, :) = mean(p, 1);
    stdP(nTissue, :) = std(p, 0, 1);
end

%%
meanTisP = mean(meanP);
stdTisP = std(meanP);

for i = 1:nCombPar
    nfig = nfig + 1;
    figure(nfig)
    nsubplot = 1;
    for nTissue = tissuesVascDens'
        subplot(nrowPlot, ncolPlot, nsubplot)
        nsubplot = nsubplot + 1;
        hold on
        output1 = cell2mat(meanOutput1(nTissue, i));
        output2 = cell2mat(meanOutput2(nTissue, i));
        plot(output1(2:end, 2), output2(2:end, 2))
        plot(output1(2:end, 2), meanTisP(1) * output1(2:end, 2) +...
            meanTisP(2))
        hold off
        grid on
        xlabel(char(outputNames(selOut1)))
        ylabel(char(outputNames(selOut2)))
        title(['Tissue ', num2str(nTissue), ' - ', num2str(par(i))...
            ' Gy'])
    end
    
    nfig = nfig + 1;
    figure(nfig)
    nsubplot = 1;
    for nTissue = tissuesVascDens'
        subplot(nrowPlot, ncolPlot, nsubplot)
        nsubplot = nsubplot + 1;
        hold on
        output1 = cell2mat(meanOutput1(nTissue, i));
        output2 = cell2mat(meanOutput2(nTissue, i));
        plot(output1(2:end, 1), output1(2:end, 2))
        plot(output2(2:end, 1), output2(2:end, 2))
        plot(output2(2:end, 1), meanTisP(1) * output1(2:end, 2) +...
            meanTisP(2))
        hold off
        grid on
%         ylim([0, inf])
        xlabel('t (h)')
        ylabel(char(outputNames(selOut2)))
        title(['Tissue ', num2str(nTissue), ' - ', num2str(par(i)), ' Gy'])
        %         legend({char(outputNames(selOut1)), char(outputNames(selOut2)),...
        %             ['Fitted ', char(outputNames(selOut2))]})
    end
end