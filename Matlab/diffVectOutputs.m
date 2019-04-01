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
% path = '../../Carlos/Results/Diff/Ang_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/Res_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/AngRes_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/HypNec_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/Ang_10x3Sim_noHypNec/Tissue';
% path = '../../Carlos/Results/Diff/Arrest_Dose_5Val_5Rep/Tissue';
% path = '../../Carlos/Results/Diff/Oxy_Dose_5Val_5Rep/Tissue';
path = '../../Carlos/Results/Diff/OxyNoHypNec_Dose_5Val_5Rep/Tissue';

nTissues = 21;
nOut = 16;
P = 5;

% withoutN = 'No angiogenesis';
% withN = 'Angiogenesis';
% withoutN = 'No healthy cell division';
% withN = 'Healthy cell division';
% withoutN = 'No angiogenesis and no healthy cell division';
% withN = 'Angiogenesis and healthy cell division';
% withoutN = 'No hypoxic necrosis';
% withN = 'Hypoxic necrosis';
% withoutN = 'No arrest';
% withN = 'Arrest';
% withoutN = 'No oxyegenation (no hypoxic necrosis)';
% withN = 'Oxygenation';
withoutN = 'No oxyegenation (no hypoxic necrosis)';
withN = 'Oxygenation (no hypoxic necrosis)';

fileNames = {'/tumDens_', '/tumVol_', '/vascDens_', '/vascDens_'...
    '/vascDens_', '/killedCell_', '/hypDens_', '/pO2Stat_'...
    '/pO2Stat_', '/vegfStat_', '/vegfStat_', '/cycle_', '/cycle_'...
    '/cycle_', '/cycle_', '/cycle_'};

outputNames = {'tumour density', 'tumour volume', 'vascular density'...
    'pre-existing vascular density', 'neo-created vascular density'...
    'killed cells', 'hypoxic density', 'median pO2', 'mean pO2'...
    'mean VEGF concentration', 'median VEGF concentration'...
    'G1 distribution', 'S distribution', 'G2 distribution'...
    'M distribution', 'G0 distribution'};

outputCol = [2, 2, 2, 3, 4, 2, 2, 2, 3, 2, 3, 2, 3, 4, 5, 6];

%%
% This block studies the correlation between the values of a given output
% for X_0_X and X_1_X simulations. For every combination of parameters
% the following indicators are calculated:
% - sim: <u0, u1> / (||u0|| * ||u1||),
% - lsd: Sum of u0i^2 - u1i^2,
% - R: the correlation coefficients,
% - p: the polynomial of degree 1 fitting the curve u0 = u1.
% - normInf: ||u0 - u1||_inf

% For the case of nTissues = 0, these indicators are calculated
% independently for each tissue.

% The following intermidiate variables are considered:
% - par: a matrix containing the combinations of parameters simulated. Each
% row corresponds to a simulation and each column, to a parameter,
% - meanOutput0i: a matrix containing the mean values for P repetitions of
% the output for X_0_X simulations for the combination of parameters i,
% - meanOutput1i: a matrix containing the mean values for P repetitions of
% the output for X_1_X simulations for the combination of parameters i,
% - meanOutput0: a cell array containing the mean values for P repetitions
% of the output for X_0_X simulations for all the combination of
% parameters,
% - meanOutput1: a cell array containing the mean values for P repetitions
% of the output for X_0_X simulations for all the combination of
% parameters,
% - stdOutput1: a cell array containing the std values for P repetitions of
% the output for X_0_X simulations for all the combination of parameters,
% - stdOutput2: a cell array containing the std values for P repetitions of
% the output for X_1_X simulations for all the combination of parameters.

nTissue = input(['Select one tissue (from 1 to ', num2str(nTissues)...
    ') or all of them (0) or a mean over them (-1): ']);

selOut = input(['Select an output [tumDens (1), tumVol (2)'...
    'vascDens (3), preExVascDens (4), neoCreVascDens (5)\n'...
    'killedCells (6), hypDens (7), pO2Med (8), pO2Mean (9), '...
    'vegfMed (10), vegfMean (11), distG1 (12),\n'...
    'distS (13),  distG2 (14), distM (15) or distG0 (16)]'...
    'or quit (0): ']);

if(nTissue >= 1 && nTissue <= nTissues)
    pathTissue = [path, num2str(nTissue)];
    par = load([pathTissue, '/combPar.res']);
    nCombPar = size(par, 1);
    
    % colTTum = 1;
    % colDThres = 2;
    % colTArrest = 3;
    colDose = 1;
    tDose = unique(par(:, colDose));
    % tTTum = unique(par(:, colTTum));
    % tDThres = unique(par(:, colDThres));
    % tTArrest = unique(par(:, colTArrest));
    
    meanOutput0 = cell(1, nCombPar);
    meanOutput1 = cell(1, nCombPar);
    stdOutput0 = cell(1, nCombPar);
    stdOutput1 = cell(1, nCombPar);
    
    sim = zeros(nCombPar, 1);
    lsd = zeros(nCombPar, 1);
    R = zeros(2, 2, nCombPar);
    p = zeros(nCombPar, 2);
    normInf = zeros(nCombPar, 1);
     
    for i = 1:nCombPar
        clear output0 output1
        for j = 1:P
            temp0 = load([pathTissue, char(fileNames(selOut))...
                num2str(i - 1), '_0_', num2str(j - 1), '.res']);
            temp1 = load([pathTissue, char(fileNames(selOut))...
                num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            output0(:, :, j) = temp0(:, [1, outputCol(selOut)]);
            output1(:, :, j) = temp1(:, [1, outputCol(selOut)]);
        end
        
        meanOutput0i = mean(output0, 3);
        meanOutput1i = mean(output1, 3);
        meanOutput0(i) = {meanOutput0i};
        meanOutput1(i) = {meanOutput1i};
        
        stdOutput0(i) = {std(output0, 0, 3)};
        stdOutput1(i) = {std(output1, 0, 3)};
        
        sim(i) = meanOutput0i(:, 2)' * meanOutput1i(:, 2) /...
            (norm(meanOutput0i(:, 2)) * norm(meanOutput1i(:, 2)));
        
        lsd(i) = sum((meanOutput0i(:, 2) - meanOutput1i(:, 2)).^2);
        
        R(:, :, i) = corrcoef(meanOutput0i(:, 2), meanOutput1i(:, 2));
        
        [p(i, :), S] = polyfit(meanOutput0i(:, 2)',...
            meanOutput1i(:, 2)', 1);
        
        normInf(i) = norm(meanOutput0i(:, 2) - meanOutput1i(:, 2), 'inf');
    end
end

if(nTissue == 0)
    for k = 1:nTissues
        pathTissue = [path, num2str(k)];
        par = load([pathTissue, '/combPar.res']);
        nCombPar = size(par, 1);
        
        % colTTum = 1;
        % colDThres = 2;
        % colTArrest = 3;
        colDose = 1;
        tDose = unique(par(:, colDose));
        % tTTum = unique(par(:, colTTum));
        % tDThres = unique(par(:, colDThres));
        % tTArrest = unique(par(:, colTArrest));
        
        meanOutput1 = cell(1, nCombPar);
        meanOutput2 = cell(1, nCombPar);
        stdOutput1 = cell(1, nCombPar);
        stdOutput2 = cell(1, nCombPar);
        
        sim = zeros(nCombPar, nTissues);
        lsd = zeros(nCombPar, nTissues);
        R = zeros(2, 2, nCombPar, nTissues);
        p = zeros(nCombPar, 2, nTissues);
        normInf = zeros(nCombPar, nTissues);
        
        for i = 1:size(par, 1)
            clear output0 output1
            for j = 1:P
                temp0 = load([pathTissue, char(fileNames(selOut))...
                    num2str(i - 1), '_0_', num2str(j - 1), '.res']);
                temp1 = load([pathTissue, char(fileNames(selOut))...
                    num2str(i - 1), '_1_', num2str(j - 1), '.res']);
                output0(:, :, j) = temp0(:, [1, outputCol(selOut)]);
                output1(:, :, j) = temp1(:, [1, outputCol(selOut)]);
            end
            
            meanOutput0i = mean(output0, 3);
            meanOutput1i = mean(output1, 3);
            meanOutput0(i) = {meanOutput0i};
            meanOutput1(i) = {meanOutput1i};
            
            stdOutput0(i) = {std(output0, 0, 3)};
            stdOutput1(i) = {std(output1, 0, 3)};
            
            sim(i, k) = meanOutput0i(:, 2)' * meanOutput1i(:, 2) /...
                (norm(meanOutput0i(:, 2)) * norm(meanOutput1i(:, 2)));
            
            lsd(i, k) = sum((meanOutput0i(:, 2) - meanOutput1i(:, 2)).^2);
            
            R(:, :, i, k) = corrcoef(meanOutput0i(:, 2),...
                meanOutput1i(:, 2));
            
            [p(i, k, :), S] = polyfit(meanOutput0i(:, 2)',...
                meanOutput1i(:, 2)', 1);
            
            normInf(i, k) = norm(meanOutput0i(:, 2) -...
                meanOutput1i(:, 2), 'inf');
        end
    end
end

%%
% This block calculates the mean and std values of all the indicators 

meanSim = mean(sim);
stdSim = std(sim);
meanR = mean(R, 3);
stdR = std(R, 0, 3);
meanp = mean(p);
stdp = std(p);
meanNorm = mean(normInf);
stdNormInf = std(normInf);

%%
% This block plots the mean and std values of all the indicators as a
% function of the studied TTum values

simMeanTTum = [];
simStdTTum = [];
lsdMeanTTum = [];
lsdStdTTum = [];
normInfMeanTTum = [];
normInfStdTTum = [];
p1MeanTTum = [];
p1StdTTum = [];
p0MeanTTum = [];
p0StdTTum = [];

for i = 1:length(tTTum)
    simMeanTTum = [simMeanTTum...
        mean(sim(par(:, colTTum) == tTTum(i)))];
    simStdTTum = [simStdTTum...
        std(sim(par(:, colTTum) == tTTum(i)))];
    
    lsdMeanTTum = [lsdMeanTTum...
        mean(lsd(par(:, colTTum) == tTTum(i)))];
    lsdStdTTum = [lsdStdTTum...
        std(lsd(par(:, colTTum) == tTTum(i)))];
    
    normInfMeanTTum = [normInfMeanTTum...
        mean(normInf(par(:, colTTum) == tTTum(i)))];
    normInfStdTTum = [normInfStdTTum...
        std(normInf(par(:, colTTum) == tTTum(i)))];
    
    p1MeanTTum = [p1MeanTTum...
        mean(p(par(:, colTTum) == tTTum(i), 1))];
    p1StdTTum = [p1StdTTum...
        std(p(par(:, colTTum) == tTTum(i), 1))];
    
    p0MeanTTum = [p0MeanTTum...
        mean(p(par(:, colTTum) == tTTum(i), 2))];
    p0StdTTum = [p0StdTTum...
        std(p(par(:, colTTum) == tTTum(i), 2))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, simMeanTTum, simStdTTum,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, lsdMeanTTum, lsdStdTTum,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - LSD in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('LSD')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, normInfMeanTTum, normInfStdTTum,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in '...
    outputName])
grid on
xlabel('TTum (h)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p1MeanTTum, p1StdTTum,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('p1')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p0MeanTTum, p0StdTTum,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p0 in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('p0')

%%
%%
% This block plots the mean and std values of all the indicators as a
% function of the studied dose values

simMeanDose = [];
simStdDose = [];
lsdMeanDose = [];
lsdStdDose = [];
normInfMeanDose = [];
normInfStdDose = [];
p1MeanDose = [];
p1StdDose = [];
p0MeanDose = [];
p0StdDose = [];

for i = 1:length(tDose)
    simMeanDose = [simMeanDose...
        mean(sim(par(:, colDose) == tDose(i)))];
    simStdDose = [simStdDose...
        std(sim(par(:, colDose) == tDose(i)))];
    
    lsdMeanDose = [lsdMeanDose...
        mean(lsd(par(:, colDose) == tDose(i)))];
    lsdStdDose = [lsdStdDose...
        std(lsd(par(:, colDose) == tDose(i)))];
    
    normInfMeanDose = [normInfMeanDose...
        mean(normInf(par(:, colDose) == tDose(i)))];
    normInfStdDose = [normInfStdDose...
        std(normInf(par(:, colDose) == tDose(i)))];
    
    p1MeanDose = [p1MeanDose...
        mean(p(par(:, colDose) == tDose(i), 1))];
    p1StdDose = [p1StdDose...
        std(p(par(:, colDose) == tDose(i), 1))];
    
    p0MeanDose = [p0MeanDose...
        mean(p(par(:, colDose) == tDose(i), 2))];
    p0StdDose = [p0StdDose...
        std(p(par(:, colDose) == tDose(i), 2))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, simMeanDose, simStdDose,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, lsdMeanDose, lsdStdDose,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - LSD in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('LSD')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, normInfMeanDose, normInfStdDose,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in '...
    outputName])
grid on
xlabel('Dose (Gy)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p1MeanDose, p1StdDose,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('p1')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p0MeanDose, p0StdDose,...
    's', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p0 in in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('p0')

%%
% This block plots the mean of the indicator sim as a function of the
% studied dose and TTum values

simTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        simTTumDose(i, j) =...
            mean(sim(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(simTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Similarity in ', outputName])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
% This block plots the mean of the indicator sim as a function of the
% studied dose and DThres values

simDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        simDThresDose(i, j) =...
            mean(sim(par(:, colDThres) == tDThres(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(simDThresDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Similarity in ', outputName])
xlabel('Dose (Gy)')
ylabel('Dthres (Gy)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tDThres))
yticklabels(tDThres)

%%
% This block plots the mean of the indicator lsd as a function of the
% studied dose and TTum values

lsdTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        lsdTTumDose(i, j) =...
            mean(lsd(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(lsdTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - LSD in ', outputName])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
% This block plots the mean of the indicator normInf as a function of the
% studied dose and TTum values

normInfTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        normInfTTumDose(i, j) =...
            mean(normInf(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(normInfTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in '...
    outputName])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
% This blocks plots the i-th curves contained than meanOutput0 and
% meanOutput1 with the corresponding std.

i = 212;

mean0 = cell2mat(meanOutput0(i));
mean1 = cell2mat(meanOutput1(i));

std0 = cell2mat(stdOutput0(i));
std1 = cell2mat(stdOutput1(i));

nfig = nfig + 1;
figure(nfig)
hold on
errorbar(mean0(:, 1), mean0(:, 2), std0(:, 2));
errorbar(mean1(:, 1), mean1(:, 2), std1(:, 2));
hold off
xlabel('Time (h)')
ylabel(char(outputName))
legend(withoutN, withN, 'location', 'northwest')
