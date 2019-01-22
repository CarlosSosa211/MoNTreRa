clear all
close all

nfig = 0;

nTissue = 4;
% path = ['../../Carlos/Results/Diff_Ang_432Sim_AllTissues/Tissue'...
%     num2str(nTissue)];
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

selOut = input(['Select an output [tumDens (1), tumVol (2)'...
    'vascDens (3), preExVascDens (4), neoCreVascDens (5)\n'...
    'killedCells (6), hypDens (7), pO2Med (8), pO2Mean (9), '...
    'vegfMed (10), vegfMean (11), distG1 (12),\n'...
    'distS (13),  distG2 (14), distM (15) or distG0 (16)] or quit (0): ']);

switch selOut
    case 1
        path = [path, '/tumDens_'];
        outputCol = 2;
        outputName = 'tumour density';
    case 2
        path = [path, '/tumVol_'];
        outputCol = 2;
        outputName = 'tumour volume';
    case 3
        path = [path, '/vascDens_'];
        outputCol = 2;
        outputName = 'vascular density';
    case 4
        path = [path, '/vascDens_'];
        outputCol = 3;
        outputName = 'pre-existing vascular density';
    case 5
        path = [path, '/vascDens_'];
        outputCol = 4;
        outputName = 'neo-created vascular density';
    case 6
        path = [path, '/killedCell_'];
        outputCol = 2;
        outputName = 'killed cells';
    case 7
        path = [path, '/hypDens_'];
        outputCol = 2;
        outputName = 'hypoxic density';
    case 8
        path = [path, '/pO2Stat_'];
        outputCol = 2;
        outputName = 'median pO2';
    case 9
        path = [path, '/pO2Stat_'];
        outputCol = 3;
        outputName = 'mean pO2';
    case 10
        path = [path, '/vegfStat_'];
        outputCol = 2;
        outputName = 'mean VEGF concentration';
    case 11
        path = [path, '/vegfStat_'];
        outputCol = 3;
        outputName = 'median VEGF concentration';
    case 12
        path = [path, '/cycle_'];
        outputCol = 2;
        outputName = 'G1 distribution';
    case 13
        path = [path, '/cycle_'];
        outputCol = 3;
        outputName = 'S distribution';
    case 14
        path = [path, '/cycle_'];
        outputCol = 4;
        outputName = 'G2 distribution';
    case 15
        path = [path, '/cycle_'];
        outputCol = 5;
        outputName = 'M distribution';
    case 16
        path = [path, '/cycle_'];
        outputCol = 6;
        outputName = 'G0 distribution';
end

%%
meanOutput0 = {};
meanOutput1 = {};
stdOutput0 = {};
stdOutput1 = {};
P = 5;
for i = 1:size(par, 1)
    clear output0 output1
    for j = 1:P
        temp0 = load([path, num2str(i - 1), '_0_', num2str(j - 1)...
            '.res']);
        temp1 = load([path, num2str(i - 1), '_1_', num2str(j - 1)...
            '.res']);
        output0(:, :, j) = temp0(:, [1, outputCol]);
        output1(:, :, j) = temp1(:, [1, outputCol]);
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
    
    [p(i, :), S] = polyfit(meanOutput0i(:, 2)', meanOutput1i(:, 2)', 1);
    
    normInf(i) = norm(meanOutput0i(:, 2) - meanOutput1i(:, 2), 'inf');
    
end

%%
meanSim = mean(sim);
stdSim = std(sim);
meanR = mean(R, 3);
stdR = std(R, 0, 3);
meanp = mean(p);
stdp = std(p);
meanNorm = mean(normInf);
stdNormInf = std(normInf);

%%
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
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, lsdMeanTTum, lsdStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - LSD in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('LSD')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, normInfMeanTTum, normInfStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in '...
    outputName])
grid on
xlabel('TTum (h)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p1MeanTTum, p1StdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('p1')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p0MeanTTum, p0StdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p0 in ', outputName])
grid on
xlabel('TTum (h)')
ylabel('p0')

%%
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
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, lsdMeanDose, lsdStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - LSD in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('LSD')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, normInfMeanDose, normInfStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in '...
    outputName])
grid on
xlabel('Dose (Gy)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p1MeanDose, p1StdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('p1')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p0MeanDose, p0StdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p0 in in ', outputName])
grid on
xlabel('Dose (Gy)')
ylabel('p0')

%%
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
i = 212;

mean0 = cell2mat(meanOutput0(i));
mean1 = cell2mat(meanOutput1(i));

std0 = cell2mat(stdOutput0(i));
std1 = cell2mat(stdOutput1(i));

nfig = nfig + 1;
figure(nfig)
hold on
% plot(mean0(:, 1), mean0(:, 2));
% plot(mean1(:, 1), mean1(:, 2));
errorbar(mean0(:, 1), mean0(:, 2), std0(:, 2));
errorbar(mean1(:, 1), mean1(:, 2), std1(:, 2));
hold off
xlabel('Time (h)')
ylabel(char(outputName))
legend('No angiogenesis', 'Angiogenesis', 'location', 'northwest')
