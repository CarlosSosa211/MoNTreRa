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

nOut = 8;
diffEndTreatTumDens = load([path, '/errEndTreatTumDens.res']);
diff3MonTumDens = load([path, '/err3MonTumDens.res']);
diffRecTumDens = load([path, '/errRecTumDens.res']);
diffTumVol = load([path, '/errFinTumVol.res']);
diffIntTumDens = load([path, '/errIntTumDens.res']);
diffTimeTo95 = load([path, '/errTimeTo95.res']);
diffTimeTo99 = load([path, '/errTimeTo99.res']);
diffRecTime = load([path, '/errRecTime.res']);

par = diffEndTreatTumDens(:, 3:5);

colTTum = 1;
colDThres = 2;
% colTArrest = 3;
colDose = 3;

tTTum = unique(par(:, colTTum));
tDThres = unique(par(:, colDThres));
% tTArrest = unique(par(:, colTArrest));
tDose = unique(par(:, colDose));

%%

b = {'$endTreatTumDens$', '$3MonTumDens$', '$recTumDens$', '$tumVol$'...
    '$intTumDens$', '$timeTo95$', '$timeTo99$', '$recTime$'};

diffRel = [num2cell(mean(diffEndTreatTumDens(:, 2))), num2cell(mean(diff3MonTumDens(:, 2)))...
    num2cell(mean(diffRecTumDens(:, 2))), num2cell(mean(diffTumVol(:, 2)))...
    num2cell(mean(diffIntTumDens(:, 2))), num2cell(mean(diffTimeTo95(:, 2)))...
    num2cell(mean(diffTimeTo99(:, 2))), num2cell(mean(diffRecTime(:, 2)));
    num2cell(std(diffEndTreatTumDens(:, 2))), num2cell(std(diff3MonTumDens(:, 2)))...
    num2cell(std(diffRecTumDens(:, 2))), num2cell(std(diffTumVol(:, 2)))...
    num2cell(std(diffIntTumDens(:, 2))), num2cell(std(diffTimeTo95(:, 2)))...
    num2cell(std(diffTimeTo99(:, 2))), num2cell(std(diffRecTime(:, 2)));
    b];

nfig = nfig + 1;
figure(nfig);
hold on
bar(cell2mat(diffRel(1, :)));
errorbar(1:nOut, cell2mat(diffRel(1, :)), cell2mat(diffRel(2, :)), '.k')
hold off

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nOut)
set(ax,'XTickLabel', diffRel(3, :), 'fontsize', 20);
xtickangle(45)
ax.YGrid = 'on';
ylim([0, inf])
title(['Tissue ', num2str(nTissue), ' - Relative differences'], 'fontsize', 20)

%%
diffEndTreatTumDensMeanTTum = [];
diffEndTreatTumDensStdTTum = [];
for i = 1:length(tTTum)
    diffEndTreatTumDensMeanTTum = [diffEndTreatTumDensMeanTTum...
        mean(diffEndTreatTumDens(par(:, colTTum) == tTTum(i), 1))];
    diffEndTreatTumDensStdTTum = [diffEndTreatTumDensStdTTum...
        std(diffEndTreatTumDens(par(:, colTTum) == tTTum(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, diffEndTreatTumDensMeanTTum, diffEndTreatTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
grid on
xlabel('TTum (h)')
ylabel('Difference')

%%
diffEndTreatTumDensMeanDThres = [];
diffEndTreatTumDensStdDThres = [];
for i = 1:length(tDThres)
    diffEndTreatTumDensMeanDThres = [diffEndTreatTumDensMeanDThres...
        mean(diffEndTreatTumDens(par(:, colDThres) == tDThres(i), 1))];
    diffEndTreatTumDensStdDThres = [diffEndTreatTumDensStdDThres...
        std(diffEndTreatTumDens(par(:, colDThres) == tDThres(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDThres, diffEndTreatTumDensMeanDThres, diffEndTreatTumDensStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
grid on
xlabel('Dose threshold(Gy)')
ylabel('Difference')

%%
diffEndTreatTumDensMeanTArrest = [];
diffEndTreatTumDensStdTArrest = [];
for i = 1:length(tTArrest)
    diffEndTreatTumDensMeanTArrest = [diffEndTreatTumDensMeanTArrest...
        mean(diffEndTreatTumDens(par(:, colTArrest) == tTArrest(i), 1))];
    diffEndTreatTumDensStdTArrest = [diffEndTreatTumDensStdTArrest...
        std(diffEndTreatTumDens(par(:, colTArrest) == tTArrest(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTArrest, diffEndTreatTumDensMeanTArrest, diffEndTreatTumDensStdTArrest,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
grid on
xlabel('TArrest (h)')
ylabel('Difference')

%%
diffEndTreatTumDensMeanDose = [];
diffEndTreatTumDensStdDose = [];
for i = 1:length(tDose)
    diffEndTreatTumDensMeanDose = [diffEndTreatTumDensMeanDose...
        mean(diffEndTreatTumDens(par(:, colDose) == tDose(i), 1))];
    diffEndTreatTumDensStdDose = [diffEndTreatTumDensStdDose...
        std(diffEndTreatTumDens(par(:, colDose) == tDose(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, diffEndTreatTumDensMeanDose, diffEndTreatTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
grid on
xlabel('Dose (Gy)')
ylabel('Difference')

%%
diff3MonTumDensMeanTTum = [];
diff3MonTumDensStdTTum = [];
for i = 1:length(tTTum)
    diff3MonTumDensMeanTTum = [diff3MonTumDensMeanTTum...
        mean(diff3MonTumDens(par(:, colTTum) == tTTum(i), 1))];
    diff3MonTumDensStdTTum = [diff3MonTumDensStdTTum...
        std(diff3MonTumDens(par(:, colTTum) == tTTum(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, diff3MonTumDensMeanTTum, diff3MonTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density 3 months after the end of treat.'])
grid on
xlabel('TTum (h)')
ylabel('Difference')

%%
diff3MonTumDensMeanDThres = [];
diff3MonTumDensStdDThres = [];
for i = 1:length(tDThres)
    diff3MonTumDensMeanDThres = [diff3MonTumDensMeanDThres...
        mean(diff3MonTumDens(par(:, colDThres) == tDThres(i), 1))];
    diff3MonTumDensStdDThres = [diff3MonTumDensStdDThres...
        std(diff3MonTumDens(par(:, colDThres) == tDThres(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDThres, diff3MonTumDensMeanDThres, diff3MonTumDensStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density 3 months after the end of treat.'])
grid on
xlabel('Dose threshold(Gy)')
ylabel('Difference')

%%
diff3MonTumDensMeanDose = [];
diff3MonTumDensStdDose = [];
for i = 1:length(tDose)
    diff3MonTumDensMeanDose = [diff3MonTumDensMeanDose...
        mean(diff3MonTumDens(par(:, colDose) == tDose(i), 1))];
    diff3MonTumDensStdDose = [diff3MonTumDensStdDose...
        std(diff3MonTumDens(par(:, colDose) == tDose(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, diff3MonTumDensMeanDose, diff3MonTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density 3 months after the end of treat.'])
grid on
xlabel('Dose (Gy)')
ylabel('Difference')

%%
diffRecTimeMeanTTum = [];
diffRecTimeStdTTum = [];
for i = 1:length(tTTum)
    diffRecTimeMeanTTum = [diffRecTimeMeanTTum...
        mean(diffRecTime(par(:, colTTum) == tTTum(i), 2))];
    diffRecTimeStdTTum = [diffRecTimeStdTTum...
        std(diffRecTime(par(:, colTTum) == tTTum(i), 2))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, diffRecTimeMeanTTum, diffRecTimeStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Rel. diff. in time to recurrence'])
grid on
xlabel('TTum (h)')
ylabel('Difference')

%%
diffEndTreatTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        diffEndTreatTumDensTTumDose(i, j) =...
            mean(diffEndTreatTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffEndTreatTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffEndTreatTumDensDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        diffEndTreatTumDensDThresDose(i, j) =...
            mean(diffEndTreatTumDens(par(:, colDThres) == tDThres(i) &...
            par(:, colDose) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffEndTreatTumDensDThresDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
xlabel('Dose (Gy)')
ylabel('Dthres (Gy)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tDThres))
yticklabels(tDThres)

%%
diffEndTreatTumDensTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        diffEndTreatTumDensTTumDThres(i, j) =...
            mean(diffEndTreatTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDThres) == tDThres(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffEndTreatTumDensTTumDThres, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
xlabel('DThres (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffEndTreatTumDensTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        diffEndTreatTumDensTTumDThres(i, j) =...
            mean(diffEndTreatTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDThres) == tDThres(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffEndTreatTumDensTTumDThres, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Abs. diff. in tumour density at the end of treat.'])
xlabel('DThres (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffIntTumDensTArrestDose = zeros(length(tTArrest), length(tDose));
for i = 1:length(tTArrest)
    for j = 1:length(tDose)
        diffIntTumDensTArrestDose(i, j) =...
            mean(diffIntTumDens(par(:, colTArrest) == tTArrest(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffIntTumDensTArrestDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Rel. diff. in integral of tumour density'])
xlabel('Dose (Gy)')
ylabel('TArrest (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTArrest))
yticklabels(tTArrest)

%%
diffIntTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        diffIntTumDensTTumDose(i, j) =...
            mean(diffIntTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffIntTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Rel. diff. in integral of tumour density'])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffTimeTo95TTumDens = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        diffTimeTo95TTumDose(i, j) =...
            mean(diffTimeTo95(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffTimeTo95TTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Rel. diff. in time to kill 95% of tumour cells'])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
diffRecTimeTTumDens = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        diffRecTimeTTumDose(i, j) =...
            mean(diffRecTime(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(diffRecTimeTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Rel. diff. in time to recurrence'])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
nSim = 432;
ndiff = 0.5 * nSim;

for i = 1:ndiff
    clear tumDens
    tumDens(:, :, 1) = load([path, '/tumDens_', num2str(2 * (i - 1)), '.res']);
    tumDens(:, :, 2) = load([path, '/tumDens_', num2str(2 * i - 1), '.res']);
    
    simTime = size(tumDens, 1);
    endTreatTime = simTime - 120;
    
    endTreatTumDens(i, 1) = tumDens(endTreatTime, 2, 1);
    endTreatTumDens(i, 2) = tumDens(endTreatTime, 2, 2);
    
    simTumDens(i) = tumDens(:, 2, 1)' * tumDens(:, 2, 2) /...
        (norm(tumDens(:, 2, 1)) * norm(tumDens(:, 2, 2)));
    
    RTumDens(:, :, i) = corrcoef(tumDens(:, 2, 1), tumDens(:, 2, 2));
    
    [pTumDens(i, :), S] = polyfit(tumDens(:, 2, 1)', tumDens(:, 2, 2)', 1);
    
    normInfTumDens(i) = norm(tumDens(:, 2, 1) - tumDens(:, 2, 2), 'inf');
end

%%
meanSimTumDens = mean(simTumDens);
stdSimTumDens = std(simTumDens);
meanRTumDens = mean(RTumDens, 3);
stdRTumDens = std(RTumDens, 0, 3);
meanpTumDens = mean(pTumDens);
stdpTumDens = std(pTumDens);
meannormTumDens = mean(normInfTumDens);
stdnormInfTumDens = std(normInfTumDens);

%%
simTumDensMeanTTum = [];
simTumDensStdTTum = [];
normInfTumDensMeanTTum = [];
normInfTumDensStdTTum = [];
p1TumDensMeanTTum = [];
p1TumDensStdTTum = [];
p0TumDensMeanTTum = [];
p0TumDensStdTTum = [];

for i = 1:length(tTTum)
    simTumDensMeanTTum = [simTumDensMeanTTum...
        mean(simTumDens(par(:, colTTum) == tTTum(i)))];
    simTumDensStdTTum = [simTumDensStdTTum...
        std(simTumDens(par(:, colTTum) == tTTum(i)))];
    
    normInfTumDensMeanTTum = [normInfTumDensMeanTTum...
        mean(normInfTumDens(par(:, colTTum) == tTTum(i)))];
    normInfTumDensStdTTum = [normInfTumDensStdTTum...
        std(normInfTumDens(par(:, colTTum) == tTTum(i)))];
    
    p1TumDensMeanTTum = [p1TumDensMeanTTum...
        mean(pTumDens(par(:, colTTum) == tTTum(i), 1))];
    p1TumDensStdTTum = [p1TumDensStdTTum...
        std(pTumDens(par(:, colTTum) == tTTum(i), 1))];
    
    p0TumDensMeanTTum = [p0TumDensMeanTTum...
        mean(pTumDens(par(:, colTTum) == tTTum(i), 2))];
    p0TumDensStdTTum = [p0TumDensStdTTum...
        std(pTumDens(par(:, colTTum) == tTTum(i), 2))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, simTumDensMeanTTum, simTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in tumour density'])
grid on
xlabel('TTum (h)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, normInfTumDensMeanTTum, normInfTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in tumour density'])
grid on
xlabel('TTum (h)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p1TumDensMeanTTum, p1TumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in tumour density'])
grid on
xlabel('TTum (h)')
ylabel('p1')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p0TumDensMeanTTum, p0TumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p0 in tumour density'])
grid on
xlabel('TTum (h)')
ylabel('p0')

%%
simTumDensMeanDose = [];
simTumDensStdDose = [];
normInfTumDensMeanDose = [];
normInfTumDensStdDose = [];
p1TumDensMeanDose = [];
p1TumDensStdDose = [];
p0TumDensMeanDose = [];
p0TumDensStdDose = [];

for i = 1:length(tDose)
    simTumDensMeanDose = [simTumDensMeanDose...
        mean(simTumDens(par(:, colDose) == tDose(i)))];
    simTumDensStdDose = [simTumDensStdDose...
        std(simTumDens(par(:, colDose) == tDose(i)))];
    
    normInfTumDensMeanDose = [normInfTumDensMeanDose...
        mean(normInfTumDens(par(:, colDose) == tDose(i)))];
    normInfTumDensStdDose = [normInfTumDensStdDose...
        std(normInfTumDens(par(:, colDose) == tDose(i)))];
    
    p1TumDensMeanDose = [p1TumDensMeanDose...
        mean(pTumDens(par(:, colDose) == tDose(i), 1))];
    p1TumDensStdDose = [p1TumDensStdDose...
        std(pTumDens(par(:, colDose) == tDose(i), 1))];
    
    p0TumDensMeanDose = [p0TumDensMeanDose...
        mean(pTumDens(par(:, colDose) == tDose(i), 2))];
    p0TumDensStdDose = [p0TumDensStdDose...
        std(pTumDens(par(:, colDose) == tDose(i), 2))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, simTumDensMeanDose, simTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in tumour density'])
grid on
xlabel('Dose (Gy)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, normInfTumDensMeanDose, normInfTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in tumour density'])
grid on
xlabel('Dose (Gy)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p1TumDensMeanDose, p1TumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in tumour density'])
grid on
xlabel('Dose (Gy)')
ylabel('p1')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p0TumDensMeanDose, p0TumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p0 in tumour density'])
grid on
xlabel('Dose (Gy)')
ylabel('p0')

%%
simTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        simTumDensTTumDose(i, j) =...
            mean(simTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(simTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Similarity in tumour density'])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
simTumDensDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        simTumDensDThresDose(i, j) =...
            mean(simTumDens(par(:, colDThres) == tDThres(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(simTumDensDThresDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Similarity in tumour density'])
xlabel('Dose (Gy)')
ylabel('Dthres (Gy)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tDThres))
yticklabels(tDThres)

%%
normInfTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        normInfTumDensTTumDose(i, j) =...
            mean(normInfTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j)));
    end
end

nfig = nfig + 1;
figure(nfig)
image(normInfTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in tumour density'])
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
nSim = 2592;
ndiff = 0.5 * nSim;

% hold on

for i = 1:ndiff
    clear tumVol vascDens killedCells cycle hypDens pO2Stat vegfStat
    
    %     tumVol(:, :, 1) = load([path, '/tumVol_', num2str(2 * (i - 1)), '.res']);
    %     tumVol(:, :, 2) = load([path, '/tumVol_', num2str(2 * i - 1), '.res']);
    %     diffTumVol(i) = tumVol(:, 2, 1)' * tumVol(:, 2, 2) /...
    %         (norm(tumVol(:, 2, 1)) * norm(tumVol(:, 2, 2)));
    
    %     cycle(:, :, 1) = load([path, '/cycle_', num2str(2 * (i - 1)), '.res']);
    %     cycle(:, :, 2) = load([path, '/cycle_', num2str(2 * i - 1), '.res']);
    %     diffG1(i) = cycle(:, 2, 1)' * cycle(:, 2, 2) /...
    %         (norm(cycle(:, 2, 1)) * norm(cycle(:, 2, 2)));
    %     RG1(:, :, i) = corrcoef(cycle(:, 2, 1), cycle(:, 2, 2));
    %     pG1(i, :) = polyfit(cycle(:, 2, 1), cycle(:, 2, 2), 1);
    %     plot(cycle(:, 2, 1), cycle(:, 2, 2))
    
    %     hypDens(:, :, 1) = load([path, '/hypDens_', num2str(2 * (i - 1)), '.res']);
    %     hypDens(:, :, 2) = load([path, '/hypDens_', num2str(2 * i - 1), '.res']);
    %     diffHypDens(i) = hypDens(:, 2, 1)' * hypDens(:, 2, 2) /...
    %         (norm(hypDens(:, 2, 1)) * norm(hypDens(:, 2, 2)));
    %     RHypDens(:, :, i) = corrcoef(hypDens(:, 2, 1), hypDens(:, 2, 2));
    %     pHypDens(i, :) = polyfit(hypDens(:, 2, 1), hypDens(:, 2, 2), 1);
    %     plot(hypDens(:, 2, 1), hypDens(:, 2, 2))
    
    pO2Stat(:, :, 1) = load([path, '/pO2Stat_', num2str(2 * (i - 1)), '.res']);
    pO2Stat(:, :, 2) = load([path, '/pO2Stat_', num2str(2 * i - 1), '.res']);
    diffpO2Mean(i) = pO2Stat(:, 2, 1)' * pO2Stat(:, 2, 2) /...
        (norm(pO2Stat(:, 2, 1)) * norm(pO2Stat(:, 2, 2)));
    diffpO2Med(i) = pO2Stat(:, 3, 1)' * pO2Stat(:, 3, 2) /...
        (norm(pO2Stat(:, 3, 1)) * norm(pO2Stat(:, 3, 2)));
    normInfpO2Mean(i) = norm(pO2Stat(:, 2, 1) - pO2Stat(:, 2, 2), 'inf');
    normInfpO2Med(i) = norm(pO2Stat(:, 3, 1) - pO2Stat(:, 3, 2), 'inf');
    RpO2Mean(:, :, i) = corrcoef(pO2Stat(:, 2, 1), pO2Stat(:, 2, 2));
    RpO2Med(:, :, i) = corrcoef(pO2Stat(:, 3, 1), pO2Stat(:, 3, 2));
    
end

% hold off
%
% plot(tumDens(:, 2, 1)', tumDens(:, 2, 2)','bo')
% hold on
% plot(tumDens(:, 2, 1)', tumDensFit(ndiff, :), 'r-')
% plot(tumDens(:, 2, 1)', tumDensFit(ndiff, :) + 2 * dTumDens(ndiff, :),'m--',...
%     tumDens(:, 2, 1), tumDensFit(ndiff, :) - 2 * dTumDens(ndiff, :), 'm--')
% title('Linear Fit of Data with 95% Prediction Interval')
% legend('Data','Linear Fit','95% Prediction Interval')
% meandiffTumDens = mean(diffTumDens)
% stddiffTumDens = std(diffTumDens)
% meanRTumDens = mean(RTumDens, 3)
% stdRTumDens = std(RTumDens, 0, 3)
% meanpTumDens = mean(pTumDens)
% stdpTumDens = std(pTumDens)

% meandiffHypDens = mean(diffHypDens)
% stddiffHypDens = std(diffHypDens)
% meanRHypDens = mean(RHypDens, 3)
% stdRHypDens = std(RHypDens, 0, 3)
% meanpHypDens = mean(pHypDens)
% stdpHypDens = std(pHypDens)

% meandiffG1 = mean(diffG1)
% stddiffG1 = std(diffG1)
% meanRG1 = mean(RG1, 3)
% stdRG1 = std(RG1, 0, 3)
% meanpG1 = mean(pG1)
% stdpG1 = std(pG1)

meandiffpO2Mean = mean(diffpO2Mean)
stddiffpO2Mean = std(diffpO2Mean)
meanRpO2Mean = mean(RpO2Mean, 3)
stdRpO2Mean = std(RpO2Mean, 0, 3)
% meanppO2Mean = mean(ppO2Mean)
% stdppO2Mean = std(ppO2Mean)
meandiffpO2Med = mean(diffpO2Med)
stddiffpO2Med = std(diffpO2Med)
meanRpO2Med = mean(RpO2Med, 3)
stdRpO2Med = std(RpO2Med, 0, 3)
meannormInfpO2Mean = mean(normInfpO2Mean)
stdnormInfpO2Mean = std(normInfpO2Mean)
meannormInfpO2Med = mean(normInfpO2Med)
stdnormInfpO2Med = std(normInfpO2Med)

%%
nSim = 2592;
ndiff = 0.5 * nSim;

nfig = nfig + 1;
figure(nfig)
hold on
for i = 1:ndiff
    clear pO2Stat
    
    pO2Stat(:, :, 1) = load([path, '/pO2Stat_', num2str(2 * (i - 1)), '.res']);
    pO2Stat(:, :, 2) = load([path, '/pO2Stat_', num2str(2 * i - 1), '.res']);
    
    simpO2Med(i) = pO2Stat(:, 2, 1)' * pO2Stat(:, 2, 2) /...
        (norm(pO2Stat(:, 2, 1)) * norm(pO2Stat(:, 2, 2)));
    simpO2Mean(i) = pO2Stat(:, 3, 1)' * pO2Stat(:, 3, 2) /...
        (norm(pO2Stat(:, 3, 1)) * norm(pO2Stat(:, 3, 2)));
    
    RpO2Med(:, :, i) = corrcoef(pO2Stat(:, 2, 1), pO2Stat(:, 2, 2));
    RpO2Mean(:, :, i) = corrcoef(pO2Stat(:, 3, 1), pO2Stat(:, 3, 2));
    
    [ppO2Med(i, :), S] = polyfit(pO2Stat(:, 2, 1)', pO2Stat(:, 2, 2)', 1);
    [ppO2Mean(i, :), S] = polyfit(pO2Stat(:, 3, 1)', pO2Stat(:, 3, 2)', 1);
    
    normInfpO2Med(i) = norm(pO2Stat(:, 2, 1) - pO2Stat(:, 2, 2), 'inf');
    normInfpO2Mean(i) = norm(pO2Stat(:, 3, 1) - pO2Stat(:, 3, 2), 'inf');
    
    plot(pO2Stat(:, 2, 1)', pO2Stat(:, 2, 2)', 'o')
end

hold off
%%
meanSimpO2Med = mean(simpO2Med);
stdSimpO2Med = std(simpO2Med);
meanRpO2Med = mean(RpO2Med, 3);
stdRpO2Med = std(RpO2Med, 0, 3);
meanppO2Med = mean(ppO2Med);
stdppO2Med = std(ppO2Med);
meannormpO2Med = mean(normInfpO2Med);
stdnormInfpO2Med = std(normInfpO2Med);

%%
simpO2MedMeanTTum = [];
simpO2MedStdTTum = [];
normInfpO2MedMeanTTum = [];
normInfpO2MedStdTTum = [];
p1pO2MedMeanTTum = [];
p1pO2MedStdTTum = [];

for i = 1:length(tTTum)
    simpO2MedMeanTTum = [simpO2MedMeanTTum...
        mean(simpO2Med(par(:, colTTum) == tTTum(i)))];
    simpO2MedStdTTum = [simpO2MedStdTTum...
        std(simpO2Med(par(:, colTTum) == tTTum(i)))];
    
    normInfpO2MedMeanTTum = [normInfpO2MedMeanTTum...
        mean(normInfpO2Med(par(:, colTTum) == tTTum(i)))];
    normInfpO2MedStdTTum = [normInfpO2MedStdTTum...
        std(normInfpO2Med(par(:, colTTum) == tTTum(i)))];
    
    p1pO2MedMeanTTum = [p1pO2MedMeanTTum...
        mean(ppO2Med(par(:, colTTum) == tTTum(i), 1))];
    p1pO2MedStdTTum = [p1pO2MedStdTTum...
        std(ppO2Med(par(:, colTTum) == tTTum(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, simpO2MedMeanTTum, simpO2MedStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in pO2 med'])
grid on
xlabel('TTum (h)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, normInfpO2MedMeanTTum, normInfpO2MedStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in pO2 med'])
grid on
xlabel('TTum (h)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p1pO2MedMeanTTum, p1pO2MedStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in pO2 med'])
grid on
xlabel('TTum (h)')
ylabel('p1')


%%
simpO2MedMeanDose = [];
simpO2MedStdDose = [];
normInfpO2MedMeanDose = [];
normInfpO2MedStdDose = [];
p1pO2MedMeanDose = [];
p1pO2MedStdDose = [];

for i = 1:length(tDose)
    simpO2MedMeanDose = [simpO2MedMeanDose...
        mean(simpO2Med(par(:, colDose) == tDose(i)))];
    simpO2MedStdDose = [simpO2MedStdDose...
        std(simpO2Med(par(:, colDose) == tDose(i)))];
    
    normInfpO2MedMeanDose = [normInfpO2MedMeanDose...
        mean(normInfpO2Med(par(:, colDose) == tDose(i)))];
    normInfpO2MedStdDose = [normInfpO2MedStdDose...
        std(normInfpO2Med(par(:, colDose) == tDose(i)))];
    
    p1pO2MedMeanDose = [p1pO2MedMeanDose...
        mean(ppO2Med(par(:, colDose) == tDose(i), 1))];
    p1pO2MedStdDose = [p1pO2MedStdDose...
        std(ppO2Med(par(:, colDose) == tDose(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, simpO2MedMeanDose, simpO2MedStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in pO2 med'])
grid on
xlabel('Dose (Gy)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, normInfpO2MedMeanDose, normInfpO2MedStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in pO2 med'])
grid on
xlabel('Dose (Gy)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, p1pO2MedMeanDose, p1pO2MedStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in pO2 med'])
grid on
xlabel('Dose (Gy)')
ylabel('p1')

%%
nSim = 2592;
ndiff = 0.5 * nSim;

for i = 1:ndiff
    clear vegfStat
    
    vegfStat(:, :, 1) = load([path, '/vegfStat_', num2str(2 * (i - 1)), '.res']);
    vegfStat(:, :, 2) = load([path, '/vegfStat_', num2str(2 * i - 1), '.res']);
    
    simVegfMean(i) = vegfStat(:, 2, 1)' * vegfStat(:, 2, 2) /...
        (norm(vegfStat(:, 2, 1)) * norm(vegfStat(:, 2, 2)));
    simVegfMed(i) = vegfStat(:, 3, 1)' * vegfStat(:, 3, 2) /...
        (norm(vegfStat(:, 3, 1)) * norm(vegfStat(:, 3, 2)));
    
    RVegfMean(:, :, i) = corrcoef(vegfStat(:, 2, 1), vegfStat(:, 2, 2));
    RVegfMed(:, :, i) = corrcoef(vegfStat(:, 3, 1), vegfStat(:, 3, 2));
    
    [pVegfMean(i, :), S] = polyfit(vegfStat(:, 2, 1)', vegfStat(:, 2, 2)', 1);
    [pVegfMed(i, :), S] = polyfit(vegfStat(:, 3, 1)', vegfStat(:, 3, 2)', 1);
    
    normInfVegfMean(i) = norm(vegfStat(:, 2, 1) - vegfStat(:, 2, 2), 'inf');
    normInfVegfMed(i) = norm(vegfStat(:, 3, 1) - vegfStat(:, 3, 2), 'inf');
end

%%
meanSimVegfMean = mean(simVegfMean);
stdSimVegfMean = std(simVegfMean);
meanRVegfMean = mean(RVegfMean, 3);
stdRVegfMean = std(RVegfMean, 0, 3);
meanpVegfMean = mean(pVegfMean);
stdpVegfMean = std(pVegfMean);
meannormVegfMean = mean(normInfVegfMean);
stdnormInfVegfMean = std(normInfVegfMean);

%%
simVegfMeanMeanTTum = [];
simVegfMeanStdTTum = [];
normInfVegfMeanMeanTTum = [];
normInfVegfMeanStdTTum = [];
p1VegfMeanMeanTTum = [];
p1VegfMeanStdTTum = [];

for i = 1:length(tTTum)
    simVegfMeanMeanTTum = [simVegfMeanMeanTTum...
        mean(simVegfMean(par(:, colTTum) == tTTum(i)))];
    simVegfMeanStdTTum = [simVegfMeanStdTTum...
        std(simVegfMean(par(:, colTTum) == tTTum(i)))];
    
    normInfVegfMeanMeanTTum = [normInfVegfMeanMeanTTum...
        mean(normInfVegfMean(par(:, colTTum) == tTTum(i)))];
    normInfVegfMeanStdTTum = [normInfVegfMeanStdTTum...
        std(normInfVegfMean(par(:, colTTum) == tTTum(i)))];
    
    p1VegfMeanMeanTTum = [p1VegfMeanMeanTTum...
        mean(pVegfMean(par(:, colTTum) == tTTum(i), 1))];
    p1VegfMeanStdTTum = [p1VegfMeanStdTTum...
        std(pVegfMean(par(:, colTTum) == tTTum(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, simVegfMeanMeanTTum, simVegfMeanStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Similarity in VEGF mean'])
grid on
xlabel('TTum (h)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, normInfVegfMeanMeanTTum, normInfVegfMeanStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - Infitiny norm of abs. diff. in VEGF mean'])
grid on
xlabel('TTum (h)')
ylabel('||diff.||_\infty')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, p1VegfMeanMeanTTum, p1VegfMeanStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title(['Tissue ', num2str(nTissue), ' - p1 in VEGF mean'])
grid on
xlabel('TTum (h)')
ylabel('p1')

%%
%     [tumDensFit(i, :), dTumDens(i)] = polyval(pTumDens(i, :), tumDens(:, 2, 1)', S);
% todos los tumDensFit no son de la misma longitud por lo que no se pueden
% guardar como filas de una misma matriz
%     polyval(pTumDens(i, :), tumDens(:, 2, 1)', S);
%     plot(tumDens(:, 2, 1)', tumDens(:, 2, 2)', 'o')

%%
nfig = nfig + 1;
figure(nfig)
hold on
histogram(endTreatTumDens(:, 1), 25);
histogram(endTreatTumDens(:, 2), 25);
hold off
title(['Tissue ', num2str(nTissue), ' - Tumour density at the end of treat.']);
legend('No angiogenesis', 'Angiogenesis')
%
% p = ranksum(endTreatTumDens(:, 1),endTreatTumDens(:, 2));
dEndTreatTumDens = endTreatTumDens(:, 1)' - endTreatTumDens(:, 2)';
medianDEndTreatTumDens = median(dEndTreatTumDens)
meanDEndTreatTumDens = mean(dEndTreatTumDens)
%skewDEndTreatTumDens = skewness(dEndTreatTumDens)
binSize = 1/16;
bound = 3;
lBound = - bound - 0.5 * binSize;
rBound = bound + 0.5 * binSize;
edges = [min([dEndTreatTumDens, lBound]), lBound:binSize:rBound, max([dEndTreatTumDens, rBound])];
nfig = nfig + 1;
figure(nfig)
histogram(dEndTreatTumDens, edges);
title(['Tissue ', num2str(nTissue), ' - Difference in tumour density at the end of treat.'])

%%
nfig = nfig + 1;
figure(nfig)
hold on
plot(endTreatTumDens(:, 1));
plot(endTreatTumDens(:, 2));
hold off
legend('No angiogenesis', 'Angiogenesis')

%%
countGreat = 0;
countLess = 0;
countEqual = 0;
s = sign(dEndTreatTumDens);
ipositif = sum(s == 1)
inegatif = sum(s == -1)
countGreat
countLess
countEqual
median(endTreatTumDens(:, 1))
median(endTreatTumDens(:, 2))
mean(endTreatTumDens(:, 1))
mean(endTreatTumDens(:, 2))