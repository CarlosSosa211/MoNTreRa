clear all
close all

nfig = 0;

% path = '../../Carlos/Results/Err_Ang_250Sim';
% path = '../../Carlos/Results/Err_AlphaBeta_432Sim';
% path = '../../Carlos/Results/Err_AlphaBeta_2592Sim';
path = '../../Carlos/Results/Err_AngAlphaBeta_2592Sim';

nOut = 8;
errEndTreatTumDens = load([path, '/errEndTreatTumDens.res']);
err3MonTumDens = load([path, '/err3MonTumDens.res']);
errRecTumDens = load([path, '/errRecTumDens.res']);
errTumVol = load([path, '/errFinTumVol.res']);
errIntTumDens = load([path, '/errIntTumDens.res']);
errTimeTo95 = load([path, '/errTimeTo95.res']);
errTimeTo99 = load([path, '/errTimeTo99.res']);
errRecTime = load([path, '/errRecTime.res']);

par = errEndTreatTumDens(:, 3:6);

colTTum = 1;
colDThres = 2;
colTArrest = 3;
colDose = 4;

tTTum = unique(par(:, colTTum));
tDThres = unique(par(:, colDThres));
tTArrest = unique(par(:, colTArrest));
tDose = unique(par(:, colDose));

%%

b = {'$endTreatTumDens$', '$3MonTumDens$', '$recTumDens$', '$tumVol$'...
    '$intTumDens$', '$timeTo95$', '$timeTo99$', '$recTime$'};

errRel = [num2cell(mean(errEndTreatTumDens(:, 2))), num2cell(mean(err3MonTumDens(:, 2)))...
    num2cell(mean(errRecTumDens(:, 2))), num2cell(mean(errTumVol(:, 2)))...
    num2cell(mean(errIntTumDens(:, 2))), num2cell(mean(errTimeTo95(:, 2)))...
    num2cell(mean(errTimeTo99(:, 2))), num2cell(mean(errRecTime(:, 2)));
    num2cell(std(errEndTreatTumDens(:, 2))), num2cell(std(err3MonTumDens(:, 2)))...
    num2cell(std(errRecTumDens(:, 2))), num2cell(std(errTumVol(:, 2)))...
    num2cell(std(errIntTumDens(:, 2))), num2cell(std(errTimeTo95(:, 2)))...
    num2cell(std(errTimeTo99(:, 2))), num2cell(std(errRecTime(:, 2)));
    b];

nfig = nfig + 1;
figure(nfig);
hold on
bar(cell2mat(errRel(1, :)));
errorbar(1:nOut, cell2mat(errRel(1, :)), cell2mat(errRel(2, :)), '.k')
hold off

ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nOut)
set(ax,'XTickLabel', errRel(3, :), 'fontsize', 20);
xtickangle(45)
ax.YGrid = 'on';
ylim([0, inf])
title('Relative errors', 'fontsize', 20)

%%
errEndTreatTumDensMeanTTum = [];
errEndTreatTumDensStdTTum = [];
for i = 1:length(tTTum)
    errEndTreatTumDensMeanTTum = [errEndTreatTumDensMeanTTum...
        mean(errEndTreatTumDens(par(:, colTTum) == tTTum(i), 1))];
    errEndTreatTumDensStdTTum = [errEndTreatTumDensStdTTum...
        std(errEndTreatTumDens(par(:, colTTum) == tTTum(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, errEndTreatTumDensMeanTTum, errEndTreatTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density at the end of treat.')
grid on
xlabel('TTum (h)')
ylabel('Error')

%%
errEndTreatTumDensMeanDThres = [];
errEndTreatTumDensStdDThres = [];
for i = 1:length(tDThres)
    errEndTreatTumDensMeanDThres = [errEndTreatTumDensMeanDThres...
        mean(errEndTreatTumDens(par(:, colDThres) == tDThres(i), 1))];
    errEndTreatTumDensStdDThres = [errEndTreatTumDensStdDThres...
        std(errEndTreatTumDens(par(:, colDThres) == tDThres(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDThres, errEndTreatTumDensMeanDThres, errEndTreatTumDensStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density at the end of treat.')
grid on
xlabel('Dose threshold(Gy)')
ylabel('Error')

%%
errEndTreatTumDensMeanTArrest = [];
errEndTreatTumDensStdTArrest = [];
for i = 1:length(tTArrest)
    errEndTreatTumDensMeanTArrest = [errEndTreatTumDensMeanTArrest...
        mean(errEndTreatTumDens(par(:, colTArrest) == tTArrest(i), 1))];
    errEndTreatTumDensStdTArrest = [errEndTreatTumDensStdTArrest...
        std(errEndTreatTumDens(par(:, colTArrest) == tTArrest(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTArrest, errEndTreatTumDensMeanTArrest, errEndTreatTumDensStdTArrest,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density at the end of treat.')
grid on
xlabel('TArrest (h)')
ylabel('Error')

%%
errEndTreatTumDensMeanDose = [];
errEndTreatTumDensStdDose = [];
for i = 1:length(tDose)
    errEndTreatTumDensMeanDose = [errEndTreatTumDensMeanDose...
        mean(errEndTreatTumDens(par(:, colDose) == tDose(i), 1))];
    errEndTreatTumDensStdDose = [errEndTreatTumDensStdDose...
        std(errEndTreatTumDens(par(:, colDose) == tDose(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, errEndTreatTumDensMeanDose, errEndTreatTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density at the end of treat.')
grid on
xlabel('Dose (Gy)')
ylabel('Error')

%%
err3MonTumDensMeanTTum = [];
err3MonTumDensStdTTum = [];
for i = 1:length(tTTum)
    err3MonTumDensMeanTTum = [err3MonTumDensMeanTTum...
        mean(err3MonTumDens(par(:, colTTum) == tTTum(i), 1))];
    err3MonTumDensStdTTum = [err3MonTumDensStdTTum...
        std(err3MonTumDens(par(:, colTTum) == tTTum(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, err3MonTumDensMeanTTum, err3MonTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density 3 months after the end of treat.')
grid on
xlabel('TTum (h)')
ylabel('Error')

%%
err3MonTumDensMeanDThres = [];
err3MonTumDensStdDThres = [];
for i = 1:length(tDThres)
    err3MonTumDensMeanDThres = [err3MonTumDensMeanDThres...
        mean(err3MonTumDens(par(:, colDThres) == tDThres(i), 1))];
    err3MonTumDensStdDThres = [err3MonTumDensStdDThres...
        std(err3MonTumDens(par(:, colDThres) == tDThres(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDThres, err3MonTumDensMeanDThres, err3MonTumDensStdDThres,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density 3 months after the end of treat.')
grid on
xlabel('Dose threshold(Gy)')
ylabel('Error')

%%
err3MonTumDensMeanDose = [];
err3MonTumDensStdDose = [];
for i = 1:length(tDose)
    err3MonTumDensMeanDose = [err3MonTumDensMeanDose...
        mean(err3MonTumDens(par(:, colDose) == tDose(i), 1))];
    err3MonTumDensStdDose = [err3MonTumDensStdDose...
        std(err3MonTumDens(par(:, colDose) == tDose(i), 1))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, err3MonTumDensMeanDose, err3MonTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Abs. error in tumor density 3 months after the end of treat.')
grid on
xlabel('Dose (Gy)')
ylabel('Error')

%%
errRecTimeMeanTTum = [];
errRecTimeStdTTum = [];
for i = 1:length(tTTum)
    errRecTimeMeanTTum = [errRecTimeMeanTTum...
        mean(errRecTime(par(:, colTTum) == tTTum(i), 2))];
    errRecTimeStdTTum = [errRecTimeStdTTum...
        std(errRecTime(par(:, colTTum) == tTTum(i), 2))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, errRecTimeMeanTTum, errRecTimeStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Rel. error in time to recurrence')
grid on
xlabel('TTum (h)')
ylabel('Error')

%%
errEndTreatTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        errEndTreatTumDensTTumDose(i, j) =...
            mean(errEndTreatTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
errEndTreatTumDensDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        errEndTreatTumDensDThresDose(i, j) =...
            mean(errEndTreatTumDens(par(:, colDThres) == tDThres(i) &...
            par(:, colDose) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensDThresDose, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('Dose (Gy)')
ylabel('Dthres (Gy)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tDThres))
yticklabels(tDThres)

%%
errEndTreatTumDensTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        errEndTreatTumDensTTumDThres(i, j) =...
            mean(errEndTreatTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDThres) == tDThres(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensTTumDThres, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('DThres (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
errEndTreatTumDensTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        errEndTreatTumDensTTumDThres(i, j) =...
            mean(errEndTreatTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDThres) == tDThres(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensTTumDThres, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('DThres (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
errIntTumDensTArrestDose = zeros(length(tTArrest), length(tDose));
for i = 1:length(tTArrest)
    for j = 1:length(tDose)
        errIntTumDensTArrestDose(i, j) =...
            mean(errIntTumDens(par(:, colTArrest) == tTArrest(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errIntTumDensTArrestDose, 'CDataMapping','scaled')
colorbar
title('Rel. error in integral of tumor density')
xlabel('Dose (Gy)')
ylabel('TArrest (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTArrest))
yticklabels(tTArrest)

%%
errIntTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        errIntTumDensTTumDose(i, j) =...
            mean(errIntTumDens(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errIntTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title('Rel. error in integral of tumor density')
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
errTimeTo95TTumDens = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        errTimeTo95TTumDose(i, j) =...
            mean(errTimeTo95(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errTimeTo95TTumDose, 'CDataMapping','scaled')
colorbar
title('Rel. error in time to kill 95% of tumor cells')
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
errRecTimeTTumDens = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        errRecTimeTTumDose(i, j) =...
            mean(errRecTime(par(:, colTTum) == tTTum(i) &...
            par(:, colDose) == tDose(j), 2));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errRecTimeTTumDose, 'CDataMapping','scaled')
colorbar
title('Rel. error in time to recurrence')
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
nSim = 2592;
nErr = 0.5 * nSim;

for i = 1:nErr
    clear tumDens
    
    tumDens(:, :, 1) = load([path, '/tumDens_', num2str(2 * (i - 1)), '.res']);
    tumDens(:, :, 2) = load([path, '/tumDens_', num2str(2 * i - 1), '.res']);
    
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

for i = 1:length(tTTum)
    simTumDensMeanTTum = [simTumDensMeanTTum...
        mean(simTumDens(par(:, colTTum) == tTTum(i)))];
    simTumDensStdTTum = [simTumDensStdTTum...
        std(simTumDens(par(:, colTTum) == tTTum(i)))];
    
    normInfTumDensMeanTTum = [normInfTumDensMeanTTum...
        mean(normInfTumDens(par(:, colTTum) == tTTum(i)))];
    normInfTumDensStdTTum = [normInfTumDensStdTTum...
        std(normInfTumDens(par(:, colTTum) == tTTum(i)))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, simTumDensMeanTTum, simTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Similarity in tumor density')
grid on
xlabel('TTum (h)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tTTum, normInfTumDensMeanTTum, normInfTumDensStdTTum,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Infitiny norm of abs. error in tumor density')
grid on
xlabel('TTum (h)')
ylabel('||Error||_\infty')

%%
simTumDensMeanDose = [];
simTumDensStdDose = [];
normInfTumDensMeanDose = [];
normInfTumDensStdDose = [];

for i = 1:length(tDose)
    simTumDensMeanDose = [simTumDensMeanDose...
        mean(simTumDens(par(:, colDose) == tDose(i)))];
    simTumDensStdDose = [simTumDensStdDose...
        std(simTumDens(par(:, colDose) == tDose(i)))];
    
    normInfTumDensMeanDose = [normInfTumDensMeanDose...
        mean(normInfTumDens(par(:, colDose) == tDose(i)))];
    normInfTumDensStdDose = [normInfTumDensStdDose...
        std(normInfTumDens(par(:, colDose) == tDose(i)))];
end

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, simTumDensMeanDose, simTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Similarity in tumor density')
grid on
xlabel('Dose (Gy)')
ylabel('Similarity')

nfig = nfig + 1;
figure(nfig)
errorbar(tDose, normInfTumDensMeanDose, normInfTumDensStdDose,...
    '-s', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
title('Infitiny norm of abs. error in tumor density')
grid on
xlabel('Dose (Gy)')
ylabel('||Error||_\infty')

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
title('Similarity in tumor density')
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
title('Similarity in tumor density')
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
title('Infitiny norm of abs. error in tumor density')
xlabel('Dose (Gy)')
ylabel('TTum (h)')
xticks(1:length(tDose))
xticklabels(tDose)
yticks(1:length(tTTum))
yticklabels(tTTum)

%%
nSim = 2592;
nErr = 0.5 * nSim;

% hold on

for i = 1:nErr
    clear tumDens tumVol vascDens killedCells cycle hypDens pO2Stat vegfStat
    
    %     tumDens(:, :, 1) = load([path, '/tumDens_', num2str(2 * (i - 1)), '.res']);
    %     tumDens(:, :, 2) = load([path, '/tumDens_', num2str(2 * i - 1), '.res']);
    %     errTumDens(i) = tumDens(:, 2, 1)' * tumDens(:, 2, 2) /...
    %         (norm(tumDens(:, 2, 1)) * norm(tumDens(:, 2, 2)));
    %     RTumDens(:, :, i) = corrcoef(tumDens(:, 2, 1), tumDens(:, 2, 2));
    %     [pTumDens(i, :), S] = polyfit(tumDens(:, 2, 1)', tumDens(:, 2, 2)', 1);
    %     [tumDensFit(i, :), dTumDens(i)] = polyval(pTumDens(i, :), tumDens(:, 2, 1)', S);
    % todos los tumDensFit no son de la misma longitud por lo que no se pueden
    % guardar como filas de una misma matriz
    %     polyval(pTumDens(i, :), tumDens(:, 2, 1)', S);
    %     plot(tumDens(:, 2, 1)', tumDens(:, 2, 2)', 'o')
    
    %     tumVol(:, :, 1) = load([path, '/tumVol_', num2str(2 * (i - 1)), '.res']);
    %     tumVol(:, :, 2) = load([path, '/tumVol_', num2str(2 * i - 1), '.res']);
    %     errTumVol(i) = tumVol(:, 2, 1)' * tumVol(:, 2, 2) /...
    %         (norm(tumVol(:, 2, 1)) * norm(tumVol(:, 2, 2)));
    
    %     cycle(:, :, 1) = load([path, '/cycle_', num2str(2 * (i - 1)), '.res']);
    %     cycle(:, :, 2) = load([path, '/cycle_', num2str(2 * i - 1), '.res']);
    %     errG1(i) = cycle(:, 2, 1)' * cycle(:, 2, 2) /...
    %         (norm(cycle(:, 2, 1)) * norm(cycle(:, 2, 2)));
    %     RG1(:, :, i) = corrcoef(cycle(:, 2, 1), cycle(:, 2, 2));
    %     pG1(i, :) = polyfit(cycle(:, 2, 1), cycle(:, 2, 2), 1);
    %     plot(cycle(:, 2, 1), cycle(:, 2, 2))
    
    %     hypDens(:, :, 1) = load([path, '/hypDens_', num2str(2 * (i - 1)), '.res']);
    %     hypDens(:, :, 2) = load([path, '/hypDens_', num2str(2 * i - 1), '.res']);
    %     errHypDens(i) = hypDens(:, 2, 1)' * hypDens(:, 2, 2) /...
    %         (norm(hypDens(:, 2, 1)) * norm(hypDens(:, 2, 2)));
    %     RHypDens(:, :, i) = corrcoef(hypDens(:, 2, 1), hypDens(:, 2, 2));
    %     pHypDens(i, :) = polyfit(hypDens(:, 2, 1), hypDens(:, 2, 2), 1);
    %     plot(hypDens(:, 2, 1), hypDens(:, 2, 2))
    
    pO2Stat(:, :, 1) = load([path, '/pO2Stat_', num2str(2 * (i - 1)), '.res']);
    pO2Stat(:, :, 2) = load([path, '/pO2Stat_', num2str(2 * i - 1), '.res']);
    errpO2Mean(i) = pO2Stat(:, 2, 1)' * pO2Stat(:, 2, 2) /...
        (norm(pO2Stat(:, 2, 1)) * norm(pO2Stat(:, 2, 2)));
    errpO2Med(i) = pO2Stat(:, 3, 1)' * pO2Stat(:, 3, 2) /...
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
% plot(tumDens(:, 2, 1)', tumDensFit(nErr, :), 'r-')
% plot(tumDens(:, 2, 1)', tumDensFit(nErr, :) + 2 * dTumDens(nErr, :),'m--',...
%     tumDens(:, 2, 1), tumDensFit(nErr, :) - 2 * dTumDens(nErr, :), 'm--')
% title('Linear Fit of Data with 95% Prediction Interval')
% legend('Data','Linear Fit','95% Prediction Interval')
% meanErrTumDens = mean(errTumDens)
% stdErrTumDens = std(errTumDens)
% meanRTumDens = mean(RTumDens, 3)
% stdRTumDens = std(RTumDens, 0, 3)
% meanpTumDens = mean(pTumDens)
% stdpTumDens = std(pTumDens)

% meanErrHypDens = mean(errHypDens)
% stdErrHypDens = std(errHypDens)
% meanRHypDens = mean(RHypDens, 3)
% stdRHypDens = std(RHypDens, 0, 3)
% meanpHypDens = mean(pHypDens)
% stdpHypDens = std(pHypDens)

% meanErrG1 = mean(errG1)
% stdErrG1 = std(errG1)
% meanRG1 = mean(RG1, 3)
% stdRG1 = std(RG1, 0, 3)
% meanpG1 = mean(pG1)
% stdpG1 = std(pG1)

meanErrpO2Mean = mean(errpO2Mean)
stdErrpO2Mean = std(errpO2Mean)
meanRpO2Mean = mean(RpO2Mean, 3)
stdRpO2Mean = std(RpO2Mean, 0, 3)
% meanppO2Mean = mean(ppO2Mean)
% stdppO2Mean = std(ppO2Mean)
meanErrpO2Med = mean(errpO2Med)
stdErrpO2Med = std(errpO2Med)
meanRpO2Med = mean(RpO2Med, 3)
stdRpO2Med = std(RpO2Med, 0, 3)
meannormInfpO2Mean = mean(normInfpO2Mean)
stdnormInfpO2Mean = std(normInfpO2Mean)
meannormInfpO2Med = mean(normInfpO2Med)
stdnormInfpO2Med = std(normInfpO2Med)

% vegfStat(:, :, 1) = load([path, '/vegfStat_', num2str(2 * (i - 1)), '.res']);
% vegfStat(:, :, 2) = load([path, '/vegfStat_', num2str(2 * i - 1), '.res']);
% errVegfMean(i) = vegfStat(:, 2, 1)' * vegfStat(:, 2, 2) /...
%     (norm(vegfStat(:, 2, 1)) * norm(vegfStat(:, 2, 2)));
% errVegfMed(i) = vegfStat(:, 3, 1)' * vegfStat(:, 3, 2) /...
%     (norm(vegfStat(:, 3, 1)) * norm(vegfStat(:, 3, 2)));