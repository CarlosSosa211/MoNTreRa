clear all
close all

global nfig
nfig = 0;

path = '../../Carlos/Results/Err_Ang_250Sim';

nOut = 8;
errEndTreatTumDens = load([path, '/errEndTreatTumDens.res']);
err3MonTumDens = load([path, '/err3MonTumDens.res']);
errRecTumDens = load([path, '/errRecTumDens.res']);
errTumVol = load([path, '/errFinTumVol.res']);
errIntTumDens = load([path, '/errIntTumDens.res']);
errTimeTo95 = load([path, '/errTimeTo95.res']);
errTimeTo99 = load([path, '/errTimeTo99.res']);
errRecTime = load([path, '/errRecTime.res']);

tTTum = unique(errEndTreatTumDens(:, 3));
tDThres = unique(errEndTreatTumDens(:, 4));
tDose = unique(errEndTreatTumDens(:, 5));

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
        mean(errEndTreatTumDens(errEndTreatTumDens(:, 3) == tTTum(i), 1))];
    errEndTreatTumDensStdTTum = [errEndTreatTumDensStdTTum...
        mean(errEndTreatTumDens(errEndTreatTumDens(:, 3) == tTTum(i), 1))];
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
        mean(errEndTreatTumDens(errEndTreatTumDens(:, 4) == tDThres(i), 1))];
    errEndTreatTumDensStdDThres = [errEndTreatTumDensStdDThres...
        mean(errEndTreatTumDens(errEndTreatTumDens(:, 4) == tDThres(i), 1))];
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
errEndTreatTumDensMeanDose = [];
errEndTreatTumDensStdDose = [];
for i = 1:length(tDose)
    errEndTreatTumDensMeanDose = [errEndTreatTumDensMeanDose...
        mean(errEndTreatTumDens(errEndTreatTumDens(:, 5) == tDose(i), 1))];
    errEndTreatTumDensStdDose = [errEndTreatTumDensStdDose...
        mean(errEndTreatTumDens(errEndTreatTumDens(:, 5) == tDose(i), 1))];
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
        mean(err3MonTumDens(err3MonTumDens(:, 3) == tTTum(i), 1))];
    err3MonTumDensStdTTum = [err3MonTumDensStdTTum...
        mean(err3MonTumDens(err3MonTumDens(:, 3) == tTTum(i), 1))];
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
        mean(err3MonTumDens(err3MonTumDens(:, 4) == tDThres(i), 1))];
    err3MonTumDensStdDThres = [err3MonTumDensStdDThres...
        mean(err3MonTumDens(err3MonTumDens(:, 4) == tDThres(i), 1))];
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
        mean(err3MonTumDens(err3MonTumDens(:, 5) == tDose(i), 1))];
    err3MonTumDensStdDose = [err3MonTumDensStdDose...
        mean(err3MonTumDens(err3MonTumDens(:, 5) == tDose(i), 1))];
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
errEndTreatTumDensTTumDose = zeros(length(tTTum), length(tDose));
for i = 1:length(tTTum)
    for j = 1:length(tDose)
        errEndTreatTumDensTTumDose(i, j) =...
            mean(errEndTreatTumDens(err3MonTumDens(:, 3) == tTTum(i) &...
            err3MonTumDens(:, 5) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensTTumDose, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('TTum (h)')
ylabel('Dose (Gy)')
xticks(1:length(tTTum))
xticklabels(tTTum)
yticks(1:length(tDose))
yticklabels(tDose)

%%
errEndTreatTumDensDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        errEndTreatTumDensDThresDose(i, j) =...
            mean(errEndTreatTumDens(err3MonTumDens(:, 4) == tDThres(i) &...
            err3MonTumDens(:, 5) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensDThresDose, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('Dthres (Gy)')
ylabel('Dose (Gy)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tDose))
yticklabels(tDose)

%%
errEndTreatTumDensTTumDThres = zeros(length(tTTum), length(tDThres));
for i = 1:length(tTTum)
    for j = 1:length(tDThres)
        errEndTreatTumDensTTumDThres(i, j) =...
            mean(errEndTreatTumDens(err3MonTumDens(:, 3) == tTTum(i) &...
            err3MonTumDens(:, 4) == tDThres(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errEndTreatTumDensTTumDThres, 'CDataMapping','scaled')
colorbar
title('Abs. error in tumor density at the end of treat.')
xlabel('TTum (h)')
ylabel('DThres (Gy)')
xticks(1:length(tTTum))
xticklabels(tTTum)
yticks(1:length(tDThres))
yticklabels(tDThres)

%%
errIntTumDensDThresDose = zeros(length(tDThres), length(tDose));
for i = 1:length(tDThres)
    for j = 1:length(tDose)
        errIntTumDensDThresDose(i, j) =...
            mean(errIntTumDens(err3MonTumDens(:, 4) == tDThres(i) &...
            err3MonTumDens(:, 5) == tDose(j), 1));
    end
end

nfig = nfig + 1;
figure(nfig)
image(errIntTumDensDThresDose, 'CDataMapping','scaled')
colorbar
title('Abs. error in integral of tumor density')
xlabel('Dthres (Gy)')
ylabel('Dose (Gy)')
xticks(1:length(tDThres))
xticklabels(tDThres)
yticks(1:length(tDose))
yticklabels(tDose)
