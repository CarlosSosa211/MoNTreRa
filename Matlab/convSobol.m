clear all
close all
nfig = 0;

allTissues = [1, 2, 21, 22, 31, 41];
tN = [2, 4, 8, 16, 32, 64, 128, 256, 512];
log10N = log10(tN);
nTissues = length(allTissues);

path = '../../Carlos/Results/Sobol/1000_2TreatPar';

fileNames = {'8wTumDens', '12wTumDens', '8wTumVol', '12wTumVol' ...
    '8wIntTumDens', '12wIntTumDens', '8wIntTumVol', '12wIntTumVol' ...
    'Killed50', 'Killed80', 'Killed90', 'Killed95', 'Killed99'...
    'Killed999', 'TimeTo95', 'TimeTo99', 'Rec', 'RecTumDens'...
    'RecTumVol', 'RecTime'};
outputNames = {'Tumor density at t = 8 weeks'...
    'Tumor density at t = 12 weeks', 'Tumor area at t = 8 weeks'...
    'Tumor area at t = 12 weeks'...
    'Integral of tumor density up to t = 8 weeks'...
    'Integral of tumor density up to t = 12 weeks'...
    'Integral of tumor area up to t = 8 weeks'...
    'Integral of tumor area up to t = 12 weeks'...
    '50% of tumor cells killed', '80% of tumor cells killed'...
    '90% of tumor cells killed', '95% of tumor cells killed'...
    '99% of tumor cells killed', '99.9% of tumor cells killed'...
    'Time to kill 95% of tumor cells'...
    'Time to kill 99% of tumor cells', 'Recurrence'...
    'Tumor density at recurrence', 'Tumor area at recurrence'...
    'Recurrence time'};
nOut = 20;

selOut = input(['Select an output [8wTumDens (1), '...
    '12wTumDens (2), 8wTumVol (3),\n12wTumVol (4), '...
    '8wIntTumDens (5), 12wIntTumDens (6), 8wIntTumVol (7),\n'...
    '12wInTumVol (8), 50%killed (9), 80%killed (10), '...
    '90%killed (11),\n95%killed (12), 99%killed (13), '...
    '99.9%killed (14), timeTo95 (15), timeTo99 (16),\nrec(17), '...
    'recTumDens (18), recTumVol (19), recTime (20)\n'...
    'or all of them (-1)]: ']);

for i = 1:nTissues
    SI(:, :, i) = load([path, '/convSI', char(fileNames(selOut))...
        '_', num2str(allTissues(i)), '.res']);
    TSI(:, :, i) = load([path, '/convTSI', char(fileNames(selOut))...
        '_', num2str(allTissues(i)), '.res']);
    sol(:, :, i) = load([path, '/sobol', char(fileNames(selOut))...
        '_', num2str(allTissues(i)), '.res']);
end

K = sol(1, 1);
sol = sol(2:K + 1, :, :);

for k = 1:K
    errSI(k, :, :) = abs(SI(k, :, :) - sol(k, 1, :));
    errTSI(k, :, :) = abs(TSI(k, :, :) - sol(k, 2, :));
end

maxErrSI = max(errSI);
maxErrTSI = max(errTSI);

meanErrSI = reshape(mean(log10(maxErrSI), 3), 1, []);
stdErrSI = reshape(std(log10(maxErrSI), 0, 3), 1, []);

meanErrTSI = reshape(mean(log10(maxErrTSI), 3), 1, []);
stdErrTSI = reshape(std(log10(maxErrTSI), 0, 3), 1, []);

nfig = nfig + 1;
figure(nfig)
hold on
plot(log10N, meanErrSI, '-o', 'Linewidth', 6);
errorbar(log10N, meanErrSI, stdErrSI, 'ok', 'Linewidth', 6)
hold off
ax = gca;
ax.FontSize = 56;
title('Error SI')
xlabel('log_{10}(N)')
ylabel('log_{10}(err)')
grid on

nfig = nfig + 1;
figure(nfig)
hold on
plot(log10N, meanErrTSI, '-o', 'Linewidth', 6);
errorbar(log10N, meanErrTSI, stdErrTSI, 'ok', 'Linewidth', 6)
hold off
ax = gca;
ax.FontSize = 56;
title('Error TSI')
xlabel('log_{10}(N)')
ylabel('log_{10}(err)')
grid on
