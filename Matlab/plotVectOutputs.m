clear all
close all


nfig = 0;
%%
path = '../../Carlos/Results/180x180_0.65_0.25_0.04_1_2x37Gy_NoHypNec/';
% path = '../../Carlos/Results/40x2GyNoHypNec/';
endTreat = load([path, 'endTreatTumDens.res']);
threeMon = load([path, '3MonTumDens.res']);

%%
blue      = [46/255, 165/255, 225/255];
transBlue = [181/255, 210/255, 225/255];
red       = [255/255, 124/255, 128/255];
redTum    = [153/255, 0, 51/255];
green     = [153/255, 255/255, 102/255];
orange    = [246/255, 161/255, 27/255];
violet    = [109/255, 96/255, 214/255];
redG0     = [187/255, 81/255, 51/255];

%%
fig = figure('position', get(0, 'screensize'));

hold on
tumDens = load('../../Carlos/Results/37x2Gy/tumDens.res');
plot(tumDens(:, 1), tumDens(:, 2), 'linewidth', 6, 'color', blue)

tumDens = load('../../Carlos/Results/20x3Gy/tumDens.res');
plot(tumDens(:, 1), tumDens(:, 2), 'linewidth', 6, 'color', green)

plot([1248, 1248], [0, 1.1 * max(tumDens(:, 2))], '--k', 'Linewidth', 6)
text(1248, 1.25 * max(tumDens(2, :)), sprintf('End\n37 x 2 Gy'),...
    'fontsize', 40, 'horizontalalignment', 'center')

plot([696, 696], [0, 1.1 * max(tumDens(:, 2))], '--k', 'Linewidth', 6)
text(696, 1.25 * max(tumDens(2, :)), sprintf('End \n20 x 3 Gy'),...
    'fontsize', 40, 'horizontalalignment', 'center')


xlabel('Time (h)', 'fontsize', 56)
ylabel('Tumor density (%)', 'fontsize', 56)
xlim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
hold off

legend({'37 x 2 Gy', '20 x 3 Gy'}, 'fontsize', 40)

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])

%%
close all
fig = figure('position', get(0, 'screensize'));
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w_20x3/';
hold on

rec = [10, 11, 19, 40, 42, 48, 63, 70, 73];
nRec = 0;
nNoRec = 0;
for i = 1:76
    if any(rec == i)
        nRec = nRec +1;
        for k = 0:4
            tumVol(:, :, k + 1) = load([path, 'rep', num2str(k), '/tumVol_', num2str(i), '.res']);
        end
        tumVol(:, 2, :) = tumVol(:, 2, :) ./ tumVol(1, 2, :);
        tumVolRec(:, :, nRec) = mean(tumVol, 3);
        
    else
        nNoRec = nNoRec +1;
        for k = 0:4
            tumVol(:, :, k + 1) = load([path, 'rep', num2str(k), '/tumVol_', num2str(i), '.res']);
        end
        tumVol(:, 2, :) = tumVol(:, 2, :) ./ tumVol(1, 2, :);
        tumVolNoRec(:, :, nNoRec) = mean(tumVol, 3);
    end
end

meanTumVolRec = mean(tumVolRec, 3);
stdTumVolRec = std(tumVolRec, 0, 3);
meanTumVolNoRec = mean(tumVolNoRec, 3);
stdTumVolNoRec = std(tumVolNoRec, 0, 3);

plot(meanTumVolRec(:, 1), meanTumVolRec(:, 2), 'linewidth', 6, 'color', red)
plot(meanTumVolNoRec(:, 1), meanTumVolNoRec(:, 2), 'linewidth', 6, 'color', green)

x = [meanTumVolRec(:, 1)', fliplr(meanTumVolRec(:, 1)')];
y = [meanTumVolRec(:, 2)' - stdTumVolRec(:, 2)', fliplr(meanTumVolRec(:, 2)' + stdTumVolRec(:, 2)')];
fill(x, y, red, 'edgecolor', 'none', 'FaceAlpha', 0.1);

x = [meanTumVolNoRec(:, 1)', fliplr(meanTumVolNoRec(:, 1)')];
y = [meanTumVolNoRec(:, 2)' - stdTumVolNoRec(:, 2)', fliplr(meanTumVolNoRec(:, 2)' + stdTumVolNoRec(:, 2)')];
fill(x, y, green, 'edgecolor', 'none', 'FaceAlpha', 0.1);


plot([1400, 1400], [0, 1.1], '--k', 'Linewidth', 6)
text(1400, 1.05, '8 weeks ','fontsize', 44, 'horizontalalignment', 'right')
plot([2160, 2160], [0, 1.1], '--k', 'Linewidth', 6)
text(2160, 1.05, '12 weeks ','fontsize', 44, 'horizontalalignment', 'right')
xlabel('Time (h)', 'fontsize', 44)
ylabel('Normalized n° of tumor cells ', 'fontsize', 44)
xlim([0, inf])
ylim([0, 1.1])
yticks(0:0.2:1)
ax = gca;
ax.FontSize = 44;
grid on
hold off

legend({'Recurrence', 'No recurrence'},...
    'fontsize', 40, 'location', 'southwest')

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])

%%
close all
fig = figure('position', get(0, 'screensize'));
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w/';
hold on

rec = [10, 11, 19, 40, 42, 48, 63, 70, 73];
nRec = 0;
nNoRec = 0;
for i = 1:76
    if any(rec == i)
        nRec = nRec +1;
        for k = 0:4
            tumVol(:, :, k + 1) = load([path, 'rep', num2str(k), '/tumVol_', num2str(i), '.res']);
        end
        tumVol(:, 2, :) = tumVol(:, 2, :) ./ tumVol(1, 2, :);
        tumVolRec(:, :, nRec) = mean(tumVol, 3);
        
    elseif(i ~= 27)
        nNoRec = nNoRec + 1;
        for k = 0:4
            tumVol(:, :, k + 1) = load([path, 'rep', num2str(k), '/tumVol_', num2str(i), '.res']);
        end
        tumVol(:, 2, :) = tumVol(:, 2, :) ./ tumVol(1, 2, :);
        tumVolNoRec(:, :, nNoRec) = mean(tumVol, 3);
    end
end

meanTumVolRec = mean(tumVolRec, 3);
stdTumVolRec = std(tumVolRec, 0, 3);
meanTumVolNoRec = mean(tumVolNoRec, 3);
stdTumVolNoRec = std(tumVolNoRec, 0, 3);

meanTumVolRec(:, 1) = meanTumVolRec(:, 1)/ (24 * 7);
meanTumVolNoRec(:, 1) = meanTumVolNoRec(:, 1)/ (24 * 7);
% plot(meanTumVolRec(:, 1), meanTumVolRec(:, 2), 'linewidth', 6, 'color', red)
% plot(meanTumVolNoRec(:, 1), meanTumVolNoRec(:, 2), 'linewidth', 6, 'color', green)

% x = [meanTumVolRec(:, 1)', fliplr(meanTumVolRec(:, 1)')];
% y = [meanTumVolRec(:, 2)' - stdTumVolRec(:, 2)', fliplr(meanTumVolRec(:, 2)' + stdTumVolRec(:, 2)')];
% fill(x, y, red, 'edgecolor', 'none', 'FaceAlpha', 0.1);
%
% x = [meanTumVolNoRec(:, 1)', fliplr(meanTumVolNoRec(:, 1)')];
% y = [meanTumVolNoRec(:, 2)' - stdTumVolNoRec(:, 2)', fliplr(meanTumVolNoRec(:, 2)' + stdTumVolNoRec(:, 2)')];
% fill(x, y, green, 'edgecolor', 'none', 'FaceAlpha', 0.1);

plot(meanTumVolRec(1:225, 1), meanTumVolRec(1:225, 2), 'linewidth', 6, 'color', red)
plot(meanTumVolNoRec(1:225, 1), meanTumVolNoRec(1:225, 2), 'linewidth', 6, 'color', green)

x = [meanTumVolRec(1:225, 1)', fliplr(meanTumVolRec(1:225, 1)')];
y = [meanTumVolRec(1:225, 2)' - stdTumVolRec(1:225, 2)', fliplr(meanTumVolRec(1:225, 2)' + stdTumVolRec(1:225, 2)')];
fill(x, y, red, 'edgecolor', 'none', 'FaceAlpha', 0.3);

x = [meanTumVolNoRec(1:225, 1)', fliplr(meanTumVolNoRec(1:225, 1)')];
y = [meanTumVolNoRec(1:225, 2)' - stdTumVolNoRec(1:225, 2)', fliplr(meanTumVolNoRec(1:225, 2)' + stdTumVolNoRec(1:225, 2)')];
fill(x, y, green, 'edgecolor', 'none', 'FaceAlpha', 0.3);


% plot([1400, 1400], [0, 1.1], '--k', 'Linewidth', 6)
% text(1400, 1.05, '8 weeks ','fontsize', 44, 'horizontalalignment', 'right')
% plot([2160, 2160], [0, 1.1], '--k', 'Linewidth', 6)
% text(2160, 1.05, '12 weeks ','fontsize', 44, 'horizontalalignment', 'right')
xlabel('Time (weeks)', 'fontsize', 44)
ylabel('Normalized n° of tumor cells ', 'fontsize', 44)
xlim([0, inf])
ylim([0, 1.1])
yticks(0:0.2:1)
ax = gca;
ax.FontSize = 44;
grid on
hold off

legend({'Recurrence', 'No recurrence'},...
    'fontsize', 40, 'location', 'northeast')

outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, 'nTumCells.pdf', '-dpdf', '-r0')

%%
eightwTumVolRec = reshape(tumVolRec(225, 2, :), [9, 1])';
eightwTumVolNoRec = reshape(tumVolNoRec(225, 2, :), [66, 1])';

fid = fopen('../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w/8wTumVolNormRec.res', 'w');
for i = 1 : length(eightwTumVolRec)
    fprintf(fid, '%.3f\n', eightwTumVolRec(i));
end
fclose(fid);

fid = fopen('../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w/8wTumVolNormNoRec.res', 'w');
for i = 1 : length(eightwTumVolNoRec)
    fprintf(fid, '%.3f\n', eightwTumVolNoRec(i));
end
fclose(fid);

%%
% nfig = nfig + 1;
% figure(nfig)
pos = get(0, 'screensize');
pos(4) = 0.35 * pos(4);
fig = figure('position', pos);
treatment = 24:24:endTreat(1) - 24;
for k = 1:length(treatment)
    if mod(k, 7) == 6 || mod(k, 7) == 0
        treatment(k) = 0;
    end
end
treatment = treatment(treatment ~= 0);
hold on
stem(treatment, 0.4 * ones(1, length(treatment)), '-v', 'lineWidth', 6,...
    'markerSize', 6, 'markerFaceColor', blue)
text(0.5 * endTreat(1), 0.75, '40 x 2 Gy','fontsize', 56,...
    'horizontalalignment', 'center')
plot([endTreat(1), endTreat(1)], [0, 1],...
    '--k', 'Linewidth', 6)
text(endTreat(1), 0.75, '8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([threeMon(1), threeMon(1)], [0, 1], '--k', 'Linewidth', 6)
text(threeMon(1), 0.75, '12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
% title('Tumor density', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
xlim([0, inf])
ylim([0, 1])
ax = gca;
set(ax,'ytick',[])
ax.FontSize = 56;
grid on
hold off

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
% print(fig, 'tumDens.pdf', '-dpdf', '-r0')
saveas(fig, 'treatment.svg')
%%
tumDens = load([path, 'tumDens.res'])';
tumDens(2, :) = tumDens(2, :) * 180 * 180 * 4e-6; % 20 * 20 * 10^-6 /100

% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));

hold on
plot(tumDens(1, :), tumDens(2, :), 'linewidth', 6, 'color', blue)
area(tumDens(1, :), tumDens(2, :), 'lineStyle', 'none', 'faceColor', transBlue)
text(200, 2,'Integral', 'fontsize', 56,...
    'horizontalalignment', 'left')
plot([endTreat(1), endTreat(1)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(endTreat(1), 1.05 * max(tumDens(2, :)),'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([tumDens(1, end), tumDens(1, end)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(tumDens(1, end), 1.05 * max(tumDens(2, :)),'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
% title('Tumor density', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
% ylabel('Tumor density (%)', 'fontsize', 56)
ylabel('Tumor area (mm^2)', 'fontsize', 56)
xlim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
hold off

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
% print(fig, 'tumDens.pdf', '-dpdf', '-r0')
saveas(fig, 'tumDens.svg')

%%
cycle = load([path, 'cycle.res'])';

% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));

hold on
plot(cycle(1, :), cycle(2, :), 'linewidth', 6, 'color', blue)
plot(cycle(1, :), cycle(3, :), 'linewidth', 6, 'color', green)
plot(cycle(1, :), cycle(4, :), 'linewidth', 6, 'color', orange)
plot(cycle(1, :), cycle(5, :), 'linewidth', 6, 'color', violet)
plot(cycle(1, :), cycle(6, :), 'linewidth', 6, 'color', redG0)
plot([endTreat(1), endTreat(1)], [0, 100], '--k', 'Linewidth', 6)
text(endTreat(1), 95, '8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([cycle(1, end), cycle(1, end)], [0, 100], '--k', 'Linewidth', 6)
text(cycle(1, end), 95,'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
ax = gca;
ax.FontSize = 56;
% title('Distribution of tumor cells', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Tumor cells in each phase (%)', 'fontsize', 50)
% legend({'G1', 'S', 'G2', 'M', 'G0'}, 'location', 'northwest',...
%     'orientation', 'horizontal','fontsize', 56)
xlim([0, inf])
grid on
hold off

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
% print(fig, 'cycle.pdf', '-dpdf', '-r0')
saveas(fig, 'cycle.svg')
%%
pO2 = load([path, 'pO2Stat.res'])';

% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));

hold on
plot(pO2(1, 2:end), pO2(2, 2:end), 'linewidth', 6, 'color', blue)
% plot(pO2(1, 2:end), pO2(3, 2:end), 'linewidth', 6, 'color', green)
plot([endTreat(1), endTreat(1)], [0, 30], '--k', 'Linewidth', 6)
text(endTreat(1), 28,'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([pO2(1, end), pO2(1, end)], [0, 30], '--k', 'Linewidth', 6)
text(pO2(1, end), 28,'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
% title('pO2', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Median pO₂ (mmHg)', 'fontsize', 56)
% legend({'Median', 'Mean'}, 'location', 'southwest', 'fontsize', 56)
xlim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
hold off

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
% print(fig, 'pO2.pdf', '-dpdf', '-r0')
saveas(fig, 'pO2.svg')
%%
vegf = load([path, 'vegfStat.res'])';

% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));

hold on
plot(vegf(1, 2:end), vegf(2, 2:end), 'linewidth', 6, 'color', blue)
% plot(vegf(1, 2:end), vegf(3, 2:end), 'linewidth', 6, 'color', green)
plot([endTreat(1), endTreat(1)], [0, 15], '--k', 'Linewidth', 6)
text(endTreat(1), 14, '8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([vegf(1, end), vegf(1, end)], [0, 15], '--k', 'Linewidth', 6)
text(vegf(1, end), 14,'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
% title('VEGF concentration', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Median v (mol/μm²)', 'fontsize', 56)
% legend({'Median', 'Mean'}, 'location', 'southwest', 'fontsize', 56)
xlim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
hold off

set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, 'VEGF.pdf', '-dpdf', '-r0')
saveas(fig, 'VEGF.svg')

%%
clear all
close all

nfig = 0;
path = '../../Carlos/Results/40x2Gy_NoHypNec_TSim5months';
pathTissue = [path, '/Tissue'];
endTreat = load([pathTissue, '1/endTreatTumDens_0.res']);
threeMon = load([pathTissue, '1/3MonTumDens_0.res']);

blue   = [46/255, 165/255, 225/255];
red    = [192/255, 0, 0];
redTum = [153/255, 0, 51/255];
green  = [156/255, 204/255, 88/255];
orange = [246/255, 161/255, 27/255];
violet = [109/255, 96/255, 214/255];
redG0  = [187/255, 81/255, 51/255];

nTissues = 21;
for k = 1:nTissues
    output(:, :, k) = load([pathTissue, num2str(k), '/tumDens_0.res'])';
    %     output(:, :, i) = load([pathTissue, num2str(i), '/tumVol_0.res'])';
    %     output(:, :, i) = load([pathTissue, num2str(i), '/cycle_0.res'])';
    %     output(:, :, i) = load([pathTissue, num2str(i), '/pO2Stat_0.res'])';
    %     output(:, :, i) = load([pathTissue, num2str(i), '/vegfStat_0.res'])';
end

outputName = 'Tumor density';
outputUnit = '(%)';
% outputName = 'Tumor volume';
% outputUnit = '(mm³)';
meanOutput = mean(output, 3);
stdOutput = std(output, 0, 3);

%%
% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));
hold on
plot(permute(output(1, :, :),[2, 3, 1])/24,...
    permute(output(2, :, :),[2, 3, 1]), 'linewidth', 6)
plot([endTreat(1)/24, endTreat(1)/24], [0, 1.1 * max(output(2, :))],...
    '--k', 'Linewidth', 6)
text(endTreat(1)/24, 1.05 * max(output(2, :)),'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([threeMon(1)/24, threeMon(1)/24], [0, 1.1 * max(output(2, :))],...
    '--k', 'Linewidth', 6)
text(threeMon(1)/24, 1.05 * max(output(2, :)),'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
hold off
title(outputName, 'fontsize', 56)
xlabel('Time (days)', 'fontsize', 56)
% xlabel('Time (h)', 'fontsize', 56)
ylabel([outputName, ' ', outputUnit], 'fontsize', 56)
xlim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, 'tumDens.pdf', '-dpdf', '-r0')

fig = figure('position', get(0, 'screensize'));
for k = 1:nTissues
    subplot(3, 7, k)
    hold on
    plot(output(1, :, k)/24, output(2, :, k), 'linewidth', 4)
    plot([endTreat(1)/24, endTreat(1)/24], [0, 1.1 * max(output(2, :, k))],...
        '--k', 'Linewidth', 4)
    text(endTreat(1)/24, 1.05 * max(output(2, :, k)),'8 w ','fontsize', 16,...
        'horizontalalignment', 'right')
    plot([threeMon(1)/24, threeMon(1)/24], [0, 1.1 * max(output(2, :, k))],...
        '--k', 'Linewidth', 4)
    text(threeMon(1)/24, 1.05 * max(output(2, :, k)),'12 w ','fontsize', 16,...
        'horizontalalignment', 'left')
    hold off
    title(outputName, 'fontsize', 16)
    xlabel('Time (days)', 'fontsize', 16)
    % xlabel('Time (h)', 'fontsize', 16)
    ylabel([outputName, ' ',  outputUnit], 'fontsize', 16)
    xlim([0, inf])
    ax = gca;
    ax.FontSize = 16;
    grid on
end
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, 'subTumDens.pdf', '-dpdf', '-r0')

fig = figure('position', get(0, 'screensize'));
hold on
plot(meanOutput(1, :)/24, meanOutput(2, :), 'Linewidth', 6)
plot([endTreat(1)/24, endTreat(1)/24], [0, 1.1 * max(meanOutput(2, :))],...
    '--k', 'Linewidth', 6)
text(endTreat(1)/24, 1.05 * max(meanOutput(2, :)),'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([threeMon(1)/24, threeMon(1)/24], [0, 1.1 * max(meanOutput(2, :))],...
    '--k', 'Linewidth', 6)
text(threeMon(1)/24, 1.05 * max(meanOutput(2, :)),'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
hold off
title(outputName, 'fontsize', 56)
xlabel('Time (days)', 'fontsize', 56)
% xlabel('Time (h)', 'fontsize', 56)
ylabel([outputName, ' ',  outputUnit], 'fontsize', 56)
xlim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
print(fig, 'meanTumDens.pdf', '-dpdf', '-r0')

