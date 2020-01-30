clear all
close all


nfig = 0;
path = '../../Carlos/Results/180x180_0.65_0.25_0.04_1_2x37Gy_NoHypNec/';
endTreat = load([path, 'endTreatTumDens.res']);

blue   = [46/255, 165/255, 225/255];
red    = [192/255, 0, 0];
redTum = [153/255, 0, 51/255];
green  = [156/255, 204/255, 88/255];
orange = [246/255, 161/255, 27/255];
violet = [109/255, 96/255, 214/255];
redG0  = [187/255, 81/255, 51/255];

%%
tumDens = load([path, 'tumDens.res'])';

% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));

hold on
plot(tumDens(1, :), tumDens(2, :), 'linewidth', 6, 'color', blue)
plot([endTreat(1), endTreat(1)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(endTreat(1), 1.05 * max(tumDens(2, :)),'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([tumDens(1, end), tumDens(1, end)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(tumDens(1, end), 1.05 * max(tumDens(2, :)),'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
title('Tumor density', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Tumor density (%)', 'fontsize', 56)
xlim([0, inf])
ylim([0, inf])
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
title('Distribution of tumor cells', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Tumor cells in each phase (%)', 'fontsize', 50)
legend({'G1', 'S', 'G2', 'M', 'G0'}, 'location', 'northwest',...
    'orientation', 'horizontal','fontsize', 56)
xlim([0, inf])
ylim([0, inf])
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
plot(pO2(1, 2:end), pO2(3, 2:end), 'linewidth', 6, 'color', green)
plot([endTreat(1), endTreat(1)], [0, 30], '--k', 'Linewidth', 6)
text(endTreat(1), 28,'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([pO2(1, end), pO2(1, end)], [0, 30], '--k', 'Linewidth', 6)
text(pO2(1, end), 28,'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
title('pO2', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('pO2 (mmHg)', 'fontsize', 56)
legend({'Median', 'Mean'}, 'location', 'southwest', 'fontsize', 56)
xlim([0, inf])
ylim([0, inf])
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
plot(vegf(1, 2:end), vegf(3, 2:end), 'linewidth', 6, 'color', green)
plot([endTreat(1), endTreat(1)], [0, 15], '--k', 'Linewidth', 6)
text(endTreat(1), 14, '8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([vegf(1, end), vegf(1, end)], [0, 15], '--k', 'Linewidth', 6)
text(vegf(1, end), 14,'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
title('VEGF concentration', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('v (mol/μm²)', 'fontsize', 56)
legend({'Median', 'Mean'}, 'location', 'southwest', 'fontsize', 56)
xlim([0, inf])
ylim([0, inf])
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
path = '../../Carlos/Results/40x2Gy_NoHypNec_TSim5weeks/';
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
for i = 1:nTissues
    tumDens(:, :, i) = load([pathTissue, num2str(i), '/tumDens_0.res'])';
    %     cycle(:, :, i) = load([pathTissue, num2str(i), '/cycle_0.res'])';
    %     pO2Stat(:, :, i) = load([pathTissue, num2str(i), '/pO2Stat_0.res'])';
    %     vegfStat(:, :, i) = load([pathTissue, num2str(i), '/vegfStat_0.res'])';
end

meanTumDens = mean(tumDens, 3);
stdTumDens = std(tumDens, 0, 3);

%%
% nfig = nfig + 1;
% figure(nfig)
fig = figure('position', get(0, 'screensize'));
hold on
plot(permute(tumDens(1, :, :),[2, 3, 1]),...
    permute(tumDens(2, :, :),[2, 3, 1]), 'linewidth', 6)
plot([endTreat(1), endTreat(1)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(endTreat(1), 1.05 * max(tumDens(2, :)),'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([threeMon(1), threeMon(1)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(threeMon(1), 1.05 * max(tumDens(2, :)),'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
hold off
title('Tumor density', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Tumor density (%)', 'fontsize', 56)
xlim([0, inf])
ylim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
% print(fig, 'tumDens.pdf', '-dpdf', '-r0')
saveas(fig, 'tumDens.svg')

fig = figure('position', get(0, 'screensize'));
hold on
plot(meanTumDens(1, :), meanTumDens(2, :), 'Linewidth', 6)
% fill(xconf, yconf, 1,'facecolor', 'red', 'edgecolor', 'none',...
%     'facealpha', 0.4);
plot([endTreat(1), endTreat(1)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(endTreat(1), 1.05 * max(tumDens(2, :)),'8 w ','fontsize', 56,...
    'horizontalalignment', 'right')
plot([threeMon(1), threeMon(1)], [0, 1.1 * max(tumDens(2, :))],...
    '--k', 'Linewidth', 6)
text(threeMon(1), 1.05 * max(tumDens(2, :)),'12 w ','fontsize', 56,...
    'horizontalalignment', 'right')
hold off
title('Tumor density', 'fontsize', 56)
xlabel('Time (h)', 'fontsize', 56)
ylabel('Tumor density (%)', 'fontsize', 56)
xlim([0, inf])
ylim([0, inf])
ax = gca;
ax.FontSize = 56;
grid on
set(fig, 'units', 'inches')
pos = get(fig, 'position');
set(fig, 'paperPositionMode', 'auto', 'paperUnits', 'Inches',...
    'paperSize', [pos(3), pos(4)])
% print(fig, 'tumDens.pdf', '-dpdf', '-r0')
saveas(fig, 'meanTumDens.svg')

