clear all 
close all

global nfig
nfig = 0;

path = '../../Carlos/Results/Var_TTum';
p = 10;

nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    tumDens = load([path, '/tumDens_', num2str(i), '.res']);
    plot(tumDens(:, 1), tumDens(:, 2));
end


hold off
title('Tumor density')
xlabel('Time (h)')
ylabel('Tumor density')
legend
grid on

nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    tumVol = load([path, '/tumVol_', num2str(i), '.res']);
    plot(tumVol(:, 1), tumVol(:, 2));
end


hold off
title('Tumor volume')
xlabel('Time (h)')
ylabel('Tumor volume (mm³)')
legend
grid on


nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    vascDens = load([path, '/vascDens_', num2str(i), '.res']);
    %     killedCells = load([path, '/killedCells_', num2str(i), '.res']);
    %     cycle = load([path, '/cycle_', num2str(i), '.res']);
    %     hypDens = load([path, '/hypDens_', num2str(i), '.res']);
    %     pO2Stat = load([path, '/pO2Stat_', num2str(i), '.res']);
    %     vegfStat = load([path, '/vegfStat_', num2str(i), '.res']);
    plot(vascDens(:, 1), vascDens(:, 2));
end

hold off
title('Vascular density')
xlabel('Time (h)')
ylabel('Vascular density')
legend
grid on








