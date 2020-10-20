clear all
close all

global nfig
nfig = 0;

% path = '../../Carlos/Results/Var_TTum_2Gy_noHypNec';
% path = '../../Carlos/Results/Var_d_noHypNec';
% path = '../../Carlos/Results/Var_ang';
path = '../OutputFiles';
p = 2;

%%
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

%%
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

%%
nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    vascDens = load([path, '/vascDens_', num2str(i), '.res']);
    plot(vascDens(:, 1), vascDens(:, 2));
    plot(vascDens(:, 1), vascDens(:, 3));
    plot(vascDens(:, 1), vascDens(:, 4));
end

hold off
title('Vascular density')
xlabel('Time (h)')
ylabel('Vascular density')
legend
grid on

%%
nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    cycle = load([path, '/cycle_', num2str(i), '.res']);
    plot(cycle(:, 1), cycle(:, 2));
end

hold off
title('Evolution of the distribution of tumor cells')
xlabel('Time (h)')
ylabel('Percentage of tumor cells in phase G1')
legend
grid on

%%
nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    cycle = load([path, '/cycle_', num2str(i), '.res']);
    plot(cycle(:, 1), cycle(:, 3));
end

hold off
title('Evolution of the distribution of tumor cells')
xlabel('Time (h)')
ylabel('Percentage of tumor cells in phase S')
legend
grid on

%%
nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    cycle = load([path, '/cycle_', num2str(i), '.res']);
    plot(cycle(:, 1), cycle(:, 6));
end

hold off
title('Evolution of the distribution of tumor cells')
xlabel('Time (h)')
ylabel('Percentage of tumor cells in phase G0')
legend
grid on

%%
nfig = nfig + 1;
figure(nfig)
hold on

for i = 0:p-1
    hypDens = load([path, '/hypDens_', num2str(i), '.res']);
    plot(hypDens(:, 1), hypDens(:, 2));
end

hold off
title('Hypoxic cells density')
xlabel('Time (h)')
ylabel('Hypoxic cells density')
legend
grid on


