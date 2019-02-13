clear all
close all
nfig = 0;

% path = '../../Carlos/Results/OxyHist/';
path = '../OutputFilesGUI/';
% pathTissue = [path, '0.11_1.0/'];
pathTissue = path;
% stable = load([pathTissue, 'oxyStable.res']);
% stableImg = stable / 10;
pO2 = load([pathTissue, 'po2.res']);
% vegf = load([pathTissue, 'vegf.res']);
N = size(pO2, 1);

tissueDim = load([pathTissue, 'tissueDim.dat']);
nrow = tissueDim(1);
ncol = tissueDim(2);
nlayer = tissueDim(3);
cellSize = tissueDim(4);

pO2 = reshape(pO2, N, ncol, nrow);
pO2 = permute(pO2, [3, 2, 1]);
% vegf = reshape(vegf, N, nrow, ncol);
% vegf = permute(vegf, [3, 2, 1]);

%%
errPO2 = zeros(1, N - 1);
for i = 1:N-1
    errPO2(i) = norm(pO2(:, :, i) - pO2(:, :, i + 1), 'inf');
    % errVegf(i) = norm(vegf(:, :, i) - vegf(:, :, N - 1), 'inf');
end

log10N = log10(1:N - 1);

nfig = nfig + 1;
figure(nfig)
hold on
plot(log10N(errPO2 ~= 0), log10(errPO2((errPO2 ~= 0))))
% if(stableImg > 0)
%     plot([log10(stableImg), log10(stableImg)], [4, -4], '--')
% end
hold off
title('Error pO2')
xlabel('N')
ylabel('log(err)')
grid on

% nfig = nfig + 1;
% figure(nfig)
% hold on
% plot(log10(1:N - 1), log10(errVegf))
% plot([log10(stableImg), log10(stableImg)], [4, -4], '--')
% hold off
% title('Error VEGF')
% xlabel('N')
% ylabel('log(err)')
% grid on

%%
nfig = nfig + 1;
figure(nfig)
image(pO2(:, :, stableImg))
colorbar
title(['stable pO2 for t = ', num2str(stable), ' ms'])

nfig = nfig + 1;
figure(nfig)
image(pO2(:, :, end))
colorbar
title('pO2 for t = 10 000 ms')

nfig = nfig + 1;
figure(nfig)
image(pO2(:, :, stableImg) - pO2(:, :, end), 'CDataMapping','scaled')
colorbar
title('Error')

%%
nBins = 16;

nfig = nfig + 1;
figure(nfig)
histogram(pO2(:, :, stableImg), nBins, 'normalization', 'probability');
ylim([0, 1])
title(['stable pO2 for t = ', num2str(stable), ' ms'])

nfig = nfig + 1;
figure(nfig)
histogram(pO2(:, :, end), nBins, 'normalization', 'probability');
ylim([0, 1])
title('pO2 for t = 10 000 ms')