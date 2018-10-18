clear all
close all
% pO2 = load('../../Carlos/Results/Oxy3D/po2.dat')';
% vegf = load('../../Carlos/Results/Oxy3D/vegf.dat');
% dim = load('../../Carlos/Results/Oxy3D/tissueDim.dat');
pO2 = load('../OutputFilesGUI/po2.dat')';
vegf = load('../OutputFilesGUI/vegf.dat');
dim = load('../OutputFilesGUI/tissueDim.dat');
pO2 = reshape(pO2, dim(1), dim(2), []);
pO2 = permute(pO2, [2, 1, 3]);
vegf = reshape(vegf', dim(1), dim(2), []);
vegf = permute(vegf, [2, 1, 3]);
% for i = 1: size(pO2, 3)
%     image(pO2(:, :, i))
% end
% for i = 1: size(vegf, 3)
%     image(vegf(:, :, i))
% end
figure(1)
colormap(jet)
image(pO2(:, :, size(pO2, 3)),'CDataMapping','scaled')
colorbar
figure(2)
colormap(jet)
image(vegf(:, :, size(vegf, 3)),'CDataMapping','scaled')
colorbar
