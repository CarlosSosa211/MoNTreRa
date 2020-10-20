clear all
close all

warning('off')

pkg load image

tTissueDim = [250, 1000, 1];

initNum = input('Enter initial number of tumor cells: ');
nrow = tTissueDim(1);
ncol = tTissueDim(2);
nlayer = tTissueDim(3);

fid = fopen('tissueDim.dat', 'w');
fprintf(fid, '%i\n%i\n%i\n', nrow, ncol, nlayer);
fclose(fid);

initVes = zeros(nrow, ncol, nlayer);
fid = fopen('inVes.dat', 'w');
for k = 1:nlayer
    for i = 1:nrow
        for j = 1:ncol
            fprintf(fid, '%f\n', initVes(i, j, k));
        end
    end
end
fclose(fid);


initTum = zeros(nrow, ncol, nlayer);
count = 0;
flag = 0;
for k = 1:nlayer
    for i = floor(nrow/2):nrow
        for j = 1:ncol
            initTum(i, j, k) = 1;
            count = count + 1;
            if (count == initNum)
                flag = 1;
                break
            end
        end
        if(flag == 1)
            break
        end
    end
    if(flag == 1)
        break
    end
end

fid = fopen('inTum.dat', 'w');
for k = 1:nlayer
    for i = 1:nrow
        for j = 1:ncol
            fprintf(fid, '%f\n', initTum(i, j, k));
        end
    end
end
fclose(fid);

%%
tissuePar = load('../InputFiles/tumAreaDensADCT2w.dat');
% tissuePar = load('../InputFiles/dataOSSyn.csv');
% totalDose = load('../InputFiles/totalDose.dat');

nTissues = size(tissuePar, 1);

for i = 1: nTissues
%     fid = fopen(['../Recurrence/tissueParADCT2w_', num2str(i),...
%         '.dat'], 'w');
%     fprintf(fid, '20.0\n%.1f\n%.2f\n0.038', tissuePar(i, 1) * 1e4,...
%         tissuePar(i, 2));
%     fclose(fid);
    
%     fid = fopen(['../Recurrence/treatment20x3_', num2str(i),...
%         '.dat'], 'w');
%     fprintf(fid, '3.0\n60.0\n24.0\n0');
%         fprintf(fid, '2.0\n%.1f\n24.0\n0', totalDose(i));
    fclose(fid);
end

%%
close all
path = '../../Carlos/Results/Recurrence/simp/TTum330_alphaG1120_ADCT2w/';

fOutput = fopen([path, 'output.res'], 'w');

nTissues = 11600;

for i = 1:nTissues
    meanOutput = 0;
    for j = 0:4
        output = load([path, 'rep', num2str(j), '/8wumVol_'...
            num2str(i), '.res']);
        %    output(:, 2) = filter(b, a, output(:, 2));
        %    plot(output(:, 1), output(:, 2))
        %    line = find(output(:, 1) == 2160);
        %    fprintf(fOutput, '%f\n', output(line, 2));
        %     output = trapz(output(1:line, 1), output(1:line, 2));
%         output = output(1, 2);
        meanOutput = meanOutput + output;
    end
    meanOutput = meanOutput / 5;
    fprintf(fOutput, '%f\n', meanOutput);
end

fclose(fOutput);

%%
close all
path = '../../Carlos/Results/Recurrence/vascDensNoPref0.03_simp_ADC/';

nTissues = 76;

for i = 1:nTissues
        output = load([path, 'rep0', '/tumDens_', num2str(i), '.res']);
        fOutput = fopen([path, 'rep0/8wTumDens_', num2str(i), '.res'], 'w');
        line = find(output(:, 1) == 1440);
        fprintf(fOutput, '%f', output(line, 2));
        fclose(fOutput);
end

% for i = 1:nTissues
%         output = load([path, 'rep0', '/tumDens_', num2str(i), '.res']);
%         fOutput = fopen([path, 'rep0/12wIntTumDens_', num2str(i), '.res'], 'w');
%         line = find(output(:, 1) == 2160);
%         output = trapz(output(1:line, 1), output(1:line, 2));
%         fprintf(fOutput, '%.3f', output);
%         fclose(fOutput);
% end
