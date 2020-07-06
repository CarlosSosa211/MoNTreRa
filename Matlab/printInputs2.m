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

tissuePar = load('../InputFiles/tumAreaDensADC.dat');
totalDose = load('../InputFiles/totalDose.dat');

nTissues = size(tissuePar, 1);

for i = 1: nTissues
    fid = fopen(['../Khemara/tissuePar_', num2str(i),...
        '.dat'], 'w');
    fprintf(fid, '20.0\n%.1f\n%.2f\n0.05', tissuePar(i, 1) * 1e4,...
        tissuePar(i, 2));
    fclose(fid);
    
        fid = fopen(['../Khemara/treatment_', num2str(i),...
        '.dat'], 'w');
    fprintf(fid, '2.0\n%.1f\n24.0\n0', totalDose(i));
    fclose(fid);
end



