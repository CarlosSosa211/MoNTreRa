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
