clear all
close all
nfig = 0;
img = imread('../../Carlos/Figures/Histo/CD31_P14-1685_Bloc8-Zoom.jpg');
[nrow, ncol] = size(img(:, :, 1));
%%
nfig = nfig + 1;
figure(nfig)
imshow(img)

%%
%img = imresize(img, 0.1);
%nfig = nfig + 1;
%figure(nfig)
%imshow(img)

%%
%imgR= img(:, :, 1);
%nfig = nfig + 1;
%figure(nfig)
%imshow(imgR)

%%
%nfig = nfig + 1;
%figure(nfig)
%imshow(imgR < 80)

%%
imgHSV= rgb2hsv(img);
%nfig = nfig + 1;
%figure(nfig)
%imshow(imgHSV)

%%
imgH = imgHSV(:, :, 1);
%nfig = nfig + 1;
%figure(nfig)
%imshow(imgH)

%%
imgEnd = imgH < 0.08;
nfig = nfig + 1;
figure(nfig)
imshow(imgEnd)

SE = strel('square',2); 
imgEnd = imopen(imgEnd,SE);
imgEnd = imclose(imgEnd,SE);

nfig = nfig + 1;
figure(nfig)
imshow(imgEnd)

%%
imgEnd = imresize(imgEnd, 0.1);
nfig = nfig + 1;
figure(nfig)
imshow(imgEnd)

fid = fopen('inVes0.dat','wt');
for i = 1:size(imgEnd, 1)
    fprintf(fid,'%i\t', imgEnd(i,:));
    fprintf(fid,'\n');
end
fclose(fid)

imgTum = rand(size(imgEnd));
imgTum = imgTum > .32 & ~imgEnd;
fid = fopen('inTum0.dat','wt');
for i = 1:size(imgTum, 1)
    fprintf(fid,'%i\t', imgTum(i,:));
    fprintf(fid,'\n');
end
fclose(fid)

fid = fopen('tissueDim0.dat','wt');
fprintf(fid,'%d\n', size(imgEnd));
fprintf(fid,'%d\n',1);
fprintf(fid,'%0.1f\n', 20.0);
fclose(fid)


