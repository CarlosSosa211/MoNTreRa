nfig = 0;
img = imread('../../Carlos/Figures/TestHES.png');
[nrow, ncol] = size(img(:, :, 1));
%%
nfig = nfig + 1;
figure(nfig)
imshow(img)

%%
img(:, :, 3) = img(:, :, 1);
nfig = nfig + 1;
figure(nfig)
imshow(imgR)

%%
imgnB = img;
imgnB(:, :, 3) = zeros(nrow, ncol);
nfig = nfig + 1;
figure(nfig)
imshow(imgnB)

%%
imgG = img(:, :, 2);
nfig = nfig + 1;
figure(nfig)
imshow(imgG)


%%
imgB = img(:, :, 3);
nfig = nfig + 1;
figure(nfig)
imshow(imgB)

%%
nfig = nfig + 1;
figure(nfig)
imshow(imgR < 110)

%%
imgHSV= rgb2hsv(img);
nfig = nfig + 1;
figure(nfig)
imshow(imgHSV)

%%
imgH = imgHSV(:, :, 1);
nfig = nfig + 1;
figure(nfig)
imshow(imgH)

%%
nfig = nfig + 1;
figure(nfig)
imshow(imgH < 0.78)