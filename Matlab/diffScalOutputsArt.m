clear all
close all

%%
% This block initialises nfig and defines
% - path: the path of the ouput directory,
% - par: a matrix containing the combinations of tissue architecture
% parameters simulated. Each row corresponds to a simulation and each
% column, to a parameter,
% - fileNames: the names of the ouput files,
% - ouputNames: the names of the outputs.

nfig = 0;

% path = '../../Carlos/Results/Diff/ArtDens_5Rep/';
% par = load([path, '/combDens.res']);

% par = par * 100;
% tTumDens = unique(par(:, 1));
% tVascDens = unique(par(:, 2));
% parName = {'d_{tum}', 'd_{vasc}'};
% parUnit = {' (%)', ' (%)'};
% parUnitTitle = {'%', '%'};
% parT = {tTumDens, tVascDens};
% 
% selPar1 = 1;
% selPar2 = 2;
% parName1 = char(parName(selPar1));
% parName2 = char(parName(selPar2));
% parUnit1 = char(parUnit(selPar1));
% parUnit2 = char(parUnit(selPar2));
% parT1 = cell2mat(parT(selPar1));
% parT2 = cell2mat(parT(selPar2));


% dtum = 0.9 and dvasc = 0.11 is in fact dtum = 0.89 and dvasc = 0.11
% path = '../../Carlos/Results/Diff/Art_5Rep/';
% par = load([path, '/combDens.res']);

% par(:, 1) = par(:, 1) * 100;
% par(:, 3) = par(:, 3) * 100;
% tTumDens = unique(par(:, 1));
% tSigmaTum = unique(par(:, 2));
% tVascDens = unique(par(:, 3));
% tSigmaVasc = unique(par(:, 4));
% parName = {'d_{tum}', '\sigma_{tum}', 'd_{vasc}', '\sigma_{vasc}'};
% parUnit = {' (%)', '', ' (%)', ''};
% parUnitTitle = {'%', '', '%', ''};
% parT = {tTumDens, tSigmaTum, tVascDens, tSigmaVasc};

% selPar1 = 3;
% selPar2 = 4;
% parName1 = char(parName(selPar1));
% parName2 = char(parName(selPar2));
% parUnit1 = char(parUnit(selPar1));
% parUnit2 = char(parUnit(selPar2));
% parT1 = cell2mat(parT(selPar1));
% parT2 = cell2mat(parT(selPar2));

% selPar1 = 1;
% selPar2 = 2;
% selPar3 = 3;
% selPar4 = 4;
% parName1 = char(parName(selPar1));
% parName2 = char(parName(selPar2));
% parName3 = char(parName(selPar3));
% parName4 = char(parName(selPar4));
% parUnitTitle1 = char(parUnitTitle(selPar1));
% parUnitTitle2 = char(parUnitTitle(selPar2));
% parUnitTitle3 = char(parUnitTitle(selPar3));
% parUnitTitle4 = char(parUnitTitle(selPar4));
% parUnit1 = char(parUnit(selPar1));
% parUnit2 = char(parUnit(selPar2));
% parUnit3 = char(parUnit(selPar3));
% parUnit4 = char(parUnit(selPar4));
% parT1 = cell2mat(parT(selPar1));
% parT2 = cell2mat(parT(selPar2));
% parT3 = cell2mat(parT(selPar3));
% parT4 = cell2mat(parT(selPar4));

path = '../../Carlos/Results/Diff/ArtDensSigmaTum_5Rep/';
par = load([path, '/combDens.res']);

par(:, 1) = par(:, 1) * 100;
par(:, 3) = par(:, 3) * 100;
tTumDens = unique(par(:, 1));
tSigmaTum = unique(par(:, 2));
tVascDens = unique(par(:, 3));
parName = {'d_{tum}', '\sigma_{tum}', 'd_{vasc}'};
parUnit = {' (%)', '', ' (%)'};
parUnitTitle = {'%', '', '%'};
parT = {tTumDens, tSigmaTum, tVascDens};

% selPar1 = 3;
% selPar2 = 4;
% parName1 = char(parName(selPar1));
% parName2 = char(parName(selPar2));
% parUnit1 = char(parUnit(selPar1));
% parUnit2 = char(parUnit(selPar2));
% parT1 = cell2mat(parT(selPar1));
% parT2 = cell2mat(parT(selPar2));

selPar1 = 1;
selPar2 = 2;
selPar3 = 3;
parName1 = char(parName(selPar1));
parName2 = char(parName(selPar2));
parName3 = char(parName(selPar3));
parUnitTitle1 = char(parUnitTitle(selPar1));
parUnitTitle2 = char(parUnitTitle(selPar2));
parUnitTitle3 = char(parUnitTitle(selPar3));
parUnit1 = char(parUnit(selPar1));
parUnit2 = char(parUnit(selPar2));
parUnit3 = char(parUnit(selPar3));
parT1 = cell2mat(parT(selPar1));
parT2 = cell2mat(parT(selPar2));
parT3 = cell2mat(parT(selPar3));
nOut = 15;

fileNames = {'/endTreatTumDens.res', '/3MonTumDens.res'...
    '/finTumVol.res', '/intTumDens.res', '/killed50.res'...
    '/killed80.res', '/killed90.res', '/killed95.res', '/killed99.res'...
    '/killed999.res', '/timeTo95.res', '/timeTo99.res', '/rec.res'...
    '/recTumDens.res', '/recTime.res'};

outputNames = {'Tumour density at the end of treat. (%)'...
    'Tumour density 3 months after the beginning of treat. (%)'...
    'Final tumour volume (mm^3)', 'Integral of tumour density'...
    '50% of tumour cells killed', '80% of tumour cells killed'...
    '90% of tumour cells killed', '95% of tumour cells killed'...
    '99% of tumour cells killed', '99.9% of tumour cells killed'...
    'Time to kill 95% of tumour cells (h)'...
    'Time to kill 99% of tumour cells (h)', 'Recurrence probability'...
    'Tumour density at recurrence (%)', 'Recurrence time (h)'};

%%
% This block asks the user to select the output to be studied. Then,
% it reads the corresponding files. The following intermediate variables
% are considered:
% - output: a matrix containing the values of the scalar output(s) in
% question. For the single output case, (1 <= selOut <= nOut), each row
% corresponds to a combination of parameters. The first column corresponds
% to the mean values for nRep repetitions; the second one, to the the std
% values. For multiple outputs case, (selOut = 0), each layer correponds to
% an output.

selOut = input(['Select an output [endTreaTumDens (1), '...
    '3MonTumDens (2), finTumVol (3),\nintTumDens (4), '...
    '50%killed (5), 80%killed (6), 90%killed (7), 95%killed (8),\n'...
    '99%killed (9), 99.9%killed (10), timeTo95 (11), '...
    'timeTo99 (12), rec(13),\nrecTumDens (14), recTime (15) '...
    'or all of them (-1)]: ']);

if(selOut >= 1 && selOut <= nOut)
    output = load([path, char(fileNames(selOut))]);
    
elseif(selOut == -1)
    for i = 1:length(fileNames)
        output(:, :, i) = load([path, char(fileNames(i))]);
    end
end

%%
% This block plots the the selected output as a function of the studied
% values of parameters par1 and par2
output12 = zeros(length(parT1), length(parT2));
for i = 1:length(parT1)
    for j = 1:length(parT2)
        output12(i, j) =...
            mean(output(par(:, selPar1) == parT1(i) &...
            par(:, selPar2) == parT2(j), 1));
    end
end

% clims = [min(min(output12)), max(max(output12))];
% minColOutput12 = min(output12);
% minOutput12 = min(minColOutput12);
% maxColOutput12 = max(output12);
% maxOutput12 = max(maxColOutput12);
% clims = [min(minColOutput12(minColOutput12 > minOutput12))...
%     max(maxColOutput12(maxColOutput12 < maxOutput12))];
clims = [0, 1];

nfig = nfig + 1;
figure(nfig)
% imagesc(output12)
imagesc(output12, clims)
colorbar
title(char(outputNames(selOut)), 'fontsize', 22)
ax = gca;
ax.FontSize = 22;

xlabel([parName2, parUnit2], 'fontsize', 22)
ylabel([parName1, parUnit1], 'fontsize', 22)
xticks(1:length(parT2))
xticklabels(parT2)
yticks(1:length(parT1))
yticklabels(parT1)

%%
% This block plots the the selected output as a function of the studied
% values of parameters par1, par2 and par3
% clims = [min(output(:, 1)), max(output(:, 1))];
clims = [0, 1];
output12 = zeros(length(parT1), length(parT2));
for k = 1:length(parT3)
        for i = 1:length(parT1)
            for j = 1:length(parT2)
                output12(i, j) =...
                    output(par(:, selPar1) == parT1(i) &...
                    par(:, selPar2) == parT2(j) &... 
                    par(:, selPar3) == parT3(k), 1);
            end
        end
        
        nfig = nfig + 1;
        figure(nfig)
        imagesc(output12, clims)
        colorbar
        title([char(outputNames(selOut)), ' (', char(parName3), ' = '... 
            num2str(parT3(k)), parUnitTitle3, ')'] , 'fontsize', 22)
        ax = gca;
        ax.FontSize = 22;
        
        xlabel([parName2, parUnit2], 'fontsize', 22)
        ylabel([parName1, parUnit1], 'fontsize', 22)
        xticks(1:length(parT2))
        xticklabels(parT2)
        yticks(1:length(parT1))
        yticklabels(parT1)
end

%%
% This block plots the the selected output as a function of the studied
% values of parameters par1, par2, par3 and par4
clims = [min(output(:, 1)), max(output(:, 1))];
% clims = [0, 1];
output12 = zeros(length(parT1), length(parT2));
for k = 1:length(parT3)
    for l = 1:length(parT4)
        for i = 1:length(parT1)
            for j = 1:length(parT2)
                output12(i, j) =...
                    output(par(:, selPar1) == parT1(i) &...
                    par(:, selPar2) == parT2(j) &... 
                    par(:, selPar3) == parT3(k) &... 
                    par(:, selPar4) == parT4(l), 1);
            end
        end
        
        nfig = nfig + 1;
        figure(nfig)
        imagesc(output12, clims)
        colorbar
        title([char(outputNames(selOut)), ' (', char(parName3), ' = '... 
            num2str(parT3(k)), parUnitTitle3, ', '...
            char(parName4), ' = ', num2str(parT4(l)), parUnitTitle4...
            ')'] , 'fontsize', 22)
        ax = gca;
        ax.FontSize = 22;
        
        xlabel([parName2, parUnit2], 'fontsize', 22)
        ylabel([parName1, parUnit1], 'fontsize', 22)
        xticks(1:length(parT2))
        xticklabels(parT2)
        yticks(1:length(parT1))
        yticklabels(parT1)
    end
end
