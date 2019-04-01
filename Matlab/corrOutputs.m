clear all
close all

%%
% This block initialises nfig and defines
% - path: the path of the ouput directory,
% - nTissues: the number of tissues,
% - fileNames: the names of the ouput files,
% - ouputNames: the names of the outputs,
% - outputCol: the column number of the output in the corresponding files.

nfig = 0;
path = '../../Carlos/Results/Corr/Dose_5Val_5Rep/Tissue';
% Simulations were performed using the input files *Res.dat. To study the
% results of simulation considering all biological processes, output files
% X_1_X.res are used.
nTissues = 21;
P = 1;

fileNames = {'/tumDens','/vascDens', '/vascDens', '/vascDens'...
    '/deadDens', '/pO2Stat', '/pO2Stat'};

outputNames = {'tumour density', 'vascular density'...
    'pre-existing vascular density', 'neo-created vascular density'...
    'dead cells density', 'median pO2', 'mean pO2'};

outputCol = [2, 2, 3, 4, 2, 2, 3];

%%
% This block loads the initial vascular densities values for the 21
% tissues.

initVascDens = load('../HistSpec/initVascDens.dat');

%%
% This block studies the correlation between two time-dependent outputs
% selOut1 and selOut2. For every combination of parameters (values of
% dose), the following indicators are calculated:
% - sim: <u, v> / (||u|| * ||v||),
% - lsd: Sum of ui^2 - vi^2,
% - R: the correlation coefficients,
% - p: the polynomial of degree 1 fitting the curve u = v.

% For the case of nTissues = 0, these indicators are calculated
% independently for each tissue.

% The following intermidiate variables are considered:
% - par: a matrix containing the combinations of parameters simulated. Each
% row corresponds to a simulation and each column, to a parameter,
% - meanOutput1i: a matrix containing the mean values for P repetitions of
% output 1 for the combination of parameters i,
% - meanOutput2i: a matrix containing the mean values for P repetitions of
% output 2 for the combination of parameters i,
% - meanOutput1: a cell array containing the mean values for P repetitions
% of output 1 for all the combination of parameters,
% - meanOutput2: a cell array containing the mean values for P repetitions
% of output 2 for all the combination of parameters,
% - stdOutput1: a cell array containing the std values for P repetitions of
% output 1 for all the combination of parameters,
% - stdOutput2: a cell array containing the std values for P repetitions of
% output 2 for all the combination of parameters.

nTissue = input(['Select one tissue (from 1 to ', num2str(nTissues)...
    ') or all of them (0) or a mean over them (-1): ']);

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);


if(nTissue >= 1 && nTissue <= nTissues)
    pathTissue = [path, num2str(nTissue)];
    par = load([pathTissue, '/combPar.res']);
    nCombPar = size(par, 1);
    colDose = 1;
    tDose = unique(par(:, colDose));
    
    meanOutput1 = cell(1, nCombPar);
    meanOutput2 = cell(1, nCombPar);
    stdOutput1 = cell(1, nCombPar);
    stdOutput2 = cell(1, nCombPar);
    
    sim = zeros(nCombPar, 1);
    lsd = zeros(nCombPar, 1);
    R = zeros(2, 2, nCombPar);
    p = zeros(nCombPar, 2);
    
    for i = 1:nCombPar
        clear output1 output2
        for j = 1:P
            temp1 = load([pathTissue, char(fileNames(selOut1)), '_'...
                num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            temp2 = load([pathTissue, char(fileNames(selOut2)), '_'...
                num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            output1(:, :, j) = temp1(:, [1, outputCol(selOut1)]);
            output2(:, :, j) = temp2(:, [1, outputCol(selOut2)]);
        end
        
        meanOutput1i = mean(output1, 3);
        meanOutput2i = mean(output2, 3);
        meanOutput1(i) = {meanOutput1i};
        meanOutput2(i) = {meanOutput2i};
        
        stdOutput1(i) = {std(output1, 0, 3)};
        stdOutput2(i) = {std(output2, 0, 3)};
        
        sim(i) = meanOutput1i(2:end, 2)' * meanOutput2i(2:end, 2) /...
            (norm(meanOutput1i(2:end, 2)) * norm(meanOutput2i(2:end, 2)));
        
        lsd(i) = sum((meanOutput1i(2:end, 2) - meanOutput2i(2:end, 2)).^2);
        
        R(:, :, i) = corrcoef(meanOutput1i(2:end, 2),...
            meanOutput2i(2:end, 2));
        
        [p(i, :), S] = polyfit(meanOutput1i(2:end, 2)',...
            meanOutput2i(2:end, 2)', 1);
    end
end

if(nTissue == 0)
    par = load([path, '1/combPar.res']);
    nCombPar = size(par, 1);
    meanOutput1 = cell(1, nCombPar);
    meanOutput2 = cell(1, nCombPar);
    stdOutput1 = cell(1, nCombPar);
    stdOutput2 = cell(1, nCombPar);
    colDose = 1;
    tDose = unique(par(:, colDose));
    sim = zeros(nCombPar, nTissues);
    lsd = zeros(nCombPar, nTissues);
    R = zeros(2, 2, nCombPar, nTissues);
    p = zeros(nCombPar, 2, nTissues);
    
    for k = 1:nTissues
        pathTissue = [path, num2str(k)];
        for i = 1:nCombPar
            clear output1 output2
            for j = 1:P
                temp1 = load([pathTissue, char(fileNames(selOut1)), '_'...
                    num2str(i - 1), '_1_', num2str(j - 1), '.res']);
                temp2 = load([pathTissue, char(fileNames(selOut2)), '_'...
                    num2str(i - 1), '_1_', num2str(j - 1), '.res']);
                output1(:, :, j) = temp1(:, [1, outputCol(selOut1)]);
                output2(:, :, j) = temp2(:, [1, outputCol(selOut2)]);
            end
            
            meanOutput1i = mean(output1, 3);
            meanOutput2i = mean(output2, 3);
            meanOutput1(i) = {meanOutput1i};
            meanOutput2(i) = {meanOutput2i};
            
            stdOutput1(i) = {std(output1, 0, 3)};
            stdOutput2(i) = {std(output2, 0, 3)};
            
            sim(i, k) = meanOutput1i(2:end, 2)' *...
                meanOutput2i(2:end, 2) / (norm(meanOutput1i(2:end, 2)) *...
                norm(meanOutput2i(2:end, 2)));
            
            lsd(i, k) = sum((meanOutput1i(2:end, 2) -...
                meanOutput2i(2:end, 2)).^2);
            
            R(:, :, i, k) = corrcoef(meanOutput1i(2:end, 2),...
                meanOutput2i(2:end, 2));
            
            [p(i, :, k), S] = polyfit(meanOutput1i(2:end, 2)',...
                meanOutput2i(2:end, 2)', 1);
        end
    end
end
%%
pAllDoses = mean(p, 1);
pAllDosesStd = std(p, 0, 1);
pAllDoses = permute(pAllDoses, [3, 2, 1]); 
pAllTissues = pAllDoses ./ [initVascDens, initVascDens];

%%
% This block plots, for given combination of parameters (dose value), the
% mean values of time-dependent output 1 as a function of the mean values
% of time-dependent output 2.

dose = 1;
meanOutput1Dose = cell2mat(meanOutput1(dose));
meanOutput2Dose = cell2mat(meanOutput2(dose));

nfig = nfig + 1;
figure(nfig)
plot(meanOutput1Dose(:, 2), meanOutput2Dose(:, 2), '-o')
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

%%
% This block fits with a polynomial of degree 1 the initial values of
% time-dependent output 1 and output 2 for all the tissues. The result
% is then plotted.

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

output1 = zeros(1, nTissues);
output2 = zeros(1, nTissues);
for i = 1:nTissues
    pathTissue = [path, num2str(i)];
    temp = load([pathTissue, char(fileNames(selOut1)),'_0_1_0.res']);
    output1(i) = temp(2, outputCol(selOut1));
    temp = load([pathTissue, char(fileNames(selOut2)),'_0_1_0.res']);
    output2(i) = temp(2, outputCol(selOut2));
end

[sortOutput1, ind] = sort(output1);
sortOutput2 = output2(ind);
p = polyfit(sortOutput1, sortOutput2, 1);

nfig = nfig + 1;
figure(nfig)
hold on
plot(sortOutput1, sortOutput2, '-o')
plot(sortOutput1, p(1) * sortOutput1 + p(2))
hold off
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

%%
% This block fits with a polynomial of degree 1 the time-dependent output 1
% and output 2. The result is then plotted.

selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5),\n'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

nTissue = 16;

pathTissue = [path, num2str(nTissue)];
par = load([pathTissue, '/combPar.res']);
nCombPar = size(par, 1);

meanOutput1 = cell(1, nCombPar);
meanOutput2 = cell(1, nCombPar);
stdOutput1 = cell(1, nCombPar);
stdOutput2 = cell(1, nCombPar);

p = zeros(nCombPar, 2);

for i = 1:nCombPar
    for j = 1:P
        temp1(:, :, j) = load([pathTissue, char(fileNames(selOut1))...
            '_', num2str(i - 1), '_1_', num2str(j - 1), '.res']);
        temp2(:, :, j) = load([pathTissue, char(fileNames(selOut2))...
            '_', num2str(i - 1), '_1_', num2str(j - 1), '.res']);
    end
    temp1 = mean(temp1, 3);
    temp2 = mean(temp2, 3);
    meanOutput1(i) = {temp1};
    meanOutput2(i) = {temp2};
    
    p(i, :) = polyfit(temp1(2:end, outputCol(selOut1)),...
        temp2(2:end, outputCol(selOut2)), 1);
end

meanP = mean(p, 1);
stdP = std(p, 0, 1);

for i = 1:nCombPar
    nfig = nfig + 1;
    figure(nfig)
    hold on
    temp1 = cell2mat(meanOutput1(i));
    temp2 = cell2mat(meanOutput2(i));
    plot(temp1(2:end, outputCol(selOut1)), temp2(2:end, outputCol(selOut2)))
    plot(temp1(2:end, outputCol(selOut1)), meanP(1) *...
        temp1(2:end, outputCol(selOut1)) + meanP(2))
    hold off
    grid on
    xlabel(char(outputNames(selOut1)))
    ylabel(char(outputNames(selOut2)))
    title(['Tissue ', num2str(nTissue), ' - ', num2str(par(i)), ' Gy'])
    
    nfig = nfig + 1;
    figure(nfig)
    hold on
    plot(temp2(2:end, 1), temp2(2:end, outputCol(selOut2)))
    plot(temp2(2:end, 1), meanP(1) * temp1(2:end, outputCol(selOut1)) +...
        meanP(2))
    hold off
    grid on
    xlabel('t(h)')
    ylabel(char(outputNames(selOut2)))
    title(['Tissue ', num2str(nTissue), ' - ', num2str(par(i)), ' Gy'])
end

