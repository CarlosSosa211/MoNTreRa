clear all
close all

nfig = 0;
% path = '../../Carlos/Results/Diff_Ang_Dose_5Val_5Rep_AllTissues/Tissue';
path = '../OutputFilesGUI/';
nTissues = 21;

fileNames = {'/tumDens','/vascDens', '/vascDens', '/vascDens'...
    '/deadDens', '/pO2Stat', '/pO2Stat'};

outputNames = {'tumour density', 'vascular density'...
    'pre-existing vascular density', 'neo-created vascular density'...
    'dead cells density', 'median pO2', 'mean pO2'};

outputCol = [2, 2, 3, 4, 2, 2, 3];

%%
meanOutput1 = {};
meanOutput2 = {};
stdOutput1 = {};
stdOutput2 = {};
P = 5;

nTissue = input(['Select one tissue (from 1 to ', num2str(nTissues)...
    ') or all of them (0) or a mean over them (-1): ']);

if(nTissue >= 1 && nTissue <= nTissues)
    pathTissue = [path, num2str(nTissue)];
    par = load([pathTissue, '/combPar.res']);
    
    % colTTum = 1;
    % colDThres = 2;
    % colTArrest = 3;
    colDose = 1;
    tDose = unique(par(:, colDose));
    % tTTum = unique(par(:, colTTum));
    % tDThres = unique(par(:, colDThres));
    % tTArrest = unique(par(:, colTArrest));
    
    selOut1 = input(['Select the firt output [vascDens (1), '...
        'preExVascDens (2), neoCreVascDens (3)\n pO2Med (4), '...
        'pO2Mean (5) or quit (0): ']);
    
    selOut2 = input(['Select the second output [vascDens (1),'...
        'preExVascDens (2), neoCreVascDens (3)\n pO2Med (4),'...
        'pO2Mean (5) or quit (0): ']);
    
    
    for i = 1:size(par, 1)
        clear output1 output2
        for j = 1:P
            temp1 = load([pathTissue, char(fileNames(selOut1))...
                num2str(i - 1), '_1_', num2str(j - 1), '.res']);
            temp2 = load([pathTissue, char(fileNames(selOut2))...
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
        
        sim(i) = meanOutput1i(:, 2)' * meanOutput2i(:, 2) /...
            (norm(meanOutput1i(:, 2)) * norm(meanOutput2i(:, 2)));
        
        lsd(i) = sum((meanOutput1i(:, 2) - meanOutput2i(:, 2)).^2);
        
        R(:, :, i) = corrcoef(meanOutput1i(:, 2), meanOutput2i(:, 2));
        
        [p(i, :), S] = polyfit(meanOutput1i(:, 2)',...
            meanOutput2i(:, 2)', 1);
    end
end

%%
dose = 2;
meanOutput1Dose = cell2mat(meanOutput1(dose));
meanOutput2Dose = cell2mat(meanOutput2(dose));

nfig = nfig + 1;
figure(nfig)
plot(meanOutput1Dose(:, 2), meanOutput2Dose(:, 2))
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

%%
selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5)\n,'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5)\n,'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

output1 = zeros(1, nTissues);
output2 = zeros(1, nTissues);
for i = 1:nTissues
    pathTissue = [path, num2str(i)];
    temp = load([pathTissue, char(fileNames(selOut1)),'_0_0_0.res']);
    output1(i) = temp(2, outputCol(selOut1));
    temp = load([pathTissue, char(fileNames(selOut2)),'_0_0_0.res']);
    output2(i) = temp(2, outputCol(selOut2));
end

[sortOutput1, ind] = sort(output1);
sortpO2 = output2(ind);
p = polyfit(sortOutput1, sortpO2, 1);

nfig = nfig + 1;
figure(nfig)
hold on
plot(sortOutput1, sortpO2, '-o')
plot(sortOutput1, p(1) * sortOutput1 + p(2))
hold off
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

%%
selOut1 = input(['Select the firt output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5)\n,'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

selOut2 = input(['Select the second output [tumDens (1), vascDens (2), '...
    'preExVascDens (3), neoCreVascDens (4), deadDens (5)\n,'...
    'pO2Med (6), pO2Mean (7) or quit (0): ']);

output1 = load([path, char(fileNames(selOut1)),'.res']);
output2 = load([path, char(fileNames(selOut2)),'.res']);

p = polyfit(output1(:, outputCol(selOut1)),...
    output2(:, outputCol(selOut2)), 1);

nfig = nfig + 1;
figure(nfig)
hold on
plot(output1(:, outputCol(selOut1)), output2(:, outputCol(selOut2)), '-o')
hold off
grid on
xlabel(char(outputNames(selOut1)))
ylabel(char(outputNames(selOut2)))

nfig = nfig + 1;
figure(nfig)
hold on
plot(output1(:, 1), output1(:, outputCol(selOut1)))
plot(output2(:, 1), output2(:, outputCol(selOut2)))
hold off
grid on
