function plotAllOutputs(path, nTissue)
global nPar
global b nfig;

endTreatTumDens = load([path, '/morrisEndTreatTumDens_', num2str(nTissue), '.res']);
endTreatTumDens(:, 3) = sqrt(endTreatTumDens(:, 1).^2 + endTreatTumDens(:, 2).^2);
threeMonTumDens = load([path, '/morris3MonTumDens_', num2str(nTissue), '.res']);
threeMonTumDens(:, 3) = sqrt(threeMonTumDens(:, 1).^2 + threeMonTumDens(:, 2).^2);
recTumDens = load([path, '/morrisRecTumDens_', num2str(nTissue), '.res']);
recTumDens(:, 3) = sqrt(recTumDens(:, 1).^2 + recTumDens(:, 2).^2);
tumVol = load([path, '/morrisTumVol_', num2str(nTissue), '.res']);
tumVol(:, 3) = sqrt(tumVol(:, 1).^2 + tumVol(:, 2).^2);
intTumDens = load([path, '/morrisIntTumDens_', num2str(nTissue), '.res']);
intTumDens(:, 3) = sqrt(intTumDens(:, 1).^2 + intTumDens(:, 2).^2);
timeTo95 = load([path, '/morrisTimeTo95_', num2str(nTissue), '.res']);
timeTo95(:, 3) = sqrt(timeTo95(:, 1).^2 + timeTo95(:, 2).^2);
timeTo99 = load([path, '/morrisTimeTo99_', num2str(nTissue), '.res']);
timeTo99(:, 3) = sqrt(timeTo99(:, 1).^2 + timeTo99(:, 2).^2);
recTime = load([path, '/morrisRecTime_', num2str(nTissue), '.res']);
recTime(:, 3) = sqrt(recTime(:, 1).^2 + recTime(:, 2).^2);

maxMeanEndTreatTumDens_ = 1./max(endTreatTumDens(:, 3));
endTreatTumDens(:, 3) = endTreatTumDens(:, 3) .* maxMeanEndTreatTumDens_;
maxMean3MonTumDens_ = 1./max(threeMonTumDens(:, 3));
threeMonTumDens(:, 3) = threeMonTumDens(:, 3) .* maxMean3MonTumDens_;
maxMeanRecTumDens_ = 1./max(recTumDens(:, 3));
recTumDens(:, 3) = recTumDens(:, 3) .* maxMeanRecTumDens_;
maxMeanTumVol_ = 1./max(tumVol(:, 3));
tumVol(:, 3) = tumVol(:, 3) .* maxMeanTumVol_;
maxMeanIntTumDens_ = 1./max(intTumDens(:, 3));
intTumDens(:, 3) = intTumDens(:, 3) .* maxMeanIntTumDens_;
maxMeanTimeTo95_ = 1./max(timeTo95(:, 3));
timeTo95(:, 3) = timeTo95(:, 3) .* maxMeanTimeTo95_;
maxMeanTimeTo99_ = 1./max(timeTo99(:, 3));
timeTo99(:, 3) = timeTo99(:, 3) .* maxMeanTimeTo99_;
maxMeanRecTime_ = 1./max(recTime(:, 3));
recTime(:, 3) = recTime(:, 3) .* maxMeanRecTime_;

cOut = [num2cell(endTreatTumDens(:, 3)), num2cell(threeMonTumDens(:, 3))...
    num2cell(recTumDens(:, 3)), num2cell(tumVol(:, 3))...
    num2cell(intTumDens(:, 3)), num2cell(timeTo95(:, 3))...
    num2cell(timeTo99(:, 3)), num2cell(recTime(:, 3)),  b'];
cOut = sortrows(cOut, 1);

nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cOut(:, 1:8)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOut(:, 9));
ax.YGrid = 'on';
ylim([0, inf])
title(['Tissue ', num2str(nTissue), ' - All outputs'], 'fontsize', 20)
legend({'Tumour density at the end of treat.', 'Tumour density 3 months after the end of treat.'...
    'Tumour density at recurrence beginning', 'Final tumour volume',  'Integral of tumour density'...
    'Time to kill 95%', 'Time to kill 99%', 'Recurrence beginnning time'},...
    'location', 'northwest', 'fontsize', 20)
xtickangle(45)