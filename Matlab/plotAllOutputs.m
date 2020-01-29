function plotAllOutputs(path, nTissue)
global nPar selPar
global b nfig;

endTreatTumDens = load([path, '/morrisEndTreatTumDens_'...
    num2str(nTissue), '.res']);
endTreatTumDens = endTreatTumDens(selPar == 1, :);
endTreatTumDens(:, 3) = sqrt(endTreatTumDens(:, 1).^2 +...
    endTreatTumDens(:, 2).^2);
maxMeanEndTreatTumDens_ = 1./max(endTreatTumDens(:, 3));
endTreatTumDens(:, 3) = endTreatTumDens(:, 3) .* maxMeanEndTreatTumDens_;

threeMonTumDens = load([path, '/morris3MonTumDens_', num2str(nTissue)...
    '.res']);
threeMonTumDens = threeMonTumDens(selPar == 1, :);
threeMonTumDens(:, 3) = sqrt(threeMonTumDens(:, 1).^2 +...
    threeMonTumDens(:, 2).^2);
maxMean3MonTumDens_ = 1./max(threeMonTumDens(:, 3));
threeMonTumDens(:, 3) = threeMonTumDens(:, 3) .* maxMean3MonTumDens_;

killed99 = load([path, '/morrisKilled99_', num2str(nTissue), '.res']);
killed99 = killed99(selPar == 1, :);
killed99(:, 3) = sqrt(killed99(:, 1).^2 + killed99(:, 2).^2);
maxMeanKilled99_ = 1./max(killed99(:, 3));
killed99(:, 3) = killed99(:, 3) .* maxMeanKilled99_;

intTumDens = load([path, '/morrisIntTumDens_', num2str(nTissue), '.res']);
intTumDens = intTumDens(selPar == 1, :);
intTumDens(:, 3) = sqrt(intTumDens(:, 1).^2 + intTumDens(:, 2).^2);
maxMeanIntTumDens_ = 1./max(intTumDens(:, 3));
intTumDens(:, 3) = intTumDens(:, 3) .* maxMeanIntTumDens_;

cOut = [num2cell(endTreatTumDens(:, 3)), num2cell(threeMonTumDens(:, 3))...
    num2cell(killed99(:, 3)), num2cell(intTumDens(:, 3)),  b'];
cOut = sortrows(cOut, 3);

nfig = nfig + 1;
figure(nfig);
bar(cell2mat(cOut(:, 1:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cOut(:, 5));
ax.YGrid = 'on';
ax.FontSize = 26;
ylim([0, inf])
title(['Tissue ', num2str(nTissue), ' - All outputs'], 'fontsize', 26)
legend({'Tumour density at the end of treat.'...
    'Tumour density 3 months after the end of treat.'...
    '99% of initial tumour cells killed', 'Integral of tumour density'},...
    'location', 'northwest', 'fontsize', 26)
xtickangle(45)