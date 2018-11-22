function plotMeanTimeTo95(path)
global nPar
global allTissues
global densTissues nonDensTissues;
global vascTissues nonVascTissues;
global b color nfig shape;

tissueSet = input(['Define a set [all (1), dense (2), non-dense (3), ' ...
    'vascularized (4) or non-vascularized (5) tissues]: ']);

switch tissueSet
    case 1
        tTissues = allTissues;
    case 2
        tTissues = densTissues;
    case 3
        tTissues = nonDensTissues;
    case 4
        tTissues = vascTissues;
    case 5
        tTissues = nonVascTissues;
end

nTissues = length(tTissues);
timeTo95 = zeros(nPar, 4, nTissues);

for i = 1:length(tTissues)
    timeTo95(:, 1:2, i) = load([path, '/morrisTimeTo95_', num2str(tTissues(i)), '.res']);
    timeTo95(:, 4, i) = sqrt(timeTo95(:, 1, i).^2 + timeTo95(:, 2, i).^2);
    timeTo95(:, 3, i) = timeTo95(:, 1, i).^2 ./ timeTo95(:, 4, i); 
end

meanTimeTo95 = mean(timeTo95, 3);
stdTimeTo95 = std(timeTo95, [], 3);

cTimeTo95 = [num2cell(meanTimeTo95), num2cell(stdTimeTo95), b'];
cTimeTo95 = sortrows(cTimeTo95, 4);

figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(timeTo95, 1)
    scatter(meanTimeTo95(i,1), meanTimeTo95(i,2), 500, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])], [0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])], '--k')
hold off

switch tissueSet
    case 1
        title('21 tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 2
        title('11 dense tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
        
end
axis([0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)]), 0, 1.1 * max([timeTo95(:, 1); timeTo95(:, 2)])])
grid on
legend(b, 'fontsize', 18, 'location', 'bestoutside', 'interpreter', 'latex')
xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)

figure(nfig);
nfig = nfig + 1;

hold on
hBar = bar(cell2mat(cTimeTo95(:, 3:4)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTimeTo95(:, 9), 'fontsize', 20);
ax.YGrid = 'on';
switch tissueSet
    case 1
        title('21 tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 2
        title('11 dense tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 3
        title('10 non-dense tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 4
        title('11 vascularized tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
    case 5
        title('10 non-vascularized tissues - Time to kill 95% of tumor cells', 'fontsize', 20)
        
end
xtickangle(45)

ctr = zeros(nPar, 2);
ydt = zeros(nPar, 2);
for k = 1:2
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cTimeTo95(:, 7:8)), '.k')
hold off
ylim([0, inf])
legend({'$\frac{\mu*^2}{\sqrt{\mu*^2 + \sigma^2}}$', '$\sqrt{\mu*^2 + \sigma^2}$'},...
    'fontsize', 20, 'location', 'northwest', 'interpreter', 'latex')
end
