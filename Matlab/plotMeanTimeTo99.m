function plotMeanTimeTo99(path)
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
timeTo99 = zeros(nPar, 3, nTissues);

for i = 1:length(tTissues)
    timeTo99(:, 1:2, i) = load([path, '/morrisTimeTo99_', num2str(tTissues(i)), '.res']);
    timeTo99(:, 3, i) = sqrt(timeTo99(:, 1, i).^2 + timeTo99(:, 2, i).^2);
end

meanTimeTo99 = mean(timeTo99, 3);
stdTimeTo99 = std(timeTo99, [], 3);

cTimeTo99 = [num2cell(meanTimeTo99), num2cell(stdTimeTo99), b'];
cTimeTo99 = sortrows(cTimeTo99, 3);

figure(nfig);
nfig = nfig + 1;
hold on
colormap(jet)
for i = 1 : size(timeTo99, 1)
    scatter(meanTimeTo99(i,1), meanTimeTo99(i,2), 20, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end
plot([0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])], [0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])], '--k')
legend(b, 'Location', 'bestoutside', 'Interpreter', 'Latex')
xlabel('\mu*')
ylabel('\sigma')
switch tissueSet
    case 1
        title('21 tissues - Time to kill 99% of tumor cells')
    case 2
        title('11 dense tissues - Time to kill 99% of tumor cells')
    case 3
        title('10 non-dense tissues - Time to kill 99% of tumor cells')
    case 4
        title('11 vascularized tissues - Time to kill 99% of tumor cells')
    case 5
        title('10 non-vascularized tissues - Time to kill 99% of tumor cells')
        
end
axis([0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)]), 0, 1.1 * max([timeTo99(:, 1); timeTo99(:, 2)])])
grid on
hold off

figure(nfig);
nfig = nfig + 1;
hold on
hBar = bar(cell2mat(cTimeTo99(:, 1:3)));
ax = gca;
ax.TickLabelInterpreter = 'latex';
set(ax, 'XTick', 1:nPar)
set(ax,'XTickLabel', cTimeTo99(:, 7));
ax.YGrid = 'on';
switch tissueSet
    case 1
        title('21 tissues - Time to kill 99% of tumor cells')
    case 2
        title('11 dense tissues - Time to kill 99% of tumor cells')
    case 3
        title('10 non-dense tissues - Time to kill 99% of tumor cells')
    case 4
        title('11 vascularized tissues - Time to kill 99% of tumor cells')
    case 5
        title('10 non-vascularized tissues - Time to kill 99% of tumor cells')
end
legend('\mu*', '\sigma', 'dist', 'location', 'northwest')
xtickangle(45)

ctr = zeros(nPar, 3);
ydt = zeros(nPar, 3);
for k = 1:3
    ctr(:, k) = bsxfun(@plus, hBar(1).XData, [hBar(k).XOffset]');
    ydt(:, k) = hBar(k).YData;
end

errorbar(ctr, ydt, cell2mat(cTimeTo99(:, 4:6)), '.k')
hold off
end
