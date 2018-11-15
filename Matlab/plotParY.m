function plotParY(path, par)
global nPar;
global varRange;
global b color nfig shape;

L = length(varRange);
Y = zeros(nPar, 2, L);

for i = 1:L
    Y(:, :, i) = load([path, '/morrisY_', num2str(i - 1), '.res']);
end

nfig = nfig + 1;
figure(nfig);

hold on
colormap(jet)
for i = 1:L
    scatter(Y(par, 1, i), Y(par, 2, i), 200, color(i), 'filled', shape(mod(i, length(shape)) + 1))
end

maxVal = 1.1 * max([reshape(Y(par, 1, :), 1, []), reshape(Y(par, 2, :), 1, [])]);
plot([0, maxVal], [0, maxVal], '--k')
hold off

xlabel('\mu*', 'fontsize', 20)
ylabel('\sigma', 'fontsize', 20)
titleY = strcat(string(b(par)), ' - Y');
title(titleY, 'interpreter', 'latex', 'fontsize', 20)
axis([0, maxVal, 0, maxVal])
grid on

for i = 1:L
    leg(i) = string(varRange(i));
end
legend(leg,'location', 'bestoutside', 'fontsize', 20)
end
