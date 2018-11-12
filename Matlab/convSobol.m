close all
nfig = 0;

SI = load('../OutputFiles/convSI.res');
TSI = load('../OutputFiles/convTSI.res');
sol = load('../OutputFiles/sobol.res');
K = sol(1, 1);
N = SI(:, 1);
SI = SI(:, 2:K + 1);
TSI = TSI(:, 2: K + 1);
sol = sol(2:K + 1, :);

for k = 1:K
    errSI(k, :) = abs(SI(:, k) - sol(k, 1))';
    errTSI(k, :) = abs(TSI(:, k) - sol(k, 1))';
end

nfig = nfig + 1;
figure(nfig)
plot(log10(N), log10(errSI));
title('SI')

nfig = nfig + 1;
figure(nfig)
plot(log10(N), log10(errTSI));
title('TSI')