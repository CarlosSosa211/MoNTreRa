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

sumSI = sum(SI, 2);
sumTSI = sum(TSI, 2);

log10N = log10(N);
nfig = nfig + 1;
figure(nfig)
plot(log10N, log10(errSI), '-o');
title('Error SI')
xlabel('log_{10}(N)')
ylabel('log_{10}(err)')
grid on
legend

nfig = nfig + 1;
figure(nfig)
plot(log10N, log10(errTSI), '-o');
title('Error TSI')
xlabel('log_{10}(N)')
ylabel('log_{10}(err)')
grid on
legend

nfig = nfig + 1;
figure(nfig)
plot(log10N, sumSI, '-o');
title('Sum SI')
xlabel('log_{10}(N)')
ylabel('\Sigma SI_i')
grid on

nfig = nfig + 1;
figure(nfig)
plot(log10N, sumTSI, '-o');
title('Sum TSI')
xlabel('log_{10}(N)')
ylabel('\Sigma TSI_i')
grid on