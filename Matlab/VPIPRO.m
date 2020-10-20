clear all
close all
nfig = 0;

tId = {'02', '03', '04', '05', '06', '08', '09', '10', '11', '12', '14',...
    '15', '16', '17', '18' '19', '20', '21'};
nPatients = length(tId);

tIdPSA = {'01', '02', '03', '05', '09', '15', '16', '19'};
nPatientsPSA = length(tIdPSA);

data = zeros(5, 4, nPatients);
PSA = cell(nPatientsPSA, 1);

for i = 1:nPatients
    data(:, 1:2, i) = load(['tumVolADC/tumVolADC0', char(tId{i}), '.dat']);
end

for i = 1:nPatients
    firstNz = find(data(:, 1, i));
    if(firstNz)
        data(:, 3, i) = data(:, 1, i) / data(find(data(:, 1, i), 1), 1, i);
    end
    
    firstNz = find(data(:, 2, i));
    if(firstNz)
        data(:, 4, i) = data(:, 2, i) / data(find(data(:, 2, i), 1), 2, i);
    end
end

for i = 1:nPatientsPSA
    temp = load(['PSA/PSA0', char(tIdPSA{i}), '.dat']);
    temp(:, 3) = temp(:, 2) / temp (1, 2);
    PSA{i} = temp;
end

%%
days = 0:14:56;

for i = 1:nPatients
    nfig = nfig + 1;
    figure(nfig)
    nonZeros = find(data(:, 1, i));
    plot(days(nonZeros), data(nonZeros, 1, i), '-o', 'linewidth', 2)
    xlabel('t (days)')
    ylabel('Tumor volume (mm^3)')
    title(['Tumor volume - ID0', char(tId{i})])
    
    nfig = nfig + 1;
    figure(nfig)
    nonZeros = find(data(:, 2, i));
    plot(days(nonZeros), data(nonZeros, 2, i), '-o', 'linewidth', 2)
    xlabel('t (days)')
    ylabel('ADC')
    title(['ADC - ID0', char(tId{i})])
end

%%
days = 0:14:56;

nfig = nfig + 1;
figure(nfig)
hold on
leg = {};
for i = 1:nPatients
    nonZeros = find(data(:, 1, i));
    if (length(nonZeros) > 3)
        plot(days(nonZeros), data(nonZeros, 1, i), '-o', 'linewidth', 2)
        leg(length(leg) + 1) = tId(i);
    end
end
hold off
xlabel('t (days)')
ylabel('Tumor volume (mm^3)')
title('Tumor volume')
legend(leg)
grid on

nfig = nfig + 1;
figure(nfig)
hold on
leg = {};
for i = 1:nPatients
    nonZeros = find(data(:, 2, i));
    if (length(nonZeros) > 3)
        plot(days(nonZeros), data(nonZeros, 2, i), '-o', 'linewidth', 2)
        leg(length(leg) + 1) = tId(i);
    end
end
hold off
xlabel('t (days)')
ylabel('ADC')
title('ADC')
legend(leg)
grid on

nfig = nfig + 1;
figure(nfig)
hold on
leg = {};
for i = 1:nPatients
    nonZeros = find(data(:, 1, i));
    if (length(nonZeros) > 3)
        plot(days(nonZeros), data(nonZeros, 3, i), '-o', 'linewidth', 2)
        leg(length(leg) + 1) = tId(i);
    end
end
hold off
xlabel('t (days)')
ylabel('Relative tumor volume')
title('Relative tumor volume')
legend(leg)
grid on

nfig = nfig + 1;
figure(nfig)
hold on
leg = {};
for i = 1:nPatients
    nonZeros = find(data(:, 2, i));
    if(length(nonZeros) > 3)
        plot(days(nonZeros), data(nonZeros, 4, i), '-o', 'linewidth', 2)
        leg(length(leg) + 1) = tId(i);
    end
end
hold off
xlabel('t (days)')
ylabel('Relative ADC')
title('Relative ADC')
legend(leg)
grid on

%%
nfig = nfig + 1;
figure(nfig)
hold on
for i = 1:nPatients
    nonZeros = find(data(:, 1, i) & data(:, 2, i)) ;
    plot(data(nonZeros, 1, i), data(nonZeros, 2, i), '-o', 'linewidth', 2)
end
hold off
xlabel('relative tumor volume')
ylabel('relaive ADC')
title('ADC')
grid on

%%
for i = 1:nPatients
    nfig = nfig + 1;
    figure(nfig)
    nonZeros = find(data(:, 1, i) & data(:, 2, i)) ;
    plot(data(nonZeros, 1, i), data(nonZeros, 2, i), '-o', 'linewidth', 2)
    xlabel('relative tumor volume')
    ylabel('relaive ADC')
    title('ADC')
    grid on
end

%%
for i = 1:nPatientsPSA
    nfig = nfig + 1;
    figure(nfig)
    temp = PSA{i};
    plot(temp(:, 1), temp(:, 2), '-o')
    title(['ID0', char(tIdPSA{i})])
end

%%
nfig = nfig + 1;
figure(nfig)
hold on
for i = 1:nPatientsPSA
    temp = PSA{i};
    plot(temp(:, 1), temp(:, 4), '-o', 'linewidth', 2)
end
hold off
xlabel('t (days)')
ylabel('Relative PSA')
title('PSA')
grid on



