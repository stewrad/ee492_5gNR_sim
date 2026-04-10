T = readtable('runlog_results.xlsx');  % add sheet name if needed
head(T)

mods = unique(T.Modulation, 'stable');

figure; hold on; grid on;
colors = lines(numel(mods));

for k = 1:numel(mods)
    idx = strcmp(T.Modulation, mods{k});
    snr = T.SNRdB(idx);
    thr = T.AverageThroughput_Mbps_(idx);

    [snrSorted, order] = sort(snr);
    thrSorted = thr(order);

    plot(snrSorted, thrSorted, '-o', ...
         'DisplayName', mods{k}, ...
         'Color', colors(k,:));
end

xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
title('Throughput vs SNR per Modulation');
legend('Location','best');

mods = unique(T.Modulation, 'stable');

figure; hold on; grid on;
colors = lines(numel(mods));

for k = 1:numel(mods)
    idx = strcmp(T.Modulation, mods{k});
    snr = T.SNRdB(idx);
    bler = T.BlockErrorRate_BLER__(idx);   % 0–1 or 0–100 depending on how you stored it

    [snrSorted, order] = sort(snr);
    blerSorted = bler(order);

    semilogy(snrSorted, max(blerSorted, 1e-4), '-o', ...
             'DisplayName', mods{k}, ...
             'Color', colors(k,:));
end

xlabel('SNR (dB)');
ylabel('BLER');
title('BLER vs SNR per Modulation');
legend('Location','southwest');

T = readtable('runlog_results.xlsx');

% Fix HARQ configuration for a fair comparison
cfgIdx = strcmp(T.HARQType,'Chase Combining') & ...
         T.NHARQProcesses == 16 & ...
         T.NumLayers == 2;          % adjust if needed

Tc = T(cfgIdx, :);
mods = unique(Tc.Modulation, 'stable');

% --- Throughput vs SNR (comparative) ---
figure; hold on; grid on;
colors = lines(numel(mods));

for k = 1:numel(mods)
    idx = strcmp(Tc.Modulation, mods{k});
    snr = Tc.SNRdB(idx);
    thr = Tc.AverageThroughput_Mbps_(idx);

    [snrSorted, order] = sort(snr);
    thrSorted = thr(order);

    plot(snrSorted, thrSorted, '-o', ...
        'DisplayName', mods{k}, 'Color', colors(k,:));
end

xlabel('SNR (dB)');
ylabel('Throughput (Mbps)');
title('Throughput vs SNR (Modulation comparison)');
legend('Location','best');
