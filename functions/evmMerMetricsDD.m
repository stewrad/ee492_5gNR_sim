function m = evmMerMetricsDD(rxSyms, refConst)
% Decision-directed EVM/MER metrics similar to Constellation Diagram readout.
% rxSyms: Nx1 complex symbols for ONE layer (equalized)
% refConst: Kx1 complex ideal constellation points

    % Nearest-point decisions (decision-directed reference)
    [~, idx] = min(abs(rxSyms - refConst.'), [], 2);
    ref = refConst(idx);

    err = rxSyms - ref;

    refP = abs(ref).^2;
    errP = abs(err).^2;

    % RMS and peak EVM (ratio), then percent
    evmRms = sqrt(mean(errP) / mean(refP));     % ratio (e.g., 0.18)
    evmPk  = sqrt(max(errP) / mean(refP));      % ratio

    m.RmsEVM_pct  = 100 * evmRms;
    m.PeakEVM_pct = 100 * evmPk;

    % EVM in dB as shown by MATLAB scopes (20*log10 of percent value)
    m.AvgEVM_dB  = 20*log10(evmRms);
    m.PeakEVM_dB = 20*log10(evmPk);

    % MER (Avg) in dB = 10*log10(Pref / Perr)
    m.AvgMER_dB = 10*log10(mean(refP) / mean(errP));
end