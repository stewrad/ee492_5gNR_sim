function m = evmMerMetricsDD(rxSyms, refConst)
% Decision-directed EVM/MER metrics matching comm.ConstellationDiagram readout.
%   rxSyms   : Nx1 complex equalized symbols for ONE layer
%   refConst : Kx1 complex ideal constellation points

    % Nearest-point decision
    [~, idx] = min(abs(rxSyms - refConst.'), [], 2);
    ref  = refConst(idx);
    err  = rxSyms - ref;
    refP = abs(ref).^2;
    errP = abs(err).^2;

    % RMS and peak EVM as a ratio, then percent
    evmRms = sqrt(mean(errP) / mean(refP));
    evmPk  = sqrt(max(errP)  / mean(refP));

    m.RmsEVM_pct  = 100 * evmRms;
    m.PeakEVM_pct = 100 * evmPk;

    % EVM in dB (20*log10 of the ratio, matching MATLAB scope display)
    m.AvgEVM_dB  = 20 * log10(evmRms);
    m.PeakEVM_dB = 20 * log10(evmPk);

    % MER (Avg) in dB = 10*log10(Pref / Perr)
    m.AvgMER_dB = 10 * log10(mean(refP) / mean(errP));
end