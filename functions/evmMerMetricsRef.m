function m = evmMerMetricsRef(rxSyms, txRefSyms)
% EVM/MER using actual transmitted reference symbols
% rxSyms    : Nx1 complex equalized received symbols for ONE layer
% txRefSyms : Nx1 complex transmitted ideal symbols for the same layer

    rxSyms    = rxSyms(:);
    txRefSyms = txRefSyms(:);

    if length(rxSyms) ~= length(txRefSyms)
        error('rxSyms and txRefSyms must have the same length.');
    end

    err  = rxSyms - txRefSyms;
    refP = abs(txRefSyms).^2;
    errP = abs(err).^2;

    evmRms = sqrt(mean(errP) / mean(refP));
    evmPk  = sqrt(max(errP)  / mean(refP));

    m.RmsEVM_pct  = 100 * evmRms;
    m.PeakEVM_pct = 100 * evmPk;

    m.AvgEVM_dB  = 20 * log10(max(evmRms, 1e-12));
    m.PeakEVM_dB = 20 * log10(max(evmPk, 1e-12));

    m.AvgMER_dB = 10 * log10(mean(refP) / max(mean(errP), 1e-12));
end