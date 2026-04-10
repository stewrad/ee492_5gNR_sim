function m = evmCommEVM(rxSyms, refConst)
    if length(rxSyms) ~= length(refConst)
        error('rxSyms and txRefSyms must have the same length.');
    end

    evm = comm.EVM(MaximumEVMOutputPort=true);
    [rmsEVM, maxEVM] = evm(refConst, rxSyms);

    err  = rxSyms - refConst;
    refP = abs(refConst).^2;
    errP = abs(err).^2;

    m.RmsEVM_pct  = rmsEVM;
    m.PeakEVM_pct = maxEVM;
    
    evmRatio = rmsEVM / 100;
    evmRatioPk = maxEVM / 100;
    m.AvgEVM_dB  = 20 * log10(max(evmRatio));
    m.PeakEVM_dB = 20 * log10(max(evmRatioPk));
    m.AvgMER_dB = 10 * log10(mean(refP) / max(mean(errP)));
end