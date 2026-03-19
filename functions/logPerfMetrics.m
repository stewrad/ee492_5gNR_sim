function [state, row] = logPerfMetrics(state, nSlot, trBlkSizes, blkerr, harqEntity, slotDuration_s)
%LOGPERFMETRICS Update cumulative performance metrics and create a per-slot log row.
%
% This version logs BOTH:
%   1) Attempt BLER (pre-HARQ): based on blkerr for this transmission attempt
%   2) Final BLER (post-HARQ): counts a TB failure only when HARQ times out
%
% Inputs
%   state: struct holding cumulative counters (see expected fields below)
%   nSlot: zero-based slot index
%   trBlkSizes: transport block sizes for each codeword (vector)
%   blkerr: block error flags for each codeword (vector, 0/1 or logical)
%   harqEntity: HARQEntity object (used to detect NewData/retx/timeout)
%   slotDuration_s: duration of a slot in seconds
%
% Outputs
%   state: updated cumulative counters
%   row: struct with per-slot + cumulative metrics snapshot
%
% Expected state fields (initialize in main script):
%   totalTxBits, totalRxBits, totalTransmissions, totalInitialTransmissions,
%   totalRetransmissions, blockErrors, successfulBlocks,
%   attemptsTotal, attemptsFailed

% ---------------- Per-slot tallies ----------------
txBitsThisSlot = 0;
rxBitsThisSlot = 0;
blkErrCountThisSlot = 0;
isNewTB = false;

% One transmission attempt occurred this slot (per codeword)
state.totalTransmissions = state.totalTransmissions + 1;

% -------- Attempt BLER tracking (pre-HARQ) --------
% Count every codeword attempt and its immediate success/failure
state.attemptsTotal  = state.attemptsTotal  + numel(blkerr);
state.attemptsFailed = state.attemptsFailed + sum(blkerr);

attemptBLER_thisSlot = mean(blkerr);                    % pre-HARQ BLER for this slot
attemptBLER_cum      = state.attemptsFailed / max(1,state.attemptsTotal);

% ---------------- TB-level accounting ----------------
for cwIdx = 1:numel(trBlkSizes)

    if harqEntity.NewData(cwIdx)
        isNewTB = true;
        state.totalInitialTransmissions = state.totalInitialTransmissions + 1;

        state.totalTxBits = state.totalTxBits + trBlkSizes(cwIdx);
        txBitsThisSlot = txBitsThisSlot + trBlkSizes(cwIdx);
    else
        state.totalRetransmissions = state.totalRetransmissions + 1;
    end

    if ~blkerr(cwIdx)
        state.successfulBlocks = state.successfulBlocks + 1;

        state.totalRxBits = state.totalRxBits + trBlkSizes(cwIdx);
        rxBitsThisSlot = rxBitsThisSlot + trBlkSizes(cwIdx);
    else
        blkErrCountThisSlot = blkErrCountThisSlot + 1;

        % Final BLER increments ONLY when the RV sequence times out (post-HARQ TB failure)
        if harqEntity.SequenceTimeout(cwIdx)
            state.blockErrors = state.blockErrors + 1;
        end
    end
end

% ---------------- Cumulative metrics ----------------
% Final BLER (post-HARQ)
if state.totalInitialTransmissions > 0
    finalBLER = state.blockErrors / state.totalInitialTransmissions;
else
    finalBLER = NaN;
end

% Throughput efficiency (goodput / transmitted bits)
if state.totalTxBits > 0
    throughputEfficiencyPct = 100 * (state.totalRxBits / state.totalTxBits);
else
    throughputEfficiencyPct = NaN;
end

% Per-slot and average throughput (goodput-based)
instThroughputMbps = (rxBitsThisSlot / slotDuration_s) / 1e6;

elapsedTime_s = (nSlot + 1) * slotDuration_s;
avgThroughputMbps = (state.totalRxBits / max(eps,elapsedTime_s)) / 1e6;

% ---------------- Build log row ----------------
row = struct();
row.slot = nSlot;
row.isNewTB = isNewTB;

row.txBitsThisSlot = txBitsThisSlot;
row.rxBitsThisSlot = rxBitsThisSlot;
row.blkErrCountThisSlot = blkErrCountThisSlot;

% Attempt BLER (pre-HARQ)
row.attemptBLER_thisSlot = attemptBLER_thisSlot;
row.attemptBLER_cum      = attemptBLER_cum;

% Final BLER (post-HARQ)
row.finalBLER = finalBLER;

% Cumulative counters snapshot
row.totalTransmissions = state.totalTransmissions;
row.totalInitialTransmissions = state.totalInitialTransmissions;
row.totalRetransmissions = state.totalRetransmissions;
row.successfulBlocks = state.successfulBlocks;
row.blockErrors = state.blockErrors;

row.throughputEfficiencyPct = throughputEfficiencyPct;
row.instThroughputMbps = instThroughputMbps;
row.avgThroughputMbps = avgThroughputMbps;

end