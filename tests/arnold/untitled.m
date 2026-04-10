%% plot_results.m
% Post-processing plots for runlog_results.xlsx

clear; close all; clc;
T = readtable("runlog_results.xlsx");

metrics = [
    T.BlockErrorRate_BLER__
    T.RetransmissionRate___
    T.AverageTransmissionsPerTB
    T.AverageThroughput_Mbps_
    T.SpectralEfficiency_bits_s_Hz_
];

labels = {
    "BLER"
    "Retransmission Rate (%)"
    "Avg Tx per TB"
    "Avg Throughput (Mbps)"
    "Spectral Eff (bits/s/Hz)"
};

figure
bar(metrics)
set(gca,"XTickLabel",labels,"XTickLabelRotation",30)
title("Single-Scenario Summary Metrics")
grid on

figure
bar([T.Layer1_MeanSINR, T.Layer2_MeanSINR])
set(gca,"XTickLabel",["Layer 1","Layer 2"])
ylabel("Mean SINR (dB)")
title("Per-Layer SINR (Single Scenario)")
grid on
