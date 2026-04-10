clc; clear;

snrRange = 0:2:10;
harqList = [1 2 4 8 16];

for h = 1:length(harqList)

    NHARQ = harqList(h);

    for s = 1:length(snrRange)

        SNRdB = snrRange(s);

        fprintf('\n===== HARQ=%d | SNR=%.1f dB =====\n', NHARQ, SNRdB);

        nr_dlsch_tx_rx_BASELINE(SNRdB, NHARQ);

    end
end


folder = '\\wsl.localhost\Ubuntu\home\arnold_kalala\ee492_5gNR_sim\logs';

files = dir(fullfile(folder,'*.txt'));

harqList = [1 2 4 8 16];

data = struct();

for f = 1:length(files)

    fname = files(f).name;

    % --- Extract SNR ---
    snrMatch = regexp(fname,'SNR([0-9.]+)','tokens');
    if isempty(snrMatch), continue; end
    SNR = str2double(snrMatch{1});

    % --- Extract HARQ ---
    harqMatch = regexp(fname,'NHARQ(\d+)','tokens');
    if isempty(harqMatch), continue; end
    HARQ = str2double(harqMatch{1});

    % --- Read file ---
    fullpath = fullfile(folder,fname);
    txt = fileread(fullpath);

    % --- Extract retransmission rate ---
    rateMatch = regexp(txt,'Retransmission Rate:\s*([\d.]+)','tokens');

    if ~isempty(rateMatch)
        rate = str2double(rateMatch{end}); % take FINAL value
    else
        rate = NaN;
    end

    % --- Store ---
    if ~isfield(data, sprintf('H%d',HARQ))
        data.(sprintf('H%d',HARQ)).SNR = [];
        data.(sprintf('H%d',HARQ)).rate = [];
    end

    data.(sprintf('H%d',HARQ)).SNR(end+1) = SNR;
    data.(sprintf('H%d',HARQ)).rate(end+1) = rate;

end

% =======================
% PLOT
% =======================
figure; hold on;

for h = harqList
    fieldName = sprintf('H%d',h);

    if isfield(data, fieldName)

        SNR = data.(fieldName).SNR;
        rate = data.(fieldName).rate;

        % sort
        [SNR, idx] = sort(SNR);
        rate = rate(idx);

        plot(SNR, rate, '-o', 'LineWidth', 2, ...
            'DisplayName', sprintf('HARQ=%d',h));
    end
end

grid on;
xlabel('SNR (dB)');
ylabel('Retransmission Rate (%)');
title('Retransmission Rate vs SNR (256QAM)');
legend('show','Location','northeast');