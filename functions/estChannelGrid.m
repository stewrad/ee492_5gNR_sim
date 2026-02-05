
function estChannelGrid = getInitialChannelEstimate(channel,carrier)
% Obtain an initial channel estimate for calculating the precoding matrix.
% This function assumes a perfect channel estimate

    % Clone of the channel
    chClone = channel.clone();
    chClone.release();

    % No filtering needed to get channel path gains
    chClone.ChannelFiltering = false;

    % Set channel response output type to calculate perfect channel
    % estimation
    chClone.ChannelResponseOutput = 'ofdm-response';
    
    % Get perfect channel estimate directly from the channel
    estChannelGrid = chClone(carrier);
end
