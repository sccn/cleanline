function channelData = cleanChannelEpochs(channelData,lineNoise)
% Restart the spectrum, iterative fits and overlap smoothing for each epoch.
for trial = 1:size(channelData,3)
    channelData(1,:,trial) = cleanline_optimized.removeLinesMovingWindow( ...
        channelData(1,:,trial),lineNoise);
end
end
