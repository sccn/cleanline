function [signal, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn)
% Remove sharp spectral peaks from signal using Sleppian filters
%
% Usage:
% signal = cleanLineNoise(signal)
% [signal, lineNoiseOut] = hcleanLineNoise(signal, lineNoiseIn)
%
% Parameters:
%    signal          Structure with .data and .srate fields
%    lineNoiseIn     Input structure with fields described below
%
% Structure parameters (lineNoiseIn):
%    fPassBand       Frequency band used (default [0, Fs/2] = entire band)
%    Fs 	            Sampling frequency
%    fScanBandWidth  +/- bandwidth centered on each f0 to scan for significant
%                       lines (TM)
%    lineFrequencies Line frequencies to be removed (default
%                       [60, 120, 180, 240, 300])
%    lineNoiseChannels  Channels to remove line noise from (default
%                       size(data, 1))
%    maximumIterations   Maximum times to iterate removal (default = 10)
%    p               Significance level cutoff (default = 0.01)
%    pad             FFT padding factor ( -1 corresponds to no padding,
%                       0 corresponds to padding to next highest power of 2
%                       etc.) (default is 0)
%    pnts
%    tapers          Precomputed tapers from dpss
%    taperBandWidth  Taper bandwidth (default 2 Hz)
%    taperWindowSize Taper sliding window length (default 4 sec)
%    taperWindowStep Sliding window step size (default 4 sec = no overlap)
%    tau             Window overlap smoothing factor (default 100)
%
% This function is based on code originally written by Tim Mullen in a
% package called tmullen-cleanline which is based on the chronux_2
% libraries.
%

lineNoiseOut = lineNoiseIn;
%% Remove line frequencies that are greater than Nyquist frequencies
tooLarge = lineNoiseOut.lineFrequencies >= lineNoiseOut.Fs/2;
if any(tooLarge)
    warning('cleanLineNoise:LineFrequenciesTooLarge', ...
        'Eliminating frequencies greater than half the sampling rate');
    lineNoiseOut.lineFrequencies(tooLarge) = [];
    lineNoiseOut.lineFrequencies = squeeze(lineNoiseOut.lineFrequencies);
end

%% Set up multi-taper parameters
hbw = lineNoiseOut.taperBandWidth/2;   % half-bandwidth
lineNoiseOut.taperTemplate = [hbw, lineNoiseOut.taperWindowSize, 1];
Nwin = round(lineNoiseOut.Fs*lineNoiseOut.taperWindowSize); % number of samples in window
lineNoiseOut.tapers = cleanline_optimized.checkTapers(lineNoiseOut.taperTemplate, Nwin, lineNoiseOut.Fs);
lineNoiseOut.significanceThreshold = finv(1 - lineNoiseOut.p, 2, 2*size(lineNoiseOut.tapers, 2) - 2);
Nfft = max(2^(nextpow2(Nwin) + lineNoiseOut.pad), Nwin);
[lineNoiseOut.frequencyGrid, lineNoiseOut.frequencyIndices] = ...
    getfgrid(lineNoiseOut.Fs, Nfft, lineNoiseOut.fPassBand);
Kodd = 1:2:size(lineNoiseOut.tapers, 2);
H0 = squeeze(sum(lineNoiseOut.tapers(:, Kodd, :), 1));
H0 = H0';
Nf = length(lineNoiseOut.frequencyIndices);
H0 = H0(:, :, ones(1, Nf));
H0 = permute(H0, [3, 1, 2]);
lineNoiseOut.H0 = H0;
lineNoiseOut.H0sq = sum(H0.*H0, 2);

%% Perform the calculation for each channel separately
data = double(signal.data);
chans = sort(lineNoiseOut.lineNoiseChannels);
showProgress = isfield(lineNoiseOut,'showProgress') && lineNoiseOut.showProgress;
progressQueue = [];
if showProgress
    progressWindow = waitbar(0,'Cleaning channels with CleanLine V2...', ...
        'Name','CleanLine Optimized');
    progressCleanup = onCleanup(@() delete(progressWindow)); %#ok<NASGU>
end
mode = 'serial';
if isfield(lineNoiseOut, 'executionMode')
    mode = lineNoiseOut.executionMode;
end
if strcmpi(mode, 'parallel')
    if showProgress
        progressQueue = parallel.pool.DataQueue;
        afterEach(progressQueue, ...
            @(ch) cleanline_optimized_progress(ch,progressWindow,numel(chans)));
    end
    if isempty(chans) || all(diff(chans)==1)
        parfor ch = chans
            data(ch, :) = cleanline_optimized.removeLinesMovingWindow(squeeze(data(ch, :)), lineNoiseOut);
            if showProgress, send(progressQueue,ch); end
        end
    else
        cleaned = cell(1,numel(chans));
        parfor k = 1:numel(chans)
            ch = chans(k);
            cleaned{k} = cleanline_optimized.removeLinesMovingWindow(squeeze(data(ch, :)), lineNoiseOut);
            if showProgress, send(progressQueue,ch); end
        end
        for k = 1:numel(chans)
            data(chans(k), :) = cleaned{k};
        end
    end
else
    for k = 1:numel(chans)
        ch = chans(k);
        data(ch, :) = cleanline_optimized.removeLinesMovingWindow(squeeze(data(ch, :)), lineNoiseOut);
        if showProgress, waitbar(k/numel(chans),progressWindow); end
    end
end
if showProgress, drawnow; end
signal.data = data;
clear data;
