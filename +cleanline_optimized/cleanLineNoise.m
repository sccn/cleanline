function [signal, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn)
% Remove sharp spectral peaks from signal using Sleppian filters
%
% Usage:
% signal = cleanLineNoise(signal)
% [signal, lineNoiseOut] = hcleanLineNoise(signal, lineNoiseIn)
%
% Parameters:
%    signal          Structure with .data and .srate fields. Data are
%                    channels x samples x epochs, with each epoch cleaned independently.
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
if size(signal.data,3) > 1
    warning('CleanLine:EpochedDataDiscouraged', ...
        ['Applying CleanLine to epoched EEG is strongly discouraged. ' ...
         'Short epochs limit the reliability of multitaper spectral estimation and sinusoidal fitting, ' ...
         'potentially resulting in inconsistent line-noise removal across trials. ' ...
         'Whenever possible, remove line noise from continuous EEG before epoching. ' ...
         'Processing will continue with epoch-independent estimation.']);
    epochSamples = size(signal.data,2);
    Nstep = round(lineNoiseOut.Fs*lineNoiseOut.taperWindowStep);
    if ~isfinite(Nwin) || Nwin < 1 || Nwin > epochSamples
        error('CleanLine:EpochWindowLength', ...
            'Window length must round to 1 to %d samples (one epoch). Got %g samples.', ...
            epochSamples,Nwin);
    end
    if ~isfinite(Nstep) || Nstep < 1 || Nstep > Nwin
        error('CleanLine:EpochWindowStep', ...
            'Window step must round to 1 to %d samples (one window). Got %g samples.', ...
            Nwin,Nstep);
    end
    timeBandwidth = lineNoiseOut.taperWindowSize*hbw;
    if ~isfinite(timeBandwidth) || floor(2*timeBandwidth-1) < 2 || timeBandwidth >= Nwin/2
        error('CleanLine:EpochTapers', ...
            ['A %.6g s epoch window with taper bandwidth %.6g Hz cannot support the F-test. ' ...
             'Choose a window and taper bandwidth that yield at least two tapers ' ...
             'and a time-bandwidth product below half the window sample count.'], ...
            lineNoiseOut.taperWindowSize,lineNoiseOut.taperBandWidth);
    end
end
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
            data(ch,:,:) = cleanChannelEpochs(data(ch,:,:),lineNoiseOut);
            if showProgress, send(progressQueue,ch); end
        end
    else
        cleaned = cell(1,numel(chans));
        parfor k = 1:numel(chans)
            ch = chans(k);
            cleaned{k} = cleanChannelEpochs(data(ch,:,:),lineNoiseOut);
            if showProgress, send(progressQueue,ch); end
        end
        for k = 1:numel(chans)
            data(chans(k),:,:) = cleaned{k};
        end
    end
else
    for k = 1:numel(chans)
        ch = chans(k);
        data(ch,:,:) = cleanChannelEpochs(data(ch,:,:),lineNoiseOut);
        if showProgress, waitbar(k/numel(chans),progressWindow); end
    end
end
if showProgress, drawnow; end
signal.data = data;
end

function channelData = cleanChannelEpochs(channelData,lineNoise)
% Restart the spectrum, iterative fits and overlap smoothing for each epoch.
for trial = 1:size(channelData,3)
    channelData(1,:,trial) = cleanline_optimized.removeLinesMovingWindow( ...
        channelData(1,:,trial),lineNoise);
end
end
