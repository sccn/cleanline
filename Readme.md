![Screenshot 2024-08-02 at 09 14 11](https://github.com/user-attachments/assets/4e5be390-d091-4bdf-9ddc-e6fff8b692bb)

# CleanLine

Welcome to the CleanLine plugin for EEGLAB! 

This plugin adaptively estimates and removes sinusoidal (e.g. line) noise from your ICA components
or scalp channels using multi-tapering and a Thompson F-statistic.

CleanLine is written by Tim Mullen (tim@sccn.ucsd.edu) with thanks to Makoto Miyakoshi for beta 
testing. CleanLine makes use of functions modified from the Mitra Lab's Chronux 
Toolbox (www.chronux.org).

CleanLine also makes use of the arg() functionality from Christian Kothe's BCILAB toolbox 
(sccn.ucsd.edu/wiki/BCILAB)

# Instalation

Installation of CleanLine is simple:

1) download CleanLine (if you are reading this, you have probably already completed this step)

2) Unzip the package and copy to <eeglabroot>/plugins/. Start eeglab from the Matlab command line.
   Alternately, you may add CleanLine (with subfolders) to your path and ensure EEGLAB is present 
   in the path.

3) If using the EEGLAB GUI you may start CleanLine from Tools-->CleanLine.
   Alternately, you can start CleanLine from the command line

```matlab
EEGclean = pop_cleanline(EEG);
```

   See "Command-line example" section below for command-line example and other parameters
   
5) For Help, type

```matlab
doc cleanline
```

   Or, hold mouse over any textbox, checkbox, etc in the GUI for tooltip help text.


# Theory

Sinusoidal noise can be a prominent artifact in recorded electrophysiological data. This can stem 
from AC power line fluctuations (e.g. 50/60 Hz line noise + harmonics), power suppliers (e.g. in 
medical equipment), fluorescent lights, etc. Notch filtering is generally undesirable due to 
creation of band-holes, and significant distortion of frequencies around the notch frequency (as well 
as phase distortion at other frequencies and Gibbs rippling in the time-domain). Other approaches for
sinusoidal ("line") noise reduction include adaptive regressive filtering approaches (e.g. RLS, LMS),
but these typically require a reference signal (e.g. externally-recorded noise reference), which is 
often unavailable. Blind-source separation approaches such as ICA may help mitigate line noise, but  
often fail to completely remove the noise due to spatiotemporal non-stationarities in the  noise. 

CleanLine uses an approach for line noise removal advocated by Partha Mitra and Hemant Bokil in 
"Observed Brain Dynamics" (2007), Chapter 7.3.4. 

In brief, the data is traversed by a sliding window. Within each window, the signal is transformed to 
the frequency domain using a multi-taper FFT. The complex amplitude (amplitude and phase) is thus
obtained for each frequency. Under the assumption of a deterministic sinusoid embedded in white 
noise, we can set up a regression of the multi-taper transform (spectrum) of this sinusoidal signal 
onto the multitaper spectrum of the original data at a given frequency. The regression coefficient 
is a complex number representing the complex amplitude (phase and amplitude) of the deterministic 
sinusoid. From this, a time-domain representation of the sinusoid may be constructed and subtracted 
from the data to remove the line.

Typically, one does not know the exact line frequency. For instance, in the U.S.A., power line noise 
is not guaranteed to be at exactly 60 Hz (or even to have constant phase over a given period of time).
To ameliorate this problem a Thompson F-Test may be applied to determine the statistical significance
of a non-zero coefficient in the above regression (indicating a sinusoid with significantly non-zero 
amplitude). We can then search within a narrow band around the expected location of the line for the 
frequency which maximizes this F-statistic above a significance threshold (e.g. p=0.05).

Line frequency scanning can be enabled/disabled using the 'ScanForLines' option.

Overlapping short (e.g. 2-4 second) windows can be used to adaptively estimate the frequency, phase,
and amplitude of the sinusoidal noise components (which typically change over the course of a recording 
session). The discontinuity at the point of window overlap can be smoothed using a sigmoidal function.
The example below demonstrates such a function for different smoothing factors (slope of the sigmoid).

```matlab
winsize = 4;                        % window length in seconds
winstep = 2;                        % window step in seconds (50% overlap)
Fs      = 128;                      % sampling rate
tau     = [1 10 100]';              % smoothing factors
overlap = winsize-winstep;
toverlap = -overlap/2:(1/Fs):overlap/2;

% specify the smoothing function
foverlap = 1-1./(1+exp(-repmat(tau,1,length(toverlap)).*repmat(toverlap,length(tau),1)/overlap));

% define some colours
yellow  = [255, 255, 25]/255;
red     = [255 0 0]/255;
h       = zeros(1,3+length(tau));

% plot the figure
figure;
axis([-winsize+overlap/2 winsize-overlap/2 0 1]); set(gca,'ColorOrder',[0 0 0; 0.7 0 0.8; 0 0 1],'fontsize',11);
hold on
h(1)=hlp_vrect([-winsize+overlap/2 -overlap/2], 'yscale',[0 1],'patchProperties',{'FaceColor',yellow,        'FaceAlpha',1,'EdgeColor','none','EdgeAlpha',0.5}); 
h(2)=hlp_vrect([overlap/2 winsize-overlap/2],   'yscale',[0 1],'patchProperties',{'FaceColor',red,           'FaceAlpha',1,'EdgeColor','none','EdgeAlpha',0.5});
h(3)=hlp_vrect([-overlap/2 overlap/2],          'yscale',[0 1],'patchProperties',{'FaceColor',(yellow+red)/2,'FaceAlpha',1,'EdgeColor','none','EdgeAlpha',0.5});
h(4:end) = plot(toverlap,foverlap,'linewidth',2);
plot(toverlap,1-foverlap,'linewidth',1,'linestyle','--');
hold off;
xlabel('Time (sec)'); ylabel('Smoothing weight'); 
title('Plot of window overlap smoothing function vs. time for different smoothing factors');
legend(h,[{'Window 1','Window 2','Overlap'},cellstr(num2str(tau))']);
```

The smoothing factor is determined by the 'SmoothingFactor' parameter in cleanline().

CleanLine allows you to specify the multi-taper frequency resolution by the 'Bandwidth' parameter. 
This is the width of a peak in the spectrum for a sinusoid at given frequency. Due to the time-frequency 
uncertainty principle, decreasing bandwidth increases the necessary length of the sliding window in 
order to obtain a reasonable frequency decomposition. The number of tapers, K, used by CleanLine is 
given by K=2TW-1 where T is the temporal resolution in seconds (window length) and W is the frequency 
resolution (Bandwidth) in Hz. CleanLine fixes T to be the sliding window length, so W (Bandwidth) is the
only required parameter. If the 'Verbosity' option is set to 1, then CleanLine will display the multi-taper
parameters in the command line on execution of the function.


# Tips on running CleanLine

The default options should work quite well, but parameters may need to be tweaked depending on the setup. 
If you have multiple epochs you need to make sure that your window size and step size exactly divides the 
epoch length. In other words, you do not want any sliding windows to overlap two epochs since line noise 
phase and amplitude may shift at that point, making it impossible to perform the time-domain line subtraction.
If you have relatively short epochs (e.g. < 5 sec) it is best if each window is taken to be the length
of the epoch and the step size is equal to the window length. In this way, the lines are estimated and 
removed for each epoch individually. When using the GUI, the default values for window length and step size
are automatically set to the epoch length.

If cleaning continuous, un-epoched data, then you may wish to use sliding windows of 3-4 seconds with 50% 
overlap.

If using the EEGLAB GUI, commands are stored in EEG.history so that the eegh() command will return the 
command-line function call corresponding to the last GUI execution of CleanLine.

You might find it useful to try the option ('PlotFigures',true) on a subset of channels/components to get a
sense of the performance of CleanLine for difference parameter choices (and also to identify where the most
significant lines lie in the spectrum) and then, once you are satisfied with the parameters, turn this option
off before cleaning all the remaining channels/components. NOTE: CleanLine is *considerably slower* if PlotFigures 
is enabled. If you don't care to see the visualize the final results of the cleaning operation, you may also
wish to set ('ComputeSpectralPower',false) which will speed up computation considerably.


# Typical command-line example (Multiple trials)

```matlab
% This will run cleanline on all channels, scanning for lines +/- 1 Hz around the 60 and 120 Hz frequencies. 
% Each epoch will be cleaned individually and epochs containing lines that are significantly sinusoidal at 
% the p<=0.01 level will be cleaned. 

EEG = pop_cleanline(EEG, 'Bandwidth',2,'ChanCompIndices',[1:EEG.nbchan] ,'SignalType','Channels','ComputeSpectralPower',true,'LineFrequencies',[60 120] ,'NormalizeSpectrum',false,'LineAlpha',0.01,'PaddingFactor',2,'PlotFigures',false,'ScanForLines',true,'SmoothingFactor',100,'VerbosityLevel',1,'SlidingWinLength',EEG.pnts/EEG.srate,'SlidingWinStep',EEG.pnts/EEG.srate);
```

# Help for cleanline.m (These can also be passed to pop_cleanline)

## Function input
  **LineFrequencies:**      Line noise frequencies to remove                                                                      
                        Input Range  : Unrestricted                                                                           
                        Default value: 60  120                                                                                
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **ScanForLines:**         Scan for line noise                                                                                   
                        This will scan for the exact line frequency in a narrow range around the specified LineFrequencies    
                        Input Range  : [true false]                                                                           
                        Default value: true                                                                                      
                        Input Data Type: boolean                                                                              
                                                                                                                              
  **LineAlpha:**            p-value for detection of significant sinusoid                                                                        
                        Input Range  : [0  1]                                                                                 
                        Default value: 0.01                                                                                   
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **Bandwidth:**            Bandwidth (Hz)                                                                                        
                        This is the width of a spectral peak for a sinusoid at fixed frequency. As such, this defines the     
                        multi-taper frequency resolution.                                                                     
                        Input Range  : Unrestricted                                                                           
                        Default value: 1                                                                                      
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **SignalType:**          Type of signal to clean                                                                               
                        Cleaned ICA components will be backprojected to channels. If channels are cleaned, ICA activations    
                        are reconstructed based on clean channels.                                                            
                        Possible values: 'Components','Channels'                                                              
                        Default value  : 'Components'                                                                         
                        Input Data Type: string                                                                               
                                                                                                                              
  **ChanCompIndices:**      IDs of Chans/Comps to clean                                                                           
                        Input Range  : Unrestricted                                                                           
                        Default value: 1:152                                                                                  
                        Input Data Type: any evaluable Matlab expression.                                                     
                                                                                                                              
  **SlidingWinLength:**     Sliding window length (sec)                                                                           
                        Default is the epoch length.                                                                          
                        Input Range  : [0  4]                                                                                 
                        Default value: 4                                                                                      
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **SlidingWinStep:**       Sliding window step size (sec)                                                                        
                        This determines the amount of overlap between sliding windows. Default is window length (no           
                        overlap).                                                                                             
                        Input Range  : [0  4]                                                                                 
                        Default value: 4                                                                                      
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **SmoothingFactor:**      Window overlap smoothing factor                                                                       
                        A value of 1 means (nearly) linear smoothing between adjacent sliding windows. A value of Inf means   
                        no smoothing. Intermediate values produce sigmoidal smoothing between adjacent windows.               
                        Input Range  : [1  Inf]                                                                               
                        Default value: 100                                                                                    
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **PaddingFactor:**        FFT padding factor                                                                                    
                        Signal will be zero-padded to the desired power of two greater than the sliding window length. The    
                        formula is NFFT = 2^nextpow2(SlidingWinLen*(PadFactor+1)). e.g. For SlidingWinLen = 500, if PadFactor = -1,    
                        we do not pad; if PadFactor = 0, we pad the FFT to 512 points, if PadFactor=1, we pad to 1024 points etc.                                                                                                  
                        Input Range  : [-1  Inf]                                                                              
                        Default value: 2                                                                                      
                        Input Data Type: real number (double)                                                                 
                                                                                                                              
  **ComputeSpectralPower:** Visualize Original and Cleaned Spectra                                                                
                        Original and clean spectral power will be computed and visualized at end of processing                             
                        Input Range  : [true false]                                                                           
                        Default value: true                                                                                      
                        Input Data Type: boolean                                                                              
                                                                                                                              
  **NormalizeSpectrum:**    Normalize log spectrum by detrending (not generally recommended)                                                                     
                        Input Range  : [true false]                                                                            
                        Default value: false                                                                                      
                        Input Data Type: boolean                                                                              
                                                                                                                              
  **VerboseOutput:**        Produce verbose output                                                                                
                        Input Range  : [true false]                                                                           
                        Default value: true                                                                                      
                        Input Data Type: boolean                                                                
                                                                                                                              
  **PlotFigures:**          Plot Individual Figures                                                                               
                        This will generate figures of F-statistic, spectrum, etc for each channel/comp while processing       
                        Input Range  : [true false]                                                                            
                        Default value: false                                                                                      
                        Input Data Type: boolean  
 
##  Output                Information
  **EEG**                   Cleaned EEG dataset
  
  **Sorig**                 Original multitaper spectrum for each component/channel
  
  **Sclean**                Cleaned multitaper spectrum for each component/channel
  
  **f**                     Frequencies at which spectrum is estimated in Sorig, Sclean
  
  **amps**                  Complex amplitudes of sinusoidal lines for each window (line time-series for window i can be reconstructed by creating a sinuoid with frequency f{i} and complex
amplitude amps{i})

  **freqs**                 Exact frequencies at which lines were removed for each window (cell array)
  
  **g**                     Parameter structure. Function call can be replicated exactly by calling >> cleanline(EEG,g);


# CleanLine versions
V1 - Original version by Tim Mullen

V2 - Include a rewrite by Kay Robbins with integration by Arnaud Delorme and 
testing by Makoto Miyakoshi

V2.1 - Fix menu returning automatically in EEGLAB 2024.0 and 2024.1

---

## CleanLine V2 optimization and validation (November 2026)

The [CleanLine GPT-6 Optimized project](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized), led by Makoto Miyakoshi with GPT-6 programming assistance, provides a computationally optimized implementation of **Kay Robbins' CleanLine V2** (`newversion=1`). The optimization does **not** change the line-noise removal method or claim a new signal-processing algorithm.

### Changes

- Compute the fixed F-test critical value (`finv`) once per invocation instead of repeatedly within the moving-window calculations.
- Reuse the frequency grid, selected frequency-bin indices, and taper sums across sliding windows. The detection threshold, harmonic selection, iteration limit, sinusoid fitting, and overlap smoothing are unchanged.
- Provide an explicit Optimized V2 Serial/Parallel choice in the official plugin's single CleanLine settings dialog. Historical V1 and Original V2 remain selectable.

**Attribution:** Channel-level `parfor` was already present in [Kay Robbins' PREP implementation in 2014](https://github.com/VisLab/EEG-Clean-Tools/commit/a37f73be8675b6670df9c749ba1a9eeeaa4b0b96) and was [integrated into official CleanLine V2 in 2020](https://github.com/sccn/cleanline/commit/294a6178157cd77fec25fcbb7e7ad2175076a60f) by Arnaud Delorme. Parallel execution is **not** a new algorithmic optimization in this update.

### Execution-time measurements

Measured on the same **33-channel, 256-Hz EEG** recording with MATLAB R2024b Update 3, EEGLAB 2024.0, Windows 11, and an Intel Core i9-10900K. Parallel tests used a pre-started pool of **eight workers**; the measured **12.08-second pool startup** was excluded. All values below are elapsed processing times, excluding data loading and saving.

| Recording | Original serial | Optimized serial | Original parallel | Optimized parallel |
| --- | ---: | ---: | ---: | ---: |
| 5 min | 51.47 s | 40.00 s | 15.28 s | 13.69 s |
| 30 min | 289.94 s | 226.85 s | 80.43 s | 69.94 s |

The **source-code changes alone** improved matched Serial performance by **1.28–1.29×** and matched Parallel performance by **1.12–1.15×**. Comparing Original Serial with Optimized Parallel yields an **overall 3.76× speedup at 5 min and 4.15× at 30 min**, including the pre-existing benefit of parallel processing. The 5- and 30-minute values are medians of three measured runs per condition.

### Numerical validation

The following table reproduces the output-difference results reported in the [original optimization README](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized#numerical-validation). The reference is **Original CleanLine 2.00 V2** (`newversion=1`), using the same EEG data, settings, and execution mode for each matched comparison.

| Comparison against Original V2 | Maximum absolute error | Relative L2 error | Equality |
| --- | ---: | ---: | --- |
| 33-channel cleaned EEG, all 5–30 min parallel runs | 0 | 0 | Exact sample equality |
| 33-channel removed line noise, all 5–30 min parallel runs | 0 | 0 | Exact sample equality |
| 33-channel cleaned EEG, all 5–30 min serial runs | 0 | 0 | Exact sample equality |
| 33-channel removed line noise, all 5–30 min serial runs | 0 | 0 | Exact sample equality |

The separately tested standalone candidate was compared with CleanLine **2.00 V2** using identical input EEG and processing parameters. Over the tested **5-, 10-, 15-, 20-, 25-, and 30-minute** recording lengths, the optimized and original implementations produced **exactly equal cleaned samples and removed line-noise samples** in both Serial and Parallel modes (**maximum absolute error = 0; relative L2 error = 0**). Direct comparisons with the installed public CleanLine V2 on 33-channel **5- and 30-minute** inputs likewise reported exact output equality and preserved event metadata. Five deterministic synthetic test cases also matched exactly in the reported waveform and detection quantities.

These are results from the **validated standalone V2 candidate**; they do not establish exhaustive equivalence for all data and settings, nor do they substitute for testing the final upstream-integrated plugin, including historical V1 and ICA-component functionality.

- [Complete optimization report and benchmark methodology](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized#benchmark-method)
- [Full numerical validation report](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized/blob/main/docs/VALIDATION_REPORT_20261118.md)
- [Runnable original-vs-optimized / installed-V2 comparison test](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized/blob/main/tests/test_cleanline_gpt6_full_parity.m)
- [Raw timing measurements](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized/blob/main/benchmarks/raw_runs.csv)

### Official `v2-optimized` integration

The official branch now registers one **Tools → CleanLine** menu item and
offers Historical V1, Original V2, and Optimized V2 in a single settings
window. The original `pop_cleanline` and `cleanline` entry points, command
arguments, V1 ICA component cleaning, V2 implementation, copyright notices,
and `license.txt` remain available. The optimized engine uses the plugin's
existing Chronux helpers. V2 cleans continuous channel data; select V1 for
ICA components. Parallel execution of Optimized V2 requires Parallel
Computing Toolbox.

The new scriptable entry point is:

```matlab
[EEGclean, com, report] = pop_cleanline_integrated(EEG, ...
    'Implementation','optimized', 'ExecutionMode','serial', ...
    'LineFrequencies',60, 'ChanCompIndices',1:EEG.nbchan, ...
    'SlidingWinLength',4, 'SlidingWinStep',1, ...
    'Bandwidth',2, 'TaperBandWidth',2, 'LineAlpha',0.01, ...
    'PaddingFactor',2, 'SmoothingFactor',100);
```

Existing calls to `pop_cleanline(EEG,...)` and their EEGLAB history remain
compatible. The [new integration results](docs/INTEGRATION_VALIDATION.md)
separate this branch's tests from the standalone measurements above.

Set `CLEANLINE_EEGLAB_ROOT` to your EEGLAB installation and run
`run('tests/test_cleanline_integration.m')` from this plugin folder. For the
adapted 5- and 30-minute official V2 comparison, set `CLEANLINE_TEST_SET` to
a continuous `.set` file with at least 33 channels and 30 minutes of data,
then run `run('tests/test_cleanline_gpt6_full_parity.m')`. The dataset is not
included in this repository.
