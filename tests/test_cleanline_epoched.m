function test_cleanline_epoched
% Run with CLEANLINE_EEGLAB_ROOT set and Parallel Computing Toolbox available.
pluginRoot = fileparts(fileparts(mfilename('fullpath')));
eeglabRoot = getenv('CLEANLINE_EEGLAB_ROOT');
assert(isfolder(eeglabRoot),'Set CLEANLINE_EEGLAB_ROOT to the EEGLAB folder.');
addpath(eeglabRoot);
eeglab('nogui');
addpath(genpath(pluginRoot),'-begin');
assert(strcmp(which('pop_cleanline'),fullfile(pluginRoot,'pop_cleanline.m')));
rng(263284,'twister');

args = {'SignalType','Channels','LineFrequencies',[60 62], ...
    'ComputeSpectralPower',false,'VerboseOutput',false};
warningId = 'CleanLine:EpochedDataDiscouraged';
expectedWarning = ['Applying CleanLine to epoched EEG is strongly discouraged. ' ...
    'Short epochs limit the reliability of multitaper spectral estimation and sinusoidal fitting, ' ...
    'potentially resulting in inconsistent line-noise removal across trials. ' ...
    'Whenever possible, remove line noise from continuous EEG before epoching. ' ...
    'Processing will continue with epoch-independent estimation.'];
previousWarning = warning('on',warningId);
restoreWarning = onCleanup(@() warning(previousWarning)); %#ok<NASGU>
EEG = makeDataset(250,500,'single');
poolBeforeSerial = gcp('nocreate');
[warningOutput,optimized] = evalc('pop_cleanline_integrated(EEG,args{:})');
assertEpochedWarning(warningOutput,expectedWarning);
assert(isequal(gcp('nocreate'),poolBeforeSerial),'Serial changed the current pool.');
warning('off',warningId);
options = [args,{'SlidingWinLength',2,'SlidingWinStep',2}];
[silentOutput,explicit] = evalc('pop_cleanline_integrated(EEG,options{:})');
assert(~contains(silentOutput,'Warning:'),'Disabling the warning should silence it.');
assert(isequal(optimized.data,explicit.data),'The epoched defaults should use the whole epoch.');
fprintf('D PASS: omitted window settings match an explicit 2 s length and step.\n');

pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('Processes',2);
    closePool = onCleanup(@() delete(pool)); %#ok<NASGU>
end

warning('on',warningId);
parallelOptions = [options,{'ExecutionMode','parallel'}];
[warningOutput,parallelOutput] = evalc('pop_cleanline_integrated(EEG,parallelOptions{:})');
assertEpochedWarning(warningOutput,expectedWarning);
assert(isequal(optimized.data,parallelOutput.data));
assert(isequal(gcp('nocreate'),pool));
fprintf('Warning PASS: one exact warning for each wrapper call in Serial and Parallel. Disabling it preserves identical output.\n');
warning('off',warningId);

% Epoch duration, window length and step, in seconds.
configurations = [2 2 2; 2 1.6 0.3; 4 4 4; 4 3 0.75; 6 6 6; 6 4 1.25];
originalDifference = false;
for inputClass = ["single","double"]
    chans = 1:2;
    if inputClass == "double"
        chans = [1 3];
    end
    for row = 1:size(configurations,1)
        duration = configurations(row,1);
        win = configurations(row,2);
        step = configurations(row,3);
        EEG = makeDataset(250,250*duration,inputClass);
        options = [args,{'ChanCompIndices',chans,'TaperBandWidth',3, ...
            'SlidingWinLength',win,'SlidingWinStep',step}];
        optimized = pop_cleanline_integrated(EEG,options{:});
        assert(~isequal(optimized.data,double(EEG.data)),'Test signal was unchanged.');
        independent = double(EEG.data);
        for trial = 1:EEG.trials
            epoch = pop_select(EEG,'trial',trial);
            % EEGLAB's single-precision preference can cast the selected data.
            epoch.data = EEG.data(:,:,trial);
            cleaned = pop_cleanline_integrated(epoch,options{:});
            assert(isequal(optimized.data(:,:,trial),cleaned.data),'Epoch differs from an independent call.');
            independent(:,:,trial) = cleaned.data;
        end

        changed = EEG;
        t = (0:EEG.pnts-1)/EEG.srate;
        changed.data(:,:,2) = changed.data(:,:,2) + 35*sin(2*pi*60*t+2.5);
        changedOutput = pop_cleanline_integrated(changed,options{:});
        assert(isequal(optimized.data(:,:,[1 3]),changedOutput.data(:,:,[1 3])), ...
            'Changing epoch 2 affected another epoch.');

        parallelOutput = pop_cleanline_integrated(EEG,options{:},'ExecutionMode','parallel');
        assert(isequal(optimized.data,parallelOutput.data),'Serial and Parallel differ.');
        assert(isequal(gcp('nocreate'),pool),'Parallel replaced the existing pool.');

        unselected = setdiff(1:EEG.nbchan,chans);
        assert(isequal(size(optimized.data),size(EEG.data)));
        assert(isequaln(rmfield(optimized,'data'),rmfield(EEG,'data')));
        assert(isa(optimized.data,'double'),'CleanLine V2 should return double.');
        assert(isequal(optimized.data(unselected,:,:),double(EEG.data(unselected,:,:))));
        fprintf(['B-F PASS: %s, epoch %g s, window %g s, step %g s, ' ...
            'channels %s; independent, contamination and parallel max abs errors all 0.\n'], ...
            inputClass,duration,win,step,mat2str(chans));

        if inputClass == "single"
            originalOptions = [options,{'Newversion',true}];
            original = pop_cleanline(EEG,originalOptions{:});
            originalIndependent = double(EEG.data);
            for trial = 1:EEG.trials
                epoch = pop_select(EEG,'trial',trial);
                cleaned = pop_cleanline(epoch,originalOptions{:});
                originalIndependent(:,:,trial) = cleaned.data;
            end
            originalDelta = max(abs(original.data - originalIndependent),[],'all');
            optimizedDelta = max(abs(original.data - optimized.data),[],'all');
            independentDelta = max(abs(originalIndependent - independent),[],'all');
            assert(isequal(originalIndependent,independent),'Independent Original and Optimized calls differ.');
            originalDifference = originalDifference || originalDelta > 0;
            fprintf(['H PASS: epoch %g s, window %g s, step %g s; max abs ' ...
                'Original full vs independent %.17g, Original full vs Optimized full %.17g, ' ...
                'Original independent vs Optimized full %.17g.\n'], ...
                duration,win,step,originalDelta,optimizedDelta,independentDelta);
        end
    end
end
assert(originalDifference,'Expected a difference from concatenated Original V2.');

% At this rate, samples/Fs is one ULP larger than samples*(1/Fs).
EEG = makeDataset(129,201,'single');
win = EEG.pnts/EEG.srate;
assert(win > EEG.pnts*(1/EEG.srate));
signal = struct('data',EEG.data,'srate',EEG.srate);
lineNoise = struct('lineNoiseChannels',1:EEG.nbchan,'Fs',EEG.srate, ...
    'lineFrequencies',[60 62],'p',0.01,'fScanBandWidth',2,'taperBandWidth',2, ...
    'taperWindowSize',win,'taperWindowStep',win,'tau',100,'pad',2, ...
    'fPassBand',[0 EEG.srate/2],'maximumIterations',10);
warning('on',warningId);
[warningOutput,direct] = evalc('cleanline_optimized.cleanLineNoise(signal,lineNoise)');
assertEpochedWarning(warningOutput,expectedWarning);
[warningOutput,optimized] = evalc('pop_cleanline_integrated(EEG,args{:})');
assertEpochedWarning(warningOutput,expectedWarning);
assert(all(isfinite(direct.data),'all'));
assert(isequal(direct.data,optimized.data));
for trial = 1:EEG.trials
    epoch = struct('data',signal.data(:,:,trial),'srate',signal.srate);
    [continuousLog,cleaned] = evalc('cleanline_optimized.cleanLineNoise(epoch,lineNoise)');
    assert(~contains(continuousLog,'Warning:'),'Continuous direct call emitted a warning.');
    assert(isequal(direct.data(:,:,trial),cleaned.data));
end
epoch = pop_select(EEG,'trial',1);
options = [args,{'SlidingWinLength',win,'SlidingWinStep',win}];
[continuousLog,cleaned] = evalc('pop_cleanline_integrated(epoch,options{:})');
assert(~contains(continuousLog,'Warning:'),'Continuous wrapper call emitted a warning.');
assert(isequal(direct.data(:,:,1),cleaned.data));
fprintf('B/D PASS: direct call on 201-sample epochs at 129 Hz matched independent calls exactly.\n');
fprintf('Warning PASS: one exact warning on direct epoched input, silence on continuous input.\n');

% Invoke the real dialog and accept its epoch-length defaults.
[guiOutput,history,report,warningOutput] = drive_cleanline_test_dialog( ...
    @() pop_cleanline_integrated(EEG),@(fig) checkEpochDefaults(fig,win),true);
assertEpochedWarning(warningOutput,expectedWarning);
assert(~isempty(history) && strcmp(report.implementation,'optimized') ...
    && strcmp(report.executionMode,'serial'));
warning('off',warningId);
explicit = pop_cleanline_integrated(EEG,'Implementation','optimized','LineFrequencies',60, ...
    'SlidingWinLength',win,'SlidingWinStep',win);
assert(isequal(guiOutput.data,explicit.data));
input = EEG;
eval(history);
assert(isequal(EEG.data,guiOutput.data),'Fractional epoch GUI history differs.');
EEG = input;
fprintf('Warning/D PASS: dialog defaults accepted, one exact warning, output identical to explicit settings.\n');

tooShort = makeDataset(250,250,'single');
assertError(@() pop_cleanline_integrated(tooShort,args{:}),'CleanLine:EpochTapers');
invalid = lineNoise;
invalid.taperWindowSize = win + 1/EEG.srate;
assertError(@() cleanline_optimized.cleanLineNoise(signal,invalid),'CleanLine:EpochWindowLength');
invalid = lineNoise;
invalid.taperWindowStep = 0;
assertError(@() cleanline_optimized.cleanLineNoise(signal,invalid),'CleanLine:EpochWindowStep');
fprintf('D PASS: invalid epoch taper, window length and step produce specific errors.\n');
end

function EEG = makeDataset(Fs,pnts,inputClass)
EEG = eeg_emptyset;
EEG.srate = Fs;
EEG.pnts = pnts;
EEG.nbchan = 4;
EEG.trials = 3;
EEG.chanlocs = struct('labels',{'A','B','C','D'});
EEG.etc.cleanLineTestMarker = 'metadata';
EEG.event = struct('type',{},'epoch',{},'latency',{},'urevent',{});
EEG.urevent = struct('type',{},'latency',{});
t = (0:pnts-1)/Fs;
EEG.xmax = t(end);
data = 0.05*randn(EEG.nbchan,pnts,EEG.trials);
for trial = 1:EEG.trials
    data(:,:,trial) = data(:,:,trial) + ...
        (trial+0.5)*sin(2*pi*60*t + trial*(1:EEG.nbchan)'/2);
    EEG.event(trial) = struct('type','stimulus','epoch',trial, ...
        'latency',(trial-1)*pnts+11,'urevent',trial);
    EEG.urevent(trial) = struct('type','stimulus','latency',(trial-1)*pnts+11);
end
EEG.data = data;
EEG = eeg_checkset(EEG,'eventconsistency');
EEG.data = cast(data,inputClass);
assert(numel(EEG.event)==EEG.trials && numel(EEG.urevent)==EEG.trials ...
    && numel(EEG.epoch)==EEG.trials);
end

function assertEpochedWarning(output,expected)
assert(contains(regexprep(output,'\s+',' '),['Warning: ' expected]), ...
    'Expected the full epoched warning. Captured output:\n%s',output);
assert(count(output,'Applying CleanLine to epoched EEG is strongly discouraged.') == 1, ...
    'Expected one epoched warning per dataset invocation.');
end

function checkEpochDefaults(fig,epochSeconds)
win = str2double(get(findobj(fig,'Tag','winsize'),'String'));
step = str2double(get(findobj(fig,'Tag','winstep'),'String'));
oldVersion = get(findobj(fig,'Tag','useoldversion'),'Value');
displayedSeconds = str2double(hlp_tostring(epochSeconds));
assert(win==displayedSeconds && step==displayedSeconds && oldVersion==0);
end

function assertError(invoke,identifier)
try
    invoke();
catch err
    assert(strcmp(err.identifier,identifier),'Expected %s, received %s: %s', ...
        identifier,err.identifier,err.message);
    return
end
error('Expected error %s.',identifier);
end
