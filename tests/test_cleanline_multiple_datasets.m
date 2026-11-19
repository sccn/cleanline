function test_cleanline_multiple_datasets
% Run in a fresh MATLAB session with EEGLAB and Parallel Computing Toolbox.
pluginRoot = fileparts(fileparts(mfilename('fullpath')));
eeglabRoot = getenv('CLEANLINE_EEGLAB_ROOT');
assert(isfolder(eeglabRoot),'Set CLEANLINE_EEGLAB_ROOT to the EEGLAB folder.');
addpath(eeglabRoot);
eeglab('nogui');
addpath(genpath(pluginRoot),'-begin');
assert(strcmp(which('pop_cleanline_integrated'),fullfile(pluginRoot,'pop_cleanline_integrated.m')));
rng(20261119,'twister');
quiet = {'ComputeSpectralPower',false,'VerboseOutput',false};
optimized = [{'Implementation','optimized','SignalType','Channels'},quiet];

continuous = [makeDataset('first',4,256,1536,1), ...
    makeDataset('second',4,256,1792,1),makeDataset('third',4,256,2048,1)];
continuous(2).data = double(continuous(2).data);
continuous(2).icaweights = eye(4);
continuous(2).icasphere = eye(4);
continuous(2).icawinv = eye(4);
continuous(2).icaact = continuous(2).data;
continuous(2).icachansind = 1:4;
% The apostrophe in this expression also exercises MATLAB history escaping.
serialArgs = [optimized,{'ExecutionMode','serial','ChanCompIndices','[1;3]'''}];
poolBefore = gcp('nocreate');
[serial,history,report] = pop_cleanline_integrated(continuous,serialArgs{:});
assert(isequal(gcp('nocreate'),poolBefore),'Batch Serial changed the existing pool.');
assert(isequaln(serial,independent(continuous,serialArgs)));
assertBatchMetadata(continuous,serial,[2 4]);
assertBatchReport(report,continuous,'optimized','serial',history);
assertHistory(continuous,history,serial);
column = continuous.';
[columnOut,~,columnReport] = pop_cleanline_integrated(column,serialArgs{:});
assert(isequaln(columnOut,serial.') && isequal(size(columnReport),size(column)));
grid = reshape([continuous continuous(1)],2,2);
[gridOut,~,gridReport] = pop_cleanline_integrated(grid,serialArgs{:});
assert(isequaln(gridOut,independent(grid,serialArgs)) && isequal(size(gridReport),[2 2]));
fprintf('A PASS: three continuous datasets, independent outputs exact (max abs 0); row, column and matrix shapes, history, metadata and unselected channels preserved.\n');

assert(license('test','Distrib_Computing_Toolbox'));
pool = gcp('nocreate');
assert(isempty(pool),'Run in a fresh session so historical tests leave external pools alone.');
pool = parpool('Processes',2);
closePool = onCleanup(@() deleteTestPool(pool));
parallelArgs = [serialArgs,{'ExecutionMode','parallel'}];
[parallelOut,history,report] = pop_cleanline_integrated(continuous,parallelArgs{:});
for k = 1:numel(continuous)
    assert(isequal(parallelOut(k).data,serial(k).data));
end
assert(isequaln(parallelOut,serial) && isequal(gcp('nocreate'),pool));
assertBatchReport(report,continuous,'optimized','parallel',history);
assertHistory(continuous,history,parallelOut);
assert(isequal(gcp('nocreate'),pool));
fprintf('B PASS: Parallel and Serial samples match by isequal (max abs 0); replay reused the same pool.\n');

warningId = 'CleanLine:EpochedDataDiscouraged';
previousWarning = warning('on',warningId);
restoreWarning = onCleanup(@() warning(previousWarning)); %#ok<NASGU>
epochs = [makeDataset('fractional epoch',3,129,201,3),makeDataset('long epoch',3,250,625,3)];
epochArgs = [optimized,{'LineFrequencies',[60 62],'TaperBandWidth',3,'ChanCompIndices',[1 3]}];
[console,epoched,history,report] = evalc('pop_cleanline_integrated(epochs,epochArgs{:})');
assert(count(console,'Applying CleanLine to epoched EEG is strongly discouraged.')==2);
assertBatchReport(report,epochs,'optimized','serial',history);
warning('off',warningId);
assert(isequaln(epoched,independent(epochs,epochArgs)));
assertBatchMetadata(epochs,epoched,2);
assertHistory(epochs,history,epoched);
parallelEpochs = pop_cleanline_integrated(epochs,epochArgs{:},'ExecutionMode','parallel');
assert(all(arrayfun(@(a,b) isequal(a.data,b.data),epoched,parallelEpochs)));
assert(isequaln(epoched,parallelEpochs) && isequal(gcp('nocreate'),pool));
mixed = [continuous(1),epochs(1)];
assert(isequaln(pop_cleanline_integrated(mixed,epochArgs{:}),independent(mixed,epochArgs)));
for k = 1:numel(epochs)
    duration = epochs(k).pnts/epochs(k).srate;
    for trial = 1:epochs(k).trials
        singleTrial = epochs(k);
        singleTrial.data = epochs(k).data(:,:,trial);
        singleTrial.trials = 1;
        out = pop_cleanline_integrated(singleTrial,epochArgs{:}, ...
            'SlidingWinLength',duration,'SlidingWinStep',duration);
        assert(isequal(epoched(k).data(:,:,trial),out.data));
    end
end
changed = epochs;
t = (0:changed(2).pnts-1)/changed(2).srate;
changed(2).data(:,:,2) = changed(2).data(:,:,2) + 20*cos(2*pi*60*t+0.8);
changedOut = pop_cleanline_integrated(changed,epochArgs{:});
assert(isequaln(epoched(1),changedOut(1)));
assert(isequal(epoched(2).data(:,:,[1 3]),changedOut(2).data(:,:,[1 3])));
assert(~isequal(epoched(2).data(:,:,2),changedOut(2).data(:,:,2)));
duration = epochs(1).pnts/epochs(1).srate;
sharedEpochArgs = [epochArgs,{'SlidingWinLength',duration,'SlidingWinStep',duration}];
[sharedEpochs,history] = pop_cleanline_integrated(epochs,sharedEpochArgs{:});
assert(isequaln(sharedEpochs,independent(epochs,sharedEpochArgs)));
assert(isequal(sharedEpochs(1).data,epoched(1).data) && ...
    ~isequal(sharedEpochs(2).data,epoched(2).data));
assertHistory(epochs,history,sharedEpochs);
fprintf('C PASS: different epoch lengths, omitted and shared windows, fractional history, one warning per dataset and exact isolation between datasets and epochs.\n');

heterogeneous = [makeDataset('higher rate',4,256,1536,1),makeDataset('lower rate',2,200,1600,1)];
frequencyArgs = [optimized,{'LineFrequencies',[60 120]}];
[heterogeneousOut,history,report] = pop_cleanline_integrated(heterogeneous,frequencyArgs{:});
assert(isequaln(heterogeneousOut,independent(heterogeneous,frequencyArgs)));
assert(isequal(report(1).effectiveLineFrequencies,[60 120]) && ...
    isequal(report(2).effectiveLineFrequencies,60));
assertBatchMetadata(heterogeneous,heterogeneousOut,[]);
assertHistory(heterogeneous,history,heterogeneousOut);
assertDatasetError(@() pop_cleanline_integrated(heterogeneous,frequencyArgs{:}, ...
    'ChanCompIndices',[1 3]),2,'ChanCompIndices');
assertDatasetError(@() pop_cleanline_integrated(heterogeneous,frequencyArgs{:}, ...
    'ChanCompIndices',[1 1]),1,'ChanCompIndices');
assertDatasetError(@() pop_cleanline_integrated(continuous,'Implementation','optimized'), ...
    2,'SignalType');
assertDatasetError(@() drive_cleanline_test_dialog(@() pop_cleanline_integrated(heterogeneous), ...
    @checkOptimizedDefaults,true),2,'ChanCompIndices');
assertDatasetError(@() drive_cleanline_test_dialog(@() pop_cleanline_integrated(fliplr(epochs)), ...
    @checkOptimizedDefaults,true),2,'SlidingWinLength');
differentNyquist = [heterogeneous(2),makeDataset('low Nyquist',2,100,600,1)];
assertDatasetError(@() drive_cleanline_test_dialog(@() pop_cleanline_integrated(differentNyquist), ...
    @checkOptimizedDefaults,true),2,'LineFrequencies');
diskBacked = heterogeneous;
diskBacked(2).data = 'unloaded.fdt';
assertDatasetError(@() pop_cleanline_integrated(diskBacked,frequencyArgs{:}),2,'load');
fprintf('D PASS: heterogeneous rates, recording lengths and channel counts, per-dataset Nyquist filtering, and indexed errors for incompatible settings and unloaded data.\n');

% The historical route keeps eeg_eval's storage and pool preferences.
clear closePool
eeglab_options;
assert(~option_storedisk,'Run historical fixture tests with EEGLAB disk storage disabled.');
historical = [makeDataset('old first',2,256,1536,1),makeDataset('old second',2,256,1792,1)];
for implementation = ["legacy","original"]
    expected = pop_cleanline(historical,quiet{:},'Newversion',implementation=="original");
    [out,history,report] = pop_cleanline_integrated(historical,quiet{:}, ...
        'Implementation',char(implementation));
    assert(isequaln(expected,out));
    assertBatchReport(report,historical,char(implementation),'official',history);
    assertHistory(historical,history,out);
    historicalGrid = reshape([historical historical],2,2);
    expected = pop_cleanline(historicalGrid(:),quiet{:},'Newversion',implementation=="original");
    out = pop_cleanline_integrated(historicalGrid,quiet{:},'Implementation',char(implementation));
    assert(isequaln(out,reshape(expected,size(historicalGrid))));
    fprintf('E PASS: %s vector and matrix batches match the official historical vector route exactly.\n',implementation);
end

[out,history,report] = drive_cleanline_test_dialog(@() pop_cleanline_integrated(epochs), ...
    @checkOptimizedDefaults,true);
expected = pop_cleanline_integrated(epochs, ...
    'ChanCompIndices',1:3,'SlidingWinLength',duration,'SlidingWinStep',duration);
assert(isequaln(out,expected));
assertBatchReport(report,epochs,'optimized','serial',history);
assertHistory(epochs,history,out);

invoke = @() pop_cleanline_integrated(historical);
[out,history,report] = drive_cleanline_test_dialog(invoke,@selectHistorical,true);
expected = pop_cleanline(historical,quiet{:},'Newversion',false);
assert(isequaln(out,expected));
assertBatchReport(report,historical,'legacy','official',history);
assertHistory(historical,history,out);
[out,history,report] = drive_cleanline_test_dialog(invoke,@checkOptimizedDefaults,false);
assert(isequaln(out,historical) && isempty(history) && isempty(fieldnames(report)));
fprintf('F PASS: one real dialog per invocation, Optimized default and Historical checkbox routes replay exactly, and Cancel preserves the array.\n');
end

function EEG = makeDataset(name,nChannels,Fs,pnts,nTrials)
EEG = eeg_emptyset;
EEG.setname = name;
EEG.nbchan = nChannels;
EEG.pnts = pnts;
EEG.trials = nTrials;
EEG.srate = Fs;
EEG.xmax = (pnts-1)/Fs;
EEG.chanlocs = struct('labels',cellstr("Ch"+(1:nChannels)));
EEG.etc.cleanLineTest = name;
t = (0:pnts-1)/Fs;
data = 0.05*randn(nChannels,pnts,nTrials);
EEG.event = struct('type',{},'latency',{},'epoch',{},'urevent',{});
EEG.urevent = struct('type',{},'latency',{});
for trial = 1:nTrials
    data(:,:,trial) = data(:,:,trial) + (trial+1)*sin(2*pi*60*t + trial*(1:nChannels)'/2);
    latency = (trial-1)*pnts+11;
    EEG.event(trial) = struct('type','cue','latency',latency,'epoch',trial,'urevent',trial);
    EEG.urevent(trial) = struct('type','cue','latency',latency);
end
EEG.data = single(data);
EEG = eeg_checkset(EEG,'eventconsistency');
end

function reference = independent(EEG,args)
reference = EEG;
for k = 1:numel(EEG)
    reference(k) = pop_cleanline_integrated(EEG(k),args{:});
end
end

function assertBatchMetadata(input,output,unselected)
assert(isequal(size(input),size(output)));
assert(isequaln(rmfield(input,'data'),rmfield(output,'data')));
for k = 1:numel(input)
    assert(isequal(size(input(k).data),size(output(k).data)));
    assert(isequal(double(input(k).data(unselected,:,:)),output(k).data(unselected,:,:)));
end
end

function assertBatchReport(report,EEG,implementation,mode,history)
assert(isequal(size(report),size(EEG)));
assert(all(strcmp({report.implementation},implementation)) && all(strcmp({report.executionMode},mode)));
assert(~isempty(history) && all(strcmp({report.historyCommand},history)));
end

function assertHistory(EEG,history,expected)
eval(history);
assert(isequaln(EEG,expected),'History replay must reproduce the complete array.');
end

function assertDatasetError(invoke,index,parameter)
try
    invoke();
catch err
    assert(startsWith(err.message,sprintf('Dataset %d: ',index)) && ...
        contains(err.message,parameter,'IgnoreCase',true),'Unexpected error: %s',err.message);
    return
end
error('Expected an indexed error for %s on dataset %d.',parameter,index);
end

function checkOptimizedDefaults(fig)
checkbox = findobj(fig,'Tag','useoldversion');
execution = findobj(fig,'Tag','executionmode');
assert(get(checkbox,'Value')==0 && get(execution,'Value')==1);
end

function selectHistorical(fig)
set(findobj(fig,'Tag','useoldversion'),'Value',1);
set(findobj(fig,'Tag','computepower'),'Value',0);
set(findobj(fig,'Tag','verb'),'Value',0);
end

function deleteTestPool(pool)
if isvalid(pool)
    delete(pool);
end
end
