function [EEGclean, commstr, report] = pop_cleanline_integrated(EEG, varargin)
%POP_CLEANLINE_INTEGRATED CleanLine settings dialog and scriptable engine selection.
% Original pop_cleanline and its command history remain supported unchanged.
% Optimized V2 accepts EEGLAB dataset arrays. Omitted options are resolved
% for each dataset; settings accepted in the dialog are shared by the array.
% Usage: [EEG,com] = pop_cleanline_integrated(EEG, ...
%     'Implementation','optimized','ExecutionMode','serial', ...
%     'LineFrequencies',60);
EEGclean = EEG;
commstr = '';
report = struct();
isBatch = numel(EEG) > 1;
fromDialog = isempty(varargin);
if fromDialog
    [varargin, accepted] = cleanline_integrated_dialog(EEG(1));
    if ~accepted, return; end
end
assert(mod(numel(varargin),2)==0,'Options must be name/value pairs.');
implementation = 'optimized';
executionMode = 'serial';
showProgress = false;
keep = true(1,numel(varargin));
for k=1:2:numel(varargin)
    switch lower(char(varargin{k}))
        case 'implementation'
            implementation = lower(char(varargin{k+1})); keep(k:k+1)=false;
        case 'executionmode'
            executionMode = lower(char(varargin{k+1})); keep(k:k+1)=false;
        case 'showprogress'
            showProgress = logical(varargin{k+1}); keep(k:k+1)=false;
    end
end
args = varargin(keep);
assert(ismember(implementation,{'legacy','original','optimized'}), ...
    'Implementation must be legacy, original, or optimized.');
assert(ismember(executionMode,{'serial','parallel'}), ...
    'ExecutionMode must be serial or parallel.');
if ~strcmp(implementation,'optimized')
    useOriginalV2 = strcmp(implementation,'original');
    if isBatch && fromDialog
        for k = 1:numel(EEG)
            try
                [~,~,~,~,~,~,g] = cleanline_check('EEG',EEG(k),args{:},'Newversion',useOriginalV2);
                validateBatchSettings(EEG(k),g,false,true);
            catch err
                raiseDatasetError(k,err);
            end
        end
    end
    % eeg_eval uses length and expects a vector. Preserve array shape for callers.
    historicalInput = EEG;
    if isBatch && ~isvector(EEG)
        historicalInput = EEG(:);
    end
    [EEGclean,commstr] = pop_cleanline(historicalInput,args{:},'Newversion',useOriginalV2);
    if fromDialog || isBatch
        commstr = historyCommand(implementation,executionMode,args);
    end
    EEGclean = reshape(EEGclean,size(EEG));
    report = repmat(struct('implementation',implementation, ...
        'executionMode','official','historyCommand',commstr),size(EEG));
    return
end
if strcmp(executionMode,'parallel')
    assert(license('test','Distrib_Computing_Toolbox'), ...
        'Parallel mode requires Parallel Computing Toolbox.');
end

reports = cell(size(EEG));
for k = 1:numel(EEG)
    try
        [EEGclean(k),reports{k}] = runOptimized( ...
            EEG(k),args,executionMode,showProgress,isBatch,fromDialog);
    catch err
        if ~isBatch
            rethrow(err);
        end
        raiseDatasetError(k,err);
    end
end
report = reshape([reports{:}],size(EEG));
historyArgs = args;
if isBatch && showProgress
    historyArgs = [historyArgs {'ShowProgress',true}];
end
commstr = historyCommand(implementation,executionMode,historyArgs);
[report.historyCommand] = deal(commstr);
end

function [EEGclean,report] = runOptimized(EEG,args,executionMode,showProgress,isBatch,fromDialog)
% Reuse the official argument definitions and validation for the V2 options.
[~,~,~,~,~,~,g] = cleanline_check('EEG',EEG,args{:},'Newversion',true);
if isBatch
    validateBatchSettings(EEG,g,true,fromDialog);
end
assert(strcmpi(g.sigtype,'Channels'), ...
    'V2 cleans channels; use Legacy V1 for ICA components.');
assert(g.winstep<=g.winsize,'SlidingWinStep cannot exceed SlidingWinLength.');
linefreqs = g.linefreqs;
if numel(linefreqs)==1, linefreqs=linefreqs*(1:4); end
linefreqs(linefreqs>EEG.srate/2)=[];
lineNoise = struct('lineNoiseMethod','clean', ...
    'lineNoiseChannels',g.chanlist,'Fs',EEG.srate, ...
    'lineFrequencies',linefreqs,'p',g.p, ...
    'fScanBandWidth',g.bandwidth,'taperBandWidth',g.taperbandwidth, ...
    'taperWindowSize',g.winsize,'taperWindowStep',g.winstep, ...
    'tau',g.tau,'pad',g.pad,'fPassBand',[0 EEG.srate/2], ...
    'maximumIterations',10,'executionMode',executionMode, ...
    'showProgress',showProgress);
if strcmp(executionMode,'parallel')
    pool = gcp('nocreate');
    if isempty(pool)
        try
            parpool;
        catch ME
            error('CleanLine:ParallelPool', ...
                'Could not start MATLAB parallel pool: %s',ME.message);
        end
    end
end
clock = tic;
[EEGclean,effective] = cleanline_optimized.cleanLineNoise(EEG,lineNoise);
report = struct('implementation','optimized', ...
    'executionMode',executionMode,'effectiveLineFrequencies', ...
    effective.lineFrequencies,'elapsedCoreSeconds',toc(clock));
end

function validateBatchSettings(EEG,g,optimized,shared)
if optimized && ~strcmpi(g.sigtype,'Channels')
    error('CleanLine:BatchSignalType', ...
        'SignalType must be Channels for Optimized V2. Use Historical V1 for ICA components.');
end
if optimized && (~isnumeric(EEG.data) || isempty(EEG.data))
    error('CleanLine:BatchDataNotLoaded', ...
        'Load this dataset''s samples into memory before running Optimized V2.');
end

indexLimit = EEG.nbchan;
if strcmpi(g.sigtype,'Components')
    indexLimit = size(EEG.icawinv,2);
    if indexLimit == 0
        error('CleanLine:BatchSignalType','SignalType Components requires an ICA decomposition.');
    end
end
if ~isempty(g.chanlist)
    validateattributes(g.chanlist,{'numeric'},{'vector','real','integer','>=',1,'<=',indexLimit}, ...
        mfilename,'ChanCompIndices');
    if optimized && numel(unique(g.chanlist)) ~= numel(g.chanlist)
        error('CleanLine:BatchChannelIndices','ChanCompIndices must list each channel at most once.');
    end
end

validateattributes(g.winsize,{'numeric'},{'scalar','real','finite','positive'},mfilename,'SlidingWinLength');
validateattributes(g.winstep,{'numeric'},{'scalar','real','finite','positive'},mfilename,'SlidingWinStep');
windowSamples = round(g.winsize*EEG.srate);
stepSamples = round(g.winstep*EEG.srate);
if windowSamples < 1 || windowSamples > EEG.pnts
    error('CleanLine:BatchWindowLength', ...
        'SlidingWinLength must span 1 to %d samples at %.17g Hz (one epoch or recording).', ...
        EEG.pnts,EEG.srate);
end
if stepSamples < 1 || stepSamples > EEG.pnts || (optimized && g.winstep > g.winsize)
    error('CleanLine:BatchWindowStep', ...
        'SlidingWinStep must span at least one sample, fit the dataset, and for Optimized V2 be no longer than SlidingWinLength.');
end
if shared && ~any(g.linefreqs(:) < EEG.srate/2)
    error('CleanLine:BatchLineFrequencies', ...
        'LineFrequencies must include a frequency below this dataset''s Nyquist limit (%.17g Hz).',EEG.srate/2);
end
end

function raiseDatasetError(index,cause)
identifier = cause.identifier;
if isempty(identifier)
    identifier = 'CleanLine:BatchDataset';
end
problem = MException(identifier,'Dataset %d: %s',index,cause.message);
throwAsCaller(addCause(problem,cause));
end

function commstr = historyCommand(implementation,executionMode,args)
% Keep full numeric precision for EEGLAB replay, including fractional epochs.
serialized = cell(1,numel(args)+4);
serialized(1:4) = {'''Implementation''',['''' implementation ''''], ...
    '''ExecutionMode''',['''' executionMode '''']};
for k=1:2:numel(args)
    serialized{k+4}=['''' char(args{k}) ''''];
    value=args{k+1};
    if ischar(value) || isstring(value)
        serialized{k+5}=['''' strrep(char(value),'''','''''') ''''];
    else
        serialized{k+5}=mat2str(value,17);
    end
end
commstr=['EEG = pop_cleanline_integrated(EEG, ' strjoin(serialized,', ') ');'];
end
