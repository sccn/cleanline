function [EEGclean, commstr, report] = pop_cleanline_integrated(EEG, varargin)
%POP_CLEANLINE_INTEGRATED One-menu selector for historical V1 and both V2 engines.
% Original pop_cleanline and its command history remain supported unchanged.
% Usage: [EEG,com] = pop_cleanline_integrated(EEG, ...
%     'Implementation','optimized','ExecutionMode','serial', ...
%     'LineFrequencies',60,'ChanCompIndices',1:EEG.nbchan);
EEGclean = EEG;
commstr = '';
report = struct();
if isempty(varargin)
    [varargin, accepted] = cleanline_integrated_dialog(EEG);
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
    if strcmp(implementation,'legacy')
        args = [args {'Newversion',false}];
    else
        args = [args {'Newversion',true}];
    end
    [EEGclean,commstr] = pop_cleanline(EEG,args{:});
    report = struct('implementation',implementation, ...
        'executionMode','official','historyCommand',commstr);
    return
end
if strcmp(executionMode,'parallel')
    assert(license('test','Distrib_Computing_Toolbox'), ...
        'Parallel mode requires Parallel Computing Toolbox.');
end
% Reuse the official argument definitions and validation for the V2 options.
[~,~,~,~,~,~,g] = cleanline_check('EEG',EEG,args{:},'Newversion',true);
assert(EEG.trials==1,'Optimized V2 requires continuous EEG.');
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
clock = tic;
[EEGclean,effective] = cleanline_optimized.cleanLineNoise(EEG,lineNoise);
report = struct('implementation',implementation, ...
    'executionMode',executionMode,'effectiveLineFrequencies', ...
    effective.lineFrequencies,'elapsedCoreSeconds',toc(clock));
% Serialize caller options, including the selected mode, for EEGLAB replay.
serialized = cell(1,numel(args)+4);
serialized(1:4) = {'''Implementation''','''optimized''', ...
    '''ExecutionMode''',['''' executionMode '''']};
for k=1:2:numel(args)
    serialized{k+4}=['''' char(args{k}) ''''];
    value=args{k+1};
    if ischar(value) || isstring(value)
        serialized{k+5}=['''' char(value) ''''];
    else
        serialized{k+5}=mat2str(value);
    end
end
commstr=['EEG = pop_cleanline_integrated(EEG, ' strjoin(serialized,', ') ');'];
report.historyCommand=commstr;
end
