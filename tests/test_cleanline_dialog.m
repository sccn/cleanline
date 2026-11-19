function test_cleanline_dialog
% Exercise the actual EEGLAB dialog and replay its history for each route.
pluginRoot = fileparts(fileparts(mfilename('fullpath')));
eeglabRoot = getenv('CLEANLINE_EEGLAB_ROOT');
assert(isfolder(eeglabRoot),'Set CLEANLINE_EEGLAB_ROOT to the EEGLAB folder.');
addpath(eeglabRoot);
eeglab('nogui');
addpath(genpath(pluginRoot),'-begin');
assert(strcmp(which('pop_cleanline_integrated'),fullfile(pluginRoot,'pop_cleanline_integrated.m')));

rng(20261119,'twister');
EEG = eeg_emptyset;
EEG.nbchan = 2;
EEG.srate = 250;
EEG.pnts = 2000;
EEG.trials = 1;
t = (0:EEG.pnts-1)/EEG.srate;
EEG.xmax = t(end);
EEG.data = single(0.1*randn(2,EEG.pnts) + [sin(2*pi*60*t);cos(2*pi*60*t)]);
EEG.chanlocs = struct('labels',{'A','B'});
invoke = @() pop_cleanline_integrated(EEG);
originalSpec = arg_report('rich',@cleanline_check,{'EEG',EEG});

[default,history,report] = drive_cleanline_test_dialog(invoke, ...
    @(fig) checkDefaults(fig,originalSpec,{'Channels'}),true);
assertRoute(EEG,default,history,report,'optimized','serial');
reference = pop_cleanline_integrated(EEG,'Implementation','optimized','ExecutionMode','serial');
assert(isequal(default.data,reference.data));
fprintf('GUI PASS: original settings and order, unchecked default Optimized V2.2 Serial, exact history.\n');

[legacy,history,report] = drive_cleanline_test_dialog(invoke,@selectLegacy,true);
assertRoute(EEG,legacy,history,report,'legacy','official');
expected = pop_cleanline(EEG,'Newversion',false,'LineFrequencies',60.1234567891);
assert(isequal(legacy.data,expected.data));
fprintf('GUI PASS: checked selects Historical V1 and exact replay preserves the entered frequency.\n');

[optimized,history,report,console] = drive_cleanline_test_dialog(invoke,@uncheckAndSelectV1Options,true);
assertRoute(EEG,optimized,history,report,'optimized','serial');
assert(isequal(optimized.data,default.data));
assert(contains(console,'spectrum normalization ignored') && contains(console,'plotting option ignored') ...
    && contains(console,'Scan for lines checkbox ignored'));
assert(contains(history,'''normSpectrum'', true') && contains(history,'''plotfigures'', true') ...
    && contains(history,'''scanforlines'', true'));
fprintf('GUI PASS: unchecking restores Optimized V2.2 and reports unsupported historical controls.\n');

assert(license('test','Distrib_Computing_Toolbox'),'Parallel Computing Toolbox is required.');
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('Processes',2);
    closePool = onCleanup(@() delete(pool)); %#ok<NASGU>
end
[parallelOutput,history,report] = drive_cleanline_test_dialog(invoke, ...
    @(fig) set(findobj(fig,'Tag','executionmode'),'Value',2),true);
assertRoute(EEG,parallelOutput,history,report,'optimized','parallel');
assert(isequal(parallelOutput.data,default.data) && isequal(gcp('nocreate'),pool));
fprintf('GUI PASS: Parallel selects Optimized V2.2 Parallel, reuses the pool and matches Serial exactly.\n');

[cancelled,history,report] = drive_cleanline_test_dialog(invoke, ...
    @(fig) set(findobj(fig,'Tag','useoldversion'),'Value',1),false);
assert(isequaln(cancelled,EEG) && isempty(history) && isempty(fieldnames(report)));
fprintf('GUI PASS: Cancel returns the input dataset and no history.\n');

ica = EEG;
ica.icaweights = eye(2);
ica.icasphere = eye(2);
ica.icawinv = eye(2);
ica.icaact = double(ica.data);
ica.icachansind = 1:2;
[channels,history,report] = drive_cleanline_test_dialog(@() pop_cleanline_integrated(ica), ...
    @(fig) checkDefaults(fig,originalSpec,{'Components','Channels'}),true);
assertRoute(ica,channels,history,report,'optimized','serial');
assert(isequal(channels.data,default.data));
fprintf('GUI PASS: ICA dataset offers Components and Channels, with Channels selected for Optimized V2.2.\n');

[components,history,report] = drive_cleanline_test_dialog( ...
    @() pop_cleanline_integrated(ica),@selectLegacyComponents,true);
assertRoute(ica,components,history,report,'legacy','official');
expected = pop_cleanline(ica,'Newversion',false,'SignalType','Components');
assert(isequal(components.data,expected.data));
fprintf('GUI PASS: explicit Historical V1 component cleaning matches the original API.\n');
end

function checkDefaults(fig,originalSpec,signalTypes)
keys = cellfun(@(aliases) aliases{1},{originalSpec.names},'UniformOutput',false);
assert(strcmp(keys{end},'newversion'));
keys = [keys(1:end-1) {'executionmode','useoldversion'}];
positions = zeros(1,numel(keys));
for k = 1:numel(keys)
    control = findobj(fig,'Tag',keys{k});
    assert(isscalar(control),'Missing control %s.',keys{k});
    bounds = get(control,'Position');
    positions(k) = bounds(2);
end
assert(all(diff(positions)<0),'The controls are out of order.');
labels = get(findobj(fig,'Style','text'),'String');
for k = 1:numel(originalSpec)-1
    assert(ismember(originalSpec(k).help{1},labels));
end
assert(isempty(findobj(fig,'Tag','newversion')));
assert(ismember('Use old version (CAUTION! There is a known bug, test use only)',labels));
assert(ismember('Optimized V2.2 execution mode',labels));
oldVersion = findobj(fig,'Tag','useoldversion');
assert(strcmp(get(oldVersion,'Style'),'checkbox') && get(oldVersion,'Value')==0);
execution = findobj(fig,'Tag','executionmode');
assert(isequal(cellstr(get(execution,'String')),{'Serial';'Parallel'}) && get(execution,'Value')==1);
signal = findobj(fig,'Tag','sigtype');
assert(strcmp(get(signal,'Style'),'popupmenu') && ...
    isequal(cellstr(get(signal,'String')),signalTypes(:)) && get(signal,'Value')==numel(signalTypes));
end

function selectLegacy(fig)
set(findobj(fig,'Tag','useoldversion'),'Value',1);
set(findobj(fig,'Tag','linefreqs'),'String','60.1234567891');
end

function uncheckAndSelectV1Options(fig)
checkbox = findobj(fig,'Tag','useoldversion');
set(checkbox,'Value',1);
set(checkbox,'Value',0);
set(findobj(fig,'Tag','normSpectrum'),'Value',1);
set(findobj(fig,'Tag','plotfigures'),'Value',1);
set(findobj(fig,'Tag','scanforlines'),'Value',1);
end

function selectLegacyComponents(fig)
set(findobj(fig,'Tag','useoldversion'),'Value',1);
set(findobj(fig,'Tag','sigtype'),'Value',1);
end

function assertRoute(EEG,output,history,report,implementation,executionMode)
assert(strcmp(report.implementation,implementation) && strcmp(report.executionMode,executionMode));
assert(strcmp(report.historyCommand,history) && contains(history,['''' implementation '''']));
eval(history);
assert(isequal(EEG.data,output.data),'GUI history replay differs.');
end
