function [args, accepted] = cleanline_integrated_dialog(EEG)
%CLEANLINE_INTEGRATED_DIALOG Show V2.1 settings with the V2.2 engine choices.
% Only dataset fields consumed by cleanline_check need to pass through the
% argument dialog's initial and final validation.
dataset = struct('nbchan',EEG.nbchan,'trials',EEG.trials, ...
    'pnts',EEG.pnts,'srate',EEG.srate,'icawinv',EEG.icawinv);
g = arg_guidialog(@dialogSettings,'Parameters',{'EEG',dataset}, ...
    'Title','CleanLine Options','Invoke',false);

args = {};
accepted = ~isempty(g);
if ~accepted
    return
end

implementation = 'optimized';
if g.useoldversion
    implementation = 'legacy';
else
    % Keep the historical V2 notices for settings exposed by this dialog.
    if g.normSpectrum
        disp('Warning: spectrum normalization ignored; use the old version for this option');
    end
    if g.plotfigures
        disp('Warning: plotting option ignored; use the old version for this option');
    end
    if g.scanforlines
        disp('Warning: Scan for lines checkbox ignored; Optimized V2.2 uses the scan bandwidth');
    end
end

if EEG.trials > 1
    % arg_guidialog shows numeric defaults with 15 significant digits.
    % Restore the full epoch duration when either field retains that value.
    epochSeconds = EEG.pnts/EEG.srate;
    displayedSeconds = str2double(hlp_tostring(epochSeconds));
    if isequal(g.winsize,displayedSeconds)
        g.winsize = epochSeconds;
    end
    if isequal(g.winstep,displayedSeconds)
        g.winstep = epochSeconds;
    end
end

args = [{'Implementation',implementation,'ExecutionMode',lower(g.executionmode)}, ...
    hlp_struct2varargin(g,'suppress',{'EEG','arg_direct','executionmode','useoldversion'})];
end

function g = dialogSettings(varargin)
dataset = arg_extract(varargin,'EEG',[],[]);
g = arg_define(0,varargin,arg_nogui('EEG',mandatory),{@dialogSpec,{dataset}});
end

function spec = dialogSpec(reportType,dataset)
% Start from the stock V2.1 definitions. Channels is a compatible default
% for Optimized V2.2, with the original ICA-dependent choices left visible.
spec = arg_report(reportType,@cleanline_check,{'EEG',dataset,'SignalType','Channels'});
names = cellfun(@(aliases) aliases{1},{spec.names},'UniformOutput',false);
% The bundled BCILAB subset can lack PropertyType's type inference.
[spec(cellfun(@islogical,{spec.value})).type] = deal('logical');
spec(strcmp(names,'sigtype')).type = 'char';
power = strcmp(names,'computepower');
spec(power).help{2} = [spec(power).help{2} '. This option applies to the old version.'];
spec(strcmp(names,'verb')).help{2} = 'This option controls the old version''s progress output.';
spec(strcmp(names,'newversion')) = [];

mode = arg('executionmode','Serial',{'Serial','Parallel'}, ...
    {'Optimized V2.2 execution mode','Applies when Use old version is unchecked. Serial is the default.'}, ...
    'type','char');
oldVersion = arg('useoldversion',false,[], ...
    {'Use old version','Unchecked runs Optimized V2.2 (default). Checked runs Historical V1.'}, ...
    'type','logical');
spec = [spec mode{1}(reportType,mode{2}{:}) oldVersion{1}(reportType,oldVersion{2}{:})];
% The mandatory checkbox label is longer than arg's summary length limit.
spec(end).help{1} = 'Use old version (CAUTION! There is a known bug, test use only)';
end
