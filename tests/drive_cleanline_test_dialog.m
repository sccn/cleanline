function [output,history,report,console] = drive_cleanline_test_dialog(invoke,inspect,accept)
% Exercise the real arg_guidialog controls and EEGLAB's modal return.
probe = timer('StartDelay',2,'ExecutionMode','fixedSpacing','Period',0.1, ...
    'TimerFcn',@(source,~) finishDialog(source,inspect,accept));
cleanup = onCleanup(@() deleteProbe(probe)); %#ok<NASGU>
probe.UserData = struct('figure',[],'error',[]);
start(probe);
[console,output,history,report] = evalc('invoke()');
state = probe.UserData;
if ~isempty(state.error)
    rethrow(state.error);
end
assert(~isempty(state.figure),'The CleanLine dialog was not exercised.');
end

function finishDialog(probe,inspect,accept)
fig = findall(0,'Type','figure','Name','CleanLine Options');
if isempty(fig)
    return
end
ok = findobj(fig(1),'Tag','ok');
if isempty(ok) || isempty(findobj(fig(1),'Tag','cancel'))
    return
end
state = probe.UserData;
if isequal(state.figure,fig(1))
    return
end
try
    assert(isscalar(fig));
    assert(isempty(state.figure),'CleanLine opened the settings dialog more than once.');
    state.figure = fig;
    probe.UserData = state;
    inspect(fig);
    if accept
        % This is the OK callback's signal to EEGLAB's inputgui wait.
        set(ok,'UserData','retuninginputui');
    else
        close(fig);
    end
catch err
    state.error = err;
    probe.UserData = state;
    close(fig(1));
end
end

function deleteProbe(probe)
stop(probe);
delete(probe);
end
