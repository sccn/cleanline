function [output,history,report,console] = drive_cleanline_test_dialog(invoke,inspect,accept)
% Exercise the real arg_guidialog controls and EEGLAB's modal return.
probe = timer('StartDelay',2,'ExecutionMode','fixedSpacing','Period',0.1, ...
    'TimerFcn',@(source,~) finishDialog(source,inspect,accept));
cleanup = onCleanup(@() delete(probe)); %#ok<NASGU>
start(probe);
[console,output,history,report] = evalc('invoke()');
if ~isempty(probe.UserData)
    rethrow(probe.UserData);
end
assert(strcmp(probe.Running,'off'),'The CleanLine dialog was not exercised.');
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
stop(probe);
try
    assert(isscalar(fig));
    inspect(fig);
    if accept
        % This is the OK callback's signal to EEGLAB's inputgui wait.
        set(ok,'UserData','retuninginputui');
    else
        close(fig);
    end
catch err
    probe.UserData = err;
    close(fig(1));
end
end
