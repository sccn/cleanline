function cleanline_optimized_progress(~, window, count)
% Update the optional channel progress window after a parallel worker completes.
if ishghandle(window)
    completed=getappdata(window,'cleanlineCompleted');
    if isempty(completed), completed=0; end
    completed = completed + 1;
    setappdata(window,'cleanlineCompleted',completed);
    waitbar(min(completed/count,1),window);
end
end
