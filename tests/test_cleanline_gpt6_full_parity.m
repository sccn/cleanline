% Adapted from CleanLine-GPT6-Optimized/tests/test_cleanline_gpt6_full_parity.m.
% Compare full-duration output directly against this checkout's official V2.
pluginRoot=fileparts(fileparts(mfilename('fullpath')));
eeglabRoot=getenv('CLEANLINE_EEGLAB_ROOT');
setFile=getenv('CLEANLINE_TEST_SET');
assert(isfolder(eeglabRoot) && isfile(setFile), ...
    'Set CLEANLINE_EEGLAB_ROOT and CLEANLINE_TEST_SET.');
addpath(eeglabRoot); eeglab('nogui');
addpath(genpath(pluginRoot),'-begin');
assert(strcmp(which('pop_cleanline'),fullfile(pluginRoot,'pop_cleanline.m')));
[setDir,setName,setExt]=fileparts(setFile);
EEG=pop_loadset('filename',[setName setExt],'filepath',setDir);
args={'newversion',1,'linefreqs',60,'winsize',4,'winstep',1, ...
    'bandwidth',2,'taperbandwidth',2,'p',0.01,'pad',2, ...
    'tau',100,'sigtype','Channels','chanlist',1:33, ...
    'computepower',false,'verb',false};
for minutes=[5 30]
    input=EEG;
    input.pnts=minutes*60*EEG.srate;
    assert(size(EEG.data,2)>=input.pnts && size(EEG.data,1)>=33);
    input.data=EEG.data(:,1:input.pnts);
    input.xmax=input.xmin+(input.pnts-1)/input.srate;
    input.event=EEG.event([EEG.event.latency]<=input.pnts);
    reference=pop_cleanline(input,args{:});
    [serial,~,serialReport]=pop_cleanline_integrated(input,args{:}, ...
        'Implementation','optimized','ExecutionMode','serial');
    assert(isequal(reference.data,serial.data));
    assert(isequal(input.event,serial.event));
    assert(isequal(serialReport.effectiveLineFrequencies,[60 120]));
    fprintf('%d min, 33 channels: official V2 vs optimized serial max abs error 0.\n',minutes);
    if license('test','Distrib_Computing_Toolbox')
        [parallel,~,parallelReport]=pop_cleanline_integrated(input,args{:}, ...
            'Implementation','optimized','ExecutionMode','parallel');
        assert(isequal(reference.data,parallel.data));
        assert(isequal(input.event,parallel.event));
        assert(isequal(parallelReport.effectiveLineFrequencies,[60 120]));
        fprintf('%d min, 33 channels: official V2 vs optimized parallel max abs error 0.\n',minutes);
    end
end
