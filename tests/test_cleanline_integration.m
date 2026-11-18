% Runnable synthetic parity and compatibility checks for the official plugin.
pluginRoot=fileparts(fileparts(mfilename('fullpath')));
eeglabRoot=getenv('CLEANLINE_EEGLAB_ROOT');
assert(isfolder(eeglabRoot),'Set CLEANLINE_EEGLAB_ROOT to the EEGLAB folder.');
addpath(eeglabRoot); eeglab('nogui');
addpath(genpath(pluginRoot),'-begin');
assert(strcmp(which('pop_cleanline'),fullfile(pluginRoot,'pop_cleanline.m')));
rng(20261118,'twister');
Fs=256; t=(0:Fs*30-1)/Fs;
EEG=eeg_emptyset;
EEG.data=single([sin(2*pi*8*t)+0.4*sin(2*pi*60*t); ...
    randn(1,numel(t))*0.1+sin(2*pi*60.2*t); ...
    sin(2*pi*13*t)+0.1*sin(2*pi*120*t); ...
    randn(1,numel(t))*0.1]);
EEG.nbchan=4; EEG.pnts=numel(t); EEG.trials=1; EEG.srate=Fs;
EEG.xmin=0; EEG.xmax=t(end);
EEG.chanlocs=struct('labels',{'1','2','3','4'});
args={'newversion',1,'linefreqs',60,'winsize',4,'winstep',1, ...
    'bandwidth',2,'taperbandwidth',2,'p',0.01,'pad',2, ...
    'tau',100,'sigtype','Channels','chanlist',1:4, ...
    'computepower',false,'verb',false};
reference=pop_cleanline(EEG,args{:});
[optimized,com,report]=pop_cleanline_integrated(EEG,args{:}, ...
    'Implementation','optimized','ExecutionMode','serial');
delta=max(abs(double(reference.data(:))-double(optimized.data(:))));
assert(delta==0,'Optimized serial differs from official V2.');
assert(isequal(EEG.event,optimized.event));
replayed=EEG; EEG=replayed; eval(com);
assert(isequal(EEG.data,optimized.data),'EEGLAB history replay failed.');
fprintf('Synthetic official V2 vs optimized serial: max abs error %.17g; history exact.\n',delta);
[original,~,~]=pop_cleanline_integrated(replayed,args{:},'Implementation','original');
assert(isequal(reference.data,original.data));
fprintf('Original V2 wrapper: exact.\n');
if license('test','Distrib_Computing_Toolbox')
    [parallel,~,~]=pop_cleanline_integrated(replayed,args{:}, ...
        'Implementation','optimized','ExecutionMode','parallel');
    deltaParallel=max(abs(double(reference.data(:))-double(parallel.data(:))));
    assert(deltaParallel==0,'Optimized parallel differs from official V2.');
    fprintf('Synthetic official V2 vs optimized parallel: max abs error %.17g.\n',deltaParallel);
end
assert(strcmp(which('pop_cleanline'),fullfile(pluginRoot,'pop_cleanline.m')));
fprintf('Legacy entry point remains available.\n');
legacyArgs={'newversion',0,'linefreqs',60,'winsize',4,'winstep',1, ...
    'bandwidth',2,'p',0.01,'pad',2,'sigtype','Channels', ...
    'chanlist',1,'computepower',false,'verb',false};
legacyReference=pop_cleanline(replayed,legacyArgs{:});
[legacyIntegrated,~,~]=pop_cleanline_integrated(replayed,legacyArgs{:}, ...
    'Implementation','legacy');
assert(isequal(legacyReference.data,legacyIntegrated.data));
fprintf('Historical V1 wrapper: exact.\n');
doubleInput=replayed; doubleInput.data=double(doubleInput.data);
selectedArgs=args;
selectedArgs{find(cellfun(@(x)ischar(x)&&strcmpi(x,'chanlist'),selectedArgs))+1}=[1 3];
doubleReference=doubleInput;
for ch=[1 3]
    oneArgs=selectedArgs;
    oneArgs{find(cellfun(@(x)ischar(x)&&strcmpi(x,'chanlist'),oneArgs))+1}=ch;
    one=pop_cleanline(doubleInput,oneArgs{:});
    doubleReference.data(ch,:)=one.data(ch,:);
end
doubleOptimized=pop_cleanline_integrated(doubleInput,selectedArgs{:}, ...
    'Implementation','optimized','ExecutionMode','serial');
assert(isequal(doubleReference.data,doubleOptimized.data));
assert(isequal(doubleInput.data(2,:),doubleOptimized.data(2,:)));
fprintf('Double input, nonconsecutive selected channels: max abs error 0.\n');
if license('test','Distrib_Computing_Toolbox')
    doubleParallel=pop_cleanline_integrated(doubleInput,selectedArgs{:}, ...
        'Implementation','optimized','ExecutionMode','parallel');
    assert(isequal(doubleReference.data,doubleParallel.data));
    fprintf('Nonconsecutive optimized parallel: max abs error 0.\n');
end
icaInput=replayed;
icaInput.icaact=double(icaInput.data);
icaInput.icawinv=eye(4);
icaInput.icaweights=eye(4);
icaInput.icasphere=eye(4);
icaInput.icachansind=1:4;
componentArgs=legacyArgs;
componentArgs{find(cellfun(@(x)ischar(x)&&strcmpi(x,'sigtype'),componentArgs))+1}='Components';
componentReference=pop_cleanline(icaInput,componentArgs{:});
componentIntegrated=pop_cleanline_integrated(icaInput,componentArgs{:}, ...
    'Implementation','legacy');
assert(isequal(componentReference.data,componentIntegrated.data));
fprintf('Historical V1 ICA component routing: exact.\n');
fig=figure('Visible','off');
menu=uimenu(fig,'Label','Tools','Tag','tools');
trystrs=struct('no_check','');
catchstrs=struct('new_and_hist','');
eegplugin_cleanline(fig,trystrs,catchstrs);
labels=get(get(menu,'Children'),'Label');
assert(numel(get(menu,'Children'))==1 && strcmp(labels,'CleanLine'));
close(fig);
fprintf('EEGLAB Tools menu: one CleanLine entry.\n');
