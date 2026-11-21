function test_cleanline_parallel_pool
% Run in a fresh MATLAB session with CLEANLINE_EEGLAB_ROOT set.
pluginRoot = fileparts(fileparts(mfilename('fullpath')));
eeglabRoot = getenv('CLEANLINE_EEGLAB_ROOT');
assert(isfolder(eeglabRoot),'Set CLEANLINE_EEGLAB_ROOT to the EEGLAB folder.');
addpath(eeglabRoot);
eeglab('nogui');
addpath(genpath(pluginRoot),'-begin');
assert(strcmp(which('pop_cleanline_integrated'), ...
    fullfile(pluginRoot,'pop_cleanline_integrated.m')));
assert(license('test','Distrib_Computing_Toolbox'), ...
    'This test requires Parallel Computing Toolbox.');
assert(isempty(gcp('nocreate')),'Run this test in a session without an active pool.');

settings = parallel.Settings;
originalAutoCreate = settings.Pool.AutoCreate;
restoreSettings = onCleanup(@() restoreAutoCreate(settings,originalAutoCreate)); %#ok<NASGU>
settings.Pool.AutoCreate = false;
assert(~settings.Pool.AutoCreate);
closePool = onCleanup(@deleteTestPool); %#ok<NASGU>

EEG = eeg_emptyset;
EEG.srate = 256;
EEG.trials = 1;
EEG.nbchan = 2;
EEG.pnts = 8*EEG.srate;
t = (0:EEG.pnts-1)/EEG.srate;
EEG.xmax = t(end);
EEG.data = [sin(2*pi*8*t)+sin(2*pi*60*t); cos(2*pi*60*t)];
args = {'Implementation','optimized','SignalType','Channels', ...
    'ChanCompIndices',1:2,'ComputeSpectralPower',false,'VerboseOutput',false};

[parallelOutput,history] = pop_cleanline_integrated(EEG,args{:},'ExecutionMode','parallel');
pool = gcp('nocreate');
assert(~isempty(pool),'Optimized Parallel did not start a pool with auto-start disabled.');
fprintf('A PASS: Explicit startup with auto-start disabled, %d worker(s), pool still active.\n',pool.NumWorkers);

input = EEG;
eval(history);
assert(isequal(gcp('nocreate'),pool),'Optimized Parallel replaced the existing pool.');
assert(isequal(EEG.data,parallelOutput.data),'Optimized Parallel history replay differed.');
fprintf('B PASS: Optimized Parallel reused the same pool and replayed exactly.\n');

delete(pool);
assert(isempty(gcp('nocreate')));
serialOutput = pop_cleanline_integrated(input,args{:},'ExecutionMode','serial');
assert(isempty(gcp('nocreate')),'Optimized Serial started a pool.');
assert(isequal(serialOutput.data,parallelOutput.data));
fprintf('C PASS: Optimized Serial left the pool closed and matched Parallel exactly.\n');

mockFolder = tempname;
mkdir(mockFolder);
mockFile = fullfile(mockFolder,'parpool.m');
removeMock = onCleanup(@() removePoolMock(mockFolder,mockFile)); %#ok<NASGU>
writelines(["function varargout = parpool(varargin)"; ...
    "assert(isempty(varargin),'Expected default parpool arguments.');"; ...
    "error('CleanLineTest:PoolStartup','Simulated pool startup failure.');"; ...
    "end"],mockFile);
addpath(mockFolder,'-begin');
assert(strcmp(which('parpool'),mockFile));
try
    pop_cleanline_integrated(input,args{:},'ExecutionMode','parallel');
    error('CleanLineTest:UnexpectedSuccess','Optimized Parallel returned without a pool.');
catch err
    assert(strcmp(err.identifier,'CleanLine:ParallelPool') && ...
        contains(err.message,'Simulated pool startup failure.'), ...
        'Unexpected error: %s (%s)',err.message,err.identifier);
end
assert(isempty(gcp('nocreate')));
fprintf('E PASS: Startup failure raised CleanLine:ParallelPool with the underlying message.\n');
end

function restoreAutoCreate(settings,value)
settings.Pool.AutoCreate = value;
end

function deleteTestPool
pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end
end

function removePoolMock(folder,file)
if strcmp(which('parpool'),file)
    rmpath(folder);
end
if isfile(file)
    delete(file);
end
rmdir(folder);
end
