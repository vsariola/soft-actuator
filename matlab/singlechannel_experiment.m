function [data,params,datadir] = demo_experiment(varargin)

isnatural = @(x) x == floor(x);
isBoolean = @(x) islogical(x) || isnumeric(x) && all(x(:)==0 | x(:)==1);

peekdefaults = ~isempty(varargin) && strcmp(varargin{1},'defaults');
if (peekdefaults)    
    varargin = {};
end

if ispc
    defaultdataroot = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    defaultdataroot = getenv('HOME');
end       

p = inputParser;   
addParameter(p,'resistancebits',6.5);
addParameter(p,'resistancerange',100,@isnumeric);
addParameter(p,'autosave',0,isnatural);
addParameter(p,'serialport','COM3');
addParameter(p,'devicename','Dev1');
addParameter(p,'temperature',false,isBoolean);
addParameter(p,'temperaturedevicename','Dev2');
addParameter(p,'dataroot',defaultdataroot);
addParameter(p,'n',5);
addParameter(p,'steptime',15);
parse(p,varargin{:});
params = p.Results;

datadir = sprintf('%ssinglechannel_%s\\',params.dataroot,datetime('now','Format','y_MM_d_HH_mm_ss'));

if (peekdefaults)
    data = [];
    return;
end

datafile = sprintf('%ssinglechannel_data.mat',datadir);

delete(instrfindall); % delete all old objects from memory

ictObj = icdevice('niDMM.mdd', params.devicename);
connect(ictObj);
disp(ictObj);

measurementFunction = 5;

configuration = get(ictObj, 'configuration');
invoke(configuration, 'configuremeasurementdigits', measurementFunction, params.resistancerange, params.resistancebits);
AutoTimeLimit = -1;

acquisition = get(ictObj, 'acquisition');
ohms =  invoke(acquisition, 'read', AutoTimeLimit);

fprintf('The current resistance is: %d ohm\n', ohms);

if (params.temperature)
    ss = daq.createSession('ni');
    addAnalogInputChannel(ss,params.temperaturedevicename,0, 'Thermocouple');   
    tc = ss.Channels(1);
    set(tc);
    tc.ThermocoupleType = 'J';    
    tc.Units = 'Celsius';
end

%%

prog = [0 1 0 1 0 0 0] * 0.5;

s = serial(params.serialport);
fopen(s);

mkdir(datadir);

data = [];
flushinput(s);
fgetl(s);
tic;
k = 1;
i = 1;
maxk = size(prog,2);

fgetl(s);
while(true)
    time = toc;    
    knew = floor(time / params.steptime)+1;    
    if (knew > maxk)
        break;
    end
    if (knew > k)
        fprintf('%.0f%% completed of the cycle\n',(knew-1)/maxk*100);    
    end
    k = knew;    
    c1 = uint8(prog(1,k) * 255);    
    fwrite(s,c1);
    flushinput(s);
    fgetl(s);
    line = fgetl(s);
    linevals = cell2mat(textscan(line,'%f'));
    for z = 1:params.n
        ohms =  invoke(acquisition, 'read', AutoTimeLimit);
        if (params.temperature)
            data(end+1,:) = [ohms linevals' time k ss.inputSingleScan()];
        else
            data(end+1,:) = [ohms linevals' time k];
        end        
    end
    if (mod(i,params.autosave) == 0)
        save(datafile,'data','params');
    end
    i = i + 1;
end      

fwrite(s,uint8(0));
fclose(s);

save(datafile,'data','params');
disconnect(ictObj);
delete(ictObj);
clear ictObj;