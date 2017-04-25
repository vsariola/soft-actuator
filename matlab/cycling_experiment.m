function [data,params,datadir] = cycling_experiment(varargin)

mydir = fileparts(mfilename('fullpath'));

isnatural = @(x) x == floor(x);
isBoolean = @(x)islogical(x) || isnumeric(x) && all(x(:)==0 | x(:)==1);

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
addParameter(p,'n',5,isnatural);    
addParameter(p,'minvalue',0.2,@isnumeric);
addParameter(p,'maxvalue',0.7,@isnumeric);
addParameter(p,'m',30,isnatural);        
addParameter(p,'cycles',1,isnatural);        
addParameter(p,'photofreq',3);       
addParameter(p,'resistancebits',6.5);
addParameter(p,'resistancerange',1000000,@isnumeric);
addParameter(p,'ramp',true,isBoolean);
addParameter(p,'ramppause',0,@isnumeric);
addParameter(p,'pause',0,@isnumeric);
addParameter(p,'autosave',0,isnatural);
addParameter(p,'serialport','COM3');
addParameter(p,'devicename','Dev1');
addParameter(p,'temperature',false,isBoolean);
addParameter(p,'temperaturedevicename','Dev2');
addParameter(p,'dataroot',defaultdataroot);
parse(p,varargin{:});
params = p.Results;

datadir = sprintf('%scycling_%s\\',params.dataroot,datetime('now','Format','y_MM_d_HH_mm_ss'));

if (peekdefaults)
    data = [];
    return;
end

datafile = sprintf('%scycling_data.mat',datadir);

delete(instrfindall); % delete all old objects from memory


upramp = linspace(params.minvalue,params.maxvalue,params.m+1);
downramp = linspace(params.maxvalue,params.minvalue,params.m+1);
onecycle = [upramp(2:end) downramp(2:end)];
duty_cycles = [upramp(1) repmat(onecycle,1,params.cycles)];
maxk = length(duty_cycles);

if (params.temperature)
    data = zeros(maxk*params.n,10);
else
    data = zeros(maxk*params.n,9);
end

ictObj = icdevice('niDMM.mdd', params.devicename);
connect(ictObj);
disp(ictObj);

% http://zone.ni.com/reference/en-XX/help/370384N-01/dmmcref/canidmm_attr_function/
% NIDMM_VAL_DC_VOLTS	1	DC Voltage	All
% NIDMM_VAL_AC_VOLTS	2	AC Voltage with AC Coupling	All
% NIDMM_VAL_DC_CURRENT	3	DC Current	All
% NIDMM_VAL_AC_CURRENT	4	AC Current	All
% NIDMM_VAL_2_WIRE_RES	5	2-Wire Resistance	All
% NIDMM_VAL_4_WIRE_RES	101	4-Wire Resistance	NI 4060, NI 4065, and NI 4070/4071/4072
measurementFunction = 5;

configuration = get(ictObj, 'configuration');
invoke(configuration, 'configuremeasurementdigits', measurementFunction, params.resistancerange, params.resistancebits);
% invoke(configuration, 'configurewaveformacquisition', measurementFunction, range, resolution, 100);

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

s = serial(params.serialport);
fopen(s);
fgetl(s);

mkdir(datadir);

tic;
k = 1;
for i = 1:maxk
    % Set new duty cycle
    if (i > 1 && params.ramp)
        dold = d;
        d = min(max(duty_cycles(i),0),1);   
        if (dold < d)
            ind = uint8(dold * 255):uint8(d * 255);
        else
            ind = uint8(dold * 255):-1:uint8(d * 255);
        end        
        for kk = ind
            fwrite(s,uint8(kk));
            pause(params.ramppause);
        end
    else
        d = min(max(duty_cycles(i),0),1);    
        fwrite(s,uint8(d * 255));
    end
    pause(params.pause);
    flushinput(s);
    fgetl(s);
    for j = 1:params.n
        line = fgetl(s);
        linevals = cell2mat(textscan(line,'%f'));
        ohms =  invoke(acquisition, 'read', AutoTimeLimit);
        if (params.temperature)
            data(k,:) = [d ohms linevals' toc 0 ss.inputSingleScan()];
        else
            data(k,:) = [d ohms linevals' toc 0];
        end
        k = k + 1;
    end
    if (mod(i-1,params.photofreq) == 0)
        system(sprintf('"%s\\..\\canon\\canon_shooting.exe" 0',mydir));
        system(sprintf('"%s\\..\\canon\\canon_save.exe" 1 %s',mydir,datadir));
%         listing = dir(sprintf('%s\\*.CR2',datadir));
%         if (isempty(listing))
%             listing = dir(sprintf('%s\\*.JPG',datadir));
%         end        
%         [~,dx] = sort([listing.datenum],'descend');
%         newest = listing(dx(1)).name;
%         expression = 'IMG_(\d+)';
%         tokens = regexp(newest,expression,'tokens');
%         data(k,9) = str2double(tokens{1}{1});   
    end
    fprintf('%.0f%% completed of the cycle\n',i*100/maxk);    
    if (mod(i-1,params.autosave) == 0)
        save(datafile,'data','params');
    end
end

fwrite(s,uint8(0));
fclose(s);

save(datafile,'data','params');
disconnect(ictObj);
delete(ictObj);
clear ictObj;