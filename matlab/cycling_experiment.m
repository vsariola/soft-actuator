function [data,params,datadir] = cycling_experiment(varargin)

mydir = fileparts(mfilename('fullpath'));

isnatural = @(x) x == floor(x);

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
addParameter(p,'minvalue',0.2,@isnumeric);
addParameter(p,'maxvalue',0.6,@isnumeric);
addParameter(p,'numsteps',30,isnatural); 
addParameter(p,'pressuremeasurementsperstep',30,isnatural);  
addParameter(p,'cycles',1,isnatural);        
addParameter(p,'serialport','COM4');
addParameter(p,'dataroot',defaultdataroot);
parse(p,varargin{:});
params = p.Results;

datadir = sprintf('%scycling_%s\\',params.dataroot,datetime('now','Format','y_MM_d_HH_mm_ss'));

if (peekdefaults)
    data = [];
    return;
end

analogTokPa = @(x) (x/1024.0 - 0.1)*100.0/0.8 * 6.89475729;

datafile = sprintf('%scycling_data.mat',datadir);

upramp = linspace(params.minvalue,params.maxvalue,params.numsteps+1);
downramp = linspace(params.maxvalue,params.minvalue,params.numsteps+1);
onecycle = [upramp(2:end) downramp(2:end)];
duty_cycles = [upramp(1) repmat(onecycle,1,params.cycles)];
maxk = length(duty_cycles);

data = zeros(maxk,3);

%%

s = serial(params.serialport);
fopen(s);
fgetl(s);

mkdir(datadir);

    flushinput(s);
    fgetl(s);

tic;
k = 1;
for i = 1:maxk
    % Set new duty cycle
    d = min(max(duty_cycles(i),0),1);    
    fwrite(s,uint8(d * 255));
    for j = 1:params.pressuremeasurementsperstep
        line = fgetl(s);
        linevals = cell2mat(textscan(line,'%f'));        
        data(k,:) = [d analogTokPa(linevals(1)) toc];
        k = k + 1;
    end
    fprintf('%.0f%% completed of the cycle\n',i*100/maxk);    
end

fwrite(s,uint8(0));
fclose(s);

save(datafile,'data','params');