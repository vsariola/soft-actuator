function [data,params,datadir] = demo_multichannel_experiment(varargin)

    isnatural = @(x) x == floor(x);
    isBoolean = @(x) islogical(x) || isnumeric(x) && all(x(:)==0 | x(:)==1);

    prog = [0 0 0 0;1 0 0 0;0 0 0 0;0 1 0 0;0 0 0 0;0 0 1 0;0 0 0 0;0 0 0 0;]';

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
    addParameter(p,'strainfreq',10000,@isnumeric);
    addParameter(p,'autosave',0,isnatural);
    addParameter(p,'serialport','COM3');
    addParameter(p,'devicename','cDAQ1Mod1');
    addParameter(p,'temperature',false,isBoolean);
    addParameter(p,'temperaturedevicename','Dev2');
    addParameter(p,'dataroot',defaultdataroot);
    addParameter(p,'n',5);
    addParameter(p,'steptime',30);
    parse(p,varargin{:});
    params = p.Results;

    datadir = sprintf('%sdemo_%s\\',params.dataroot,datetime('now','Format','y_MM_d_HH_mm_ss'));

    if (peekdefaults)
        data = [];
        return;
    end

    maxk = size(prog,2);
    datafile = sprintf('%sdemo_data.mat',datadir);

    delete(instrfindall); % delete all old objects from memory

    nis = daq.createSession('ni');
    nis.Rate = params.strainfreq;
    nis.NotifyWhenDataAvailableExceeds = 1000;
    nis.DurationInSeconds = params.steptime * maxk + 2;
    addChannel('ai0');
    addChannel('ai1');
    addChannel('ai2'); 
    
    function ch = addChannel(channel_name)
        ch = addAnalogInputChannel(nis,'cDAQ1Mod1',channel_name,'Bridge');
        ch.BridgeMode = 'Half';
        ch.NominalBridgeResistance = 120;        
    end
    
    straink = 1;
    straindata = zeros(nis.NumberOfScans,4);

    addlistener(nis,'DataAvailable',@saveStrain); 
 
    function saveStrain(~,event)
        newdata = [event.TimeStamps event.Data];
        remain = size(straindata,1)-straink+1;
        fits = min(size(newdata,1),remain);
        straindata(straink:(straink+fits-1),:) = newdata(1:fits,:);
        straink = straink + fits;
    end

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

    mkdir(datadir);

    data = [];
    flushinput(s);
    fgetl(s);
    tic;
    k = 1;
    i = 1;


    fgetl(s);
    startBackground(nis);
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
        c1 = uint8(prog(1,k) * 3);
        c2 = bitshift(uint8(prog(2,k) * 3),2);
        c3 = bitshift(uint8(prog(3,k) * 3),4);
        c4 = bitshift(uint8(prog(4,k) * 3),6);
        fwrite(s,c1+c2+c3+c4);
        flushinput(s);
        fgetl(s);
        line = fgetl(s);
        linevals = cell2mat(textscan(line,'%f'));
        for z = 1:params.n
            if (params.temperature)
                data(end+1,:) = [linevals' time k ss.inputSingleScan()];
            else
                data(end+1,:) = [linevals' time k];
            end        
        end
        if (mod(i,params.autosave) == 0)
            save(datafile,'data','params');
        end
        i = i + 1;
    end      

    fwrite(s,uint8(0));
    fclose(s);

    nis.wait();
    
    save(datafile,'data','straindata','params','prog');
    
end