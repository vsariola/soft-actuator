ictObj = icdevice('niDMM.mdd', 'Dev1');
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
range = 1e6; % Ohm?
resolution = 6.5;

configuration = get(ictObj, 'configuration');
invoke(configuration, 'configuremeasurementdigits', measurementFunction, range, resolution);
% invoke(configuration, 'configurewaveformacquisition', measurementFunction, range, resolution, 100);

AutoTimeLimit = -1;

acquisition = get(ictObj, 'acquisition');
ohms =  invoke(acquisition, 'read', AutoTimeLimit);

fprintf('The current ohm is: %d Ohm\n', ohms);

%%

s = serial('COM3');
fopen(s);
line = fgetl(s);

n = 50;
minds = 0.1;
maxds = 0.6;
halfm = 10;
pausetime = 10;
duty_cycles = [linspace(minds,maxds,halfm) linspace(maxds,minds,halfm)];
m = length(duty_cycles);
data = zeros(m*n,8);
tic;
k = 1;
for i = 1:m
    % Set new duty cycle
    if (i > 1)
        dold = d;
        d = min(max(duty_cycles(i),0),1);   
        if (dold < d)
            ind = uint8(dold * 255):uint8(d * 255);
        else
            ind = uint8(dold * 255):-1:uint8(d * 255);
        end        
        for kk = ind
            fwrite(s,uint8(kk));
        end
    else
        d = min(max(duty_cycles(i),0),1);    
        fwrite(s,uint8(d * 255));
    end
    pause(pausetime);
    flushinput(s);
    line = fgetl(s);
    for j = 1:n
        line = fgetl(s);
        linevals = cell2mat(textscan(line,'%f'));
        ohms =  invoke(acquisition, 'read', AutoTimeLimit);
        data(k,:) = [d ohms linevals' toc];
        k = k +1;
    end
    system('canon_shooting.exe 0');
    system('canon_save.exe 1 d:\Kuvat\');
    fprintf('%.0f%% completed of the cycle\n',i*100/m);
    ohms =  invoke(acquisition, 'read', AutoTimeLimit);
    fprintf('The current ohm is: %d Ohm\n', ohms);
end

fwrite(s,uint8(0));

fclose(s);
save('d:\Kuvat\data.mat','data');
disconnect(ictObj);
delete(ictObj);
clear ictObj;