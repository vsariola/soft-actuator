function cycling_gui
    [~,params,~] = cycling_experiment('defaults');         
    
    paramdesc = struct(...
                'n','Number of pressure/resistance measurements per data point', ...
                'minvalue','Minimum duty cycle (0 - 1)', ...
                'maxvalue','Maximum duty cycle (0 - 1)', ...
                'm','m = Number of data points in one direction (total = 2*m+1)', ...
                'cycles','Number of cycles', ...
                'photofreq','Data points after photo is taken (0 = never take a photo, 1 = after every data point)', ...
                'resistancebits','Number of bits in resistance measurement', ...
                'resistancerange','Range of the resistance measurement, in Ohms', ...
                'ramp','Ramp duty cycle (0 = no, 1 = yes)', ...
                'ramppause','If ramping, what is the pause between one step of duty cycle (seconds)', ...
                'pause','Pause after adjusting duty cycle (seconds)', ...              
                'autosave','Autosave after this many data points',...
                'serialport','Serial port to use',...
                'devicename','NI DMM Device name',...
                'temperature','Record temperature (0 = no, 1 = yes)',...
                'temperaturedevicename','NI USB TC-01 Device name',...
                'dataroot','Directory save data');
    
    while(true)
        lines = {sprintf('Are these parameters ok?\n')};
        fnames = fieldnames(params)';
        lenf = length(fnames);            
        for i = 1:lenf
            if (isfield(paramdesc,fnames{i}))
                desc = paramdesc.(fnames{i});               
            else
                desc = fnames{i};                               
            end
            lines{end+1} = sprintf('%s:',desc);
            lines{end+1} = sprintf('     %s',num2str(params.(fnames{i})));
        end           
        
        choice = questdlg(lines,'Parameters OK?','Yes','No','Quit','Yes');
        if (strcmp(choice,'Yes'))
            break;
        elseif (strcmp(choice,'Quit'))
            return;
        end
        
        choice = questdlg('Do you want to use input the parameters manually or load parameters from an old experiment?','Use new or reuse old parameters?','Input','Reuse old','Quit','Input');
        if (strcmp(choice,'Input'))            
            dlg_title = 'Experimental parameters';
            num_lines = 1;               
            prompt = cell(1,lenf);
            defaultans = cell(1,lenf);
            for i = 1:lenf
                if (isfield(paramdesc,fnames{i}))
                    prompt{i} = strcat(paramdesc.(fnames{i}),':');
                else
                    prompt{i} = strcat(fnames{i},':');
                end
                defaultans{i} = num2str(params.(fnames{i}));
            end            
            half = round(length(prompt)/2);
            answer1 = inputdlg(prompt(1:half),dlg_title,num_lines,defaultans(1:half));        
            if (isempty(answer1))
                continue;
            end
            answer2 = inputdlg(prompt((half+1):end),dlg_title,num_lines,defaultans((half+1):end)); 
             if (isempty(answer2))
                continue;
            end
            answer1 = cellfun(@toNumericIfPossible,answer1,'UniformOutput',false);
            answer2 = cellfun(@toNumericIfPossible,answer2,'UniformOutput',false);
            for i = 1:half
                params.(fnames{i}) = toNumericIfPossible(answer1{i});
            end
            for i = (half+1):lenf
                params.(fnames{i}) = toNumericIfPossible(answer2{i-half});
            end
        elseif (strcmp(choice,'Reuse old'))
            oldexpfile = uigetfile('*.mat','Select old experiment where to load parameters');
            s = load(oldexpfile);
            params = s.params;
        else
            return;
        end        
    end
    
    args = [fieldnames(params)'; struct2cell(params)'];
    args = reshape(args,1,length(fieldnames(params))*2);
    cycling_experiment(args{:});
    
    %f1 = uigetfile('*.CR2','Select picture where the actuator is most straight');
    %f2 = uigetfile('*.CR2','Select picture where the actuator is most bent');    
end