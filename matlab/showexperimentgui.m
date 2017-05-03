function showexperimentgui(paramdesc,expfunc)
    [~,params,~] = expfunc('defaults');               
    
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
    expfunc(args{:});        
end