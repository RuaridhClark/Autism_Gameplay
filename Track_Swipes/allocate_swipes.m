%% Create swipe files that identify each individual swipe
clear all
addpath("Data")

%% Trial
file_loc = NaN; % [data restricted]
load('tC_trial.mat','tC_save','name_save') % load each swipe identified 

[swipe_save] = create_swipe_files(tC_save,name_save,file_loc,'trial');

save('..\Create_adj\swipes_trial.mat','swipe_save','name_save')

%% Pretrial
file_loc = NaN; % [data restricted]
load('tC_pretrial.mat','tC_save','name_save') % load each swipe identified 

[swipe_save] = create_swipe_files(tC_save,name_save,file_loc,'pretrial');

save('..\Create_adj\swipes_pretrial.mat','swipe_save','name_save')

%%%%%%% Functions %%%%%%%
function [swipe_save] = create_swipe_files(tC_save,name_save,file_loc,option)
    m = 0; n = m;
    swipe_save = cell(1,length(tC_save));
    for i = 1:length(tC_save)
        tC=tC_save{i};
    
        if strcmp(option,'trial')
            [~,~,skip] = setup_trial(i,name_save,file_loc);
        elseif strcmp(option,'pretrial')
            [~,~,skip] = setup_pretrial(i,name_save,file_loc);
        end

        swipe = [];
        if ~skip
            n=n+1;
            swipe = cell(1,max(tC));
            for ii = 1 : max(tC)
                swipe{ii}=find(tC==ii);
            end
        end
        swipe_save{i}=swipe;
    end
end

function [filename,savename,skip] = setup_trial(i,name_save,swipe_data_loc)
    skip = 0;
    
    file_id = [name_save{i},'\',name_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([swipe_data_loc,file_id])
        filename = [swipe_data_loc,file_id];
        savename = ['subject_',name_save{i}];
    else
        filename=[];
        savename=[];
        skip = 1;
    end
end

function [filename,savename,skip] = setup_pretrial(i,name_save,swipe_data_loc)
    skip = 0;

    file_id = [name_save{i}];
    
    if isfile([swipe_data_loc,file_id])
        filename = [swipe_data_loc,file_id];
        savename = ['subject_',name_save{i}];
    else
        filename=[];
        savename=[];
        skip = 1;
    end
end