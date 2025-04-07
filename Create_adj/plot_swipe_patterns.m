%% Create adjacency matrices
% Track the zones that swipes pass through zones to create an adjacency matrix 
% linking the origin and destination zones of swipes during gameplay.

clear all

zone = defineZonePositions(); % define game zones

%%%% Trial data
swipe_data_loc = NaN; % [data restricted]
swipe_data_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';
load('swipes_trial.mat');

adj_processing(name_save,swipe_save,swipe_data_loc,zone,'trial');
%%%%

%%%% Pre-trial data
swipe_data_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\Krysiek_data\subject_data\'; %NaN; % [data restricted]
load('swipes_pretrial.mat')

adj_processing(name_save,swipe_save,swipe_data_loc,zone,'pretrial');
%%%%

%%%%%% Functions %%%%%%
function [] = adj_processing(name_save,swipe_save,swipe_data_loc,zone,option)
    load('subject_details_combine.mat')
    

    % Process each swipe
    for i = 1:length(name_save)
        if ~isempty(name_save{i})
            Ind = find(strcmp(subject_details_combine(:,1),name_save{i}));
            titlename = subject_details_combine{Ind,5};
            if strcmp(titlename,'TD')
                titlename = 'WP';
            end
    
            swipe = swipe_save{i};
        
            if strcmp(option,'trial')
                [filename,savename,skip] = setup_trial(i,name_save,swipe_data_loc);
            elseif strcmp(option,'pretrial')
                [filename,savename,skip] = setup_pretrial(i,name_save,swipe_data_loc);
                if ~strcmp(titlename,'WP')
                    skip = 1;
                end
            end
        
            if ~skip
                [adj] = adj_create(option,zone,filename,swipe,savename,titlename);    
    %ADD BACK IN             save(['..\adjs\adj_zones\',savename],'adj');  % Save adj
            end
        end
    end
end

function zone = defineZonePositions()
    % Initialize yl
    yl = 768;
    
    % Initialize zone
    zone = zeros(16,4);
    zone(1,:)=[19, 227, yl-752, yl-510];
    zone(2,:)=[234, 774, yl-678, yl-510];
    zone(3,:)=[781, 1004, yl-752, yl-510];
    zone(4,:)=[106, 101+126, yl-400-96, yl-350]; 
    zone(5,:)=[330,325+126, yl-400-96, yl-350];  
    zone(6,:)=[556, 551+126, yl-400-96, yl-350]; 
    zone(7,:)=[781, 776+126, yl-400-96, yl-350];  
    zone(8,:)=[19, 252, yl-190, yl-27];
    zone(9,:)=[257, 420, yl-190, yl-27];
    zone(10,:)=[426, 589, yl-190, yl-27];
    zone(11,:)=[595, 758, yl-190, yl-27];
    zone(12,:)=[770, 1005, yl-190, yl-27];
    zone(13,:)=[19, 309, yl-508, yl-196];
    zone(14,:)=[311, 531, yl-508, yl-196];
    zone(15,:)=[533,753, yl-508, yl-196];
    zone(16,:)=[755, 1004, yl-508, yl-196];
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

    file_id = [name_save{i},'.mat'];
    
    if isfile([swipe_data_loc,file_id])
        filename = [swipe_data_loc,file_id];
        savename = ['subject_',name_save{i}];
    else
        filename=[];
        savename=[];
        skip = 1;
    end
end

function [adj] = adj_create(option,zone,filename,swipe,savename,titlename)

    if strcmp(option,'trial')
        tab = sortrows(readtable(filename),4); % sort table according to time
        if iscell(tab.X(1)) % Convert strings to doubles for X and Y
            tab = sortrows(readtable(filename),14); % sort table according to time
            tab.X = str2double(strrep(tab.X,',','.'));
            tab.Y = str2double(strrep(tab.Y,',','.'));
        end
    elseif strcmp(option,'pretrial')
        load(filename,'subject_table');
        tab=sortrows(subject_table,8); % sort table according to time
        tab = renamevars(tab,["x","y","time"],["X","Y","Time"]);
    end

    % Reduce tab to unique values
    [~,at,~] = unique(tab.Time,'stable');
    [~,ax,~] = unique(tab.X,'stable');
    [~,ay,~] = unique(tab.Y,'stable');

    keep = unique([at;ax;ay],'first');

    tab=tab(keep,:);

    adj = zeros(16,16);
%     con_list = [];      
    figure
    for ii = 1 : 16
        hold on
        plot([zone(ii,1),zone(ii,2)],[zone(ii,3),zone(ii,3)],'k')
        plot([zone(ii,1),zone(ii,2)],[zone(ii,4),zone(ii,4)],'k')
        plot([zone(ii,1),zone(ii,1)],[zone(ii,3),zone(ii,4)],'k')
        plot([zone(ii,2),zone(ii,2)],[zone(ii,3),zone(ii,4)],'k')
    end
    % Process each swipe
    for m = 1 : length(swipe)
        k = swipe{m};
        hold on
        p = plot(tab.X(k),tab.Y(k),'c',LineWidth=2);
        scatter(tab.X(k),tab.Y(k),'m.');
    end
    for m = 1 : length(swipe)
        k = swipe{m};
        s1 = scatter(tab.X(k(1)),tab.Y(k(1)),'k.');
        s2 = scatter(tab.X(k(end)),tab.Y(k(end)),'ko');
    end
    axis off
    axis equal
    title(titlename)
    saveas(gcf,[savename,'.png'])
    close gcf
end