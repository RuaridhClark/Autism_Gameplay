%% Create adjacency matrices
% Track the zones that swipes pass through zones to create an adjacency matrix 
% linking the origin and destination zones of swipes during gameplay.

clear all

zone = defineZonePositions(); % define game zones

%%%% Trial data
swipe_data_loc = NaN; % [data restricted]
load('swipes_trial.mat');
save_adj_loc = '..\adjs\';

adj_processing(name_save,swipe_save,swipe_data_loc,zone,'trial');

%%%%

%%%% Pre-trial data
swipe_data_loc = NaN; % [data restricted]
load('swipes_pretrial.mat')

adj_processing(name_save,swipe_save,swipe_data_loc,zone,'pretrial');

%%%%

%%%%%% Functions %%%%%%
function [] = adj_processing(name_save,swipe_save,swipe_data_loc,zone,option)

    % Process each swipe
    for i = 1:length(name_save)
        swipe = swipe_save{i};
    
        if strcmp(option,'trial')
            [filename,savename,skip] = setup_trial(i,name_save,swipe_data_loc);
        elseif strcmp(option,'pretrial')
            [filename,savename,skip] = setup_pretrial(i,name_save,swipe_data_loc);
        end
    
        if ~skip
            [adj] = adj_create(option,zone,filename,swipe);    
            save([save_adj_loc,'adj_zones\',savename],'adj');  % Save adj
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

function [adj] = adj_create(option,zone,filename,swipe)

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
        tab = renamevars(tab,["x","y"],["X","Y"]);
    end

    adj = zeros(16,16);
          
    % Process each swipe
    for m = 1 : length(swipe)
        prev_n = [];
        n = 1; % start from swipe m first contact point
        while isempty(prev_n) && n <= length(swipe{m})
            k = swipe{m}(n);
            nn = [];
            for jj = 1 : 16
                ii = 17 - jj;
                if zone(ii,1) < tab.X(k) && tab.X(k) < zone(ii,2) && zone(ii,3) < tab.Y(k) && tab.Y(k) < zone(ii,4)
                    nn = ii; % current contact point is within the bounds of the zone ii
                end
            end
            
            if ~isempty(nn)
                prev_n = nn; % assign prev_n for first zone swipe enters
            end
            n = n + 1; % next swipe contact point
        end
     
        nn = [];
        n = -1; % Start from the end of swipe
        while isempty(nn) && ~isempty(prev_n)
            n = n + 1;
            k = swipe{m}(end - n);
            
            for jj = 1 : 16
                ii = 17 - jj;
                if zone(ii,1) < tab.X(k) && tab.X(k) < zone(ii,2) && zone(ii,3) < tab.Y(k) && tab.Y(k) < zone(ii,4)
                    nn = ii; % current contact point is within the bounds of the zone ii
                end
            end
            
            if ~isempty(prev_n) && ~isempty(nn)
                adj(prev_n,nn) = adj(prev_n,nn) + 1; % Increment the adjacency matrix at (prev_n, nn)
                prev_n = nn; % Update prev_n to nn
            end
        end
    end
end