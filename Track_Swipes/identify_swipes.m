%% Identify swipes, by associating swipes with touch phase starting 0 and ending 2.
% Use spline curves and linear interpolation to find the continuation of swipes, 
% when multiple swipes occur simultaneously. Also check for and correct errors in the touch count.

Data_loc = NaN; % [data folder shared with reviewers]

%% Trial
% folder_loc = [Data_loc,'\Data\PlayCare\']; % [data restricted]
% [tC_save,name_save] = track_swipes(folder_loc,'trial');
% save("Data\tC_trial_dummy.mat","name_save","tC_save")

%% Pretrial
folder_loc = [Data_loc,'\Data\Krysiek_data\subject_data\']; % [data restricted]
[tC_save,name_save] = track_swipes(folder_loc,'pretrial');
save("Data\tC_pretrial_dummy.mat","name_save","tC_save")

%%%%%%% Functions %%%%%%%
function [tC_save,name_save] = track_swipes(folder_loc,option)
    
    %% cycle through folders containing subjects
    D = dir(folder_loc);         % D is a structure
    D(1:3)=[];                  % Remove metadata rows

    %% initialise variables
    name_save=cell(1,length(D));
    tC_save=cell(1,length(D));

    for f = 1:length(D)       
        currD = D(f).name; % Get the current subdirectory name
        [sub_id,file_loc] = read_subject_file(f,D,folder_loc,option);
    
        if isfile([file_loc,sub_id])        % Check if sharing data is available for subject "currD"
            filename = [file_loc,sub_id];
            [tab, starts, ends] = read_swipe_file(filename,option);    % tab is list of recorded touches, start and end are lists of start and end swipe touches
    
            [valid,~,starts,~] = check_validity(tab,starts,ends);  % Check and remove duplicate touch recordings
    
            tC = zeros(length(tab.X),1);
            swipe_ID = 1:1:length(starts);  % set swipe IDs
            tC(starts)=swipe_ID;            % tC maps touch to swipe ID
            pos=cell(1,max(swipe_ID));      % initialise
    
            %% Either assign touch to swipe ID or update pos with position info
            for j = 1 : length(valid)
                vj = valid(j);
                if tC(vj)~=0    % if tC is assigned
                    pos{tC(vj)} = [pos{tC(vj)};tab.X(vj),tab.Y(vj)];  % append X,Y position to pos
                else            % if tC is unassigned to swipe ID, then find next touch point and assign swipe ID
                    [tC,starts,pos,swipe_ID] = find_continuation(j,tab,tC,starts,pos,swipe_ID,valid);
                end
            end
            tC_save{f}=tC;        % save list of swipe_IDs for all touch points
            name_save{f}=erase(currD,'.mat');    % save subject ID
        end
    end
end

function [sub_id,file_loc] = read_subject_file(f,D,folder_loc,option)
    currD = D(f).name;              % Get the current subdirectory name
    file_loc =[folder_loc,currD,'\'];
    if strcmp(option,'trial')
        sub_id = [currD,'.Sharing.TouchData.typed.csv'];
    elseif strcmp(option,'pretrial')
        sub_id = '';
    end
end

function [tab, starts, ends] = read_swipe_file(filename,option)

    if strcmp(option,'trial')
        tab=sortrows(readtable(filename),4);
        if iscell(tab.X(1))             % Convert strings to doubles for X and Y
            tab=sortrows(readtable(filename),14);
            tab.X=strrep(tab.X,',','.');
            tab.Y=strrep(tab.Y,',','.');
            tab.X=str2double(tab.X);
            tab.Y=str2double(tab.Y);
        end
    elseif strcmp(option,'pretrial')
        load(filename,'subject_table');
        tab=sortrows(subject_table,8);
        tab = renamevars(tab,["x","y","time","touch_phase"],["X","Y","Time","TouchPhase"]);
    end

    % Reduce tab to unique values
    [~,at,~] = unique(tab.Time,'stable');
    [~,ax,~] = unique(tab.X,'stable');
    [~,ay,~] = unique(tab.Y,'stable');

    keep = unique([at;ax;ay],'first');

    tab=tab(keep,:);

    starts = find(tab.TouchPhase==0);   % swipe starts  
    ends = find(tab.TouchPhase==3);     % swipe ends

    if starts(1)~=1     % set first touch as start (even if it doesn't have touch phase 0)
        starts = [1;starts];
    end
end

function [valid,invalid,starts,ends] = check_validity(tab,starts,ends)
    %% Check and remove duplicated touch recordings

    invalid=[];
    for j = 2 : length(tab.X)           % for each touch recording
        if (tab.X(j-1)-tab.X(j))== 0 && (tab.Y(j-1)-tab.Y(j))==0 ...
            && abs(tab.TouchPhase(j)-tab.TouchPhase(j-1))~=3 % Check for repeats and only allow them if a touch is ending
            invalid=[invalid,j];        % add touch to invalid list
            if ismember(j,starts)
                starts(starts==j)=[];   % remove from list of starts
            end
            if ismember(j,ends)
                ends(ends==j)=[];       % remove from list of ends
            end
        end
    end
    
    valid=1:length(tab.X);
    valid(ismember(valid,invalid))=[];  % remove all invalid touches from valid list
end

function [tC,starts,pos,swipe_ID] = find_continuation(j,tab,tC,starts,pos,swipe_ID,valid)
    %% initialise
    nm_srch = 750;                  % number of previous touches checked (restricted for processing speed)
    tC_check=[]; k=j; 
    dist_all=ones(nm_srch,1)*10000; % set large initial distances 
                
    %% check up to nm_srch prior touches for best match
    for jj = 1:nm_srch              
        k=k-1;                      % check previous touches for swipe continuation
        if k>0 && ~ismember(tC(valid(k)),tC_check) % if touch is available and not yet checked
            vj = valid(k);              % cycle through only valid touches
            if size(pos{tC(vj)},1)>1    % if more than one touch position in swipe
                xx=tab.X(vj);
                [xx,yy] = predict_next_pos(xx,pos,tC,vj);                   % predict next pos using splines and linear interp.
            else
                xx = tab.X(vj); yy = tab.Y(vj);
            end
            dist_all(jj) = check_dist(valid(j),xx,yy,tab);	% Predicted vs actual Euclidean distance
            tC_check =[tC_check,tC(vj)];                    % list of checked touches
        end
    end
    [I,dist_all] = find_closest_match(j,dist_all,tC,valid,swipe_ID,tab);
    
    if min(dist_all) == 10000   % if no connection found in previous nm_srch touches, create new start 
        tC(valid(j))=max(tC)+1;
        starts=[starts;valid(j)];
        pos{tC(valid(j))}=[tab.X(valid(j)),tab.Y(valid(j))];
        swipe_ID=[swipe_ID,max(swipe_ID)+1];
        disp([num2str(valid(j)),' is now a start'])
    else                        % connection found
        tC(valid(j))=tC(valid(j-I));
        pos{tC(valid(j))}=[pos{tC(valid(j))};tab.X(valid(j)),tab.Y(valid(j))];
    end

    if tab.TouchPhase(valid(j))==3
        swipe_ID(find(swipe_ID==tC(valid(j))))=[]; % remove swipe as an option if end touch phase is reached
    end
end

function [xx,yy] = predict_next_pos(xx,pos,tC,vj)
    if length(pos{tC(vj)}(:,1))>4           % swipe greater than 4 points
        x_in = pos{tC(vj)}(end-4:end,1);    % use only last 4 touch points
        y_in = pos{tC(vj)}(end-4:end,2);
    else
        x_in = pos{tC(vj)}(:,1);
        y_in = pos{tC(vj)}(:,2);
    end

    try         % use spline to predict next swipe touch position for known x position
        yy = spline(x_in,y_in,xx);
    catch       % if spline cannot be computed use linear interpolation with last two points
        xx=2*pos{tC(vj)}(end,1)-pos{tC(vj)}(end-1,1);
        yy=2*pos{tC(vj)}(end,2)-pos{tC(vj)}(end-1,2);
    end
end

function [dist] = check_dist(cmpr,xx,yy,tab) 
    dist=sqrt((tab.X(cmpr)-xx)^2+(tab.Y(cmpr)-yy)^2);
end

function [I,dist_all] = find_closest_match(j,dist_all,tC,valid,swipe_ID,tab)
    [~,I]=min(dist_all);                        % find closest match (predicted versus actual)
    while ~ismember(tC(valid(j-I)),swipe_ID) ...
            || tab.TouchPhase(valid(j-I))==3	% check if touch has already been assigned or swipe has already ended
        dist_all(I)=10000;
        [~,I]=min(dist_all);
        if min(dist_all) == 10000
            break
        end
    end
end