load('Data/tCs_all704.mat')
tC_save{106}=[];
file_init = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';


%% cycle through folders containing subjects
D = dir(file_init);         % D is a structure
D(1:3)=[];                  % Remove metadata rows


speed = cell(1,704);
for f = 1:length(D)       
    [sub_id,file_loc] = read_subject_file(f,D,file_init);
    
    if isfile([file_loc,sub_id])        % Check if sharing data is available for subject "currD"
        filename = [file_loc,sub_id];
        tab=sortrows(readtable(filename),4);
        if iscell(tab.X(1))             % Convert strings to doubles for X and Y
            tab=sortrows(readtable(filename),14);
            tab.X=strrep(tab.X,',','.');
            tab.Y=strrep(tab.Y,',','.');
            tab.X=str2double(tab.X);
            tab.Y=str2double(tab.Y);
        end

%         i=i+1;
        tC = tC_save{1,f};
        if ~isempty(tC)
            sw_speed = zeros(max(tC),1);
            for j = 1 : max(tC)
                pos=zeros(2,2); time=zeros(2,1);
                swipe = find(tC==j);
                pos(1,:)=[tab.X(swipe(1)),tab.Y(swipe(1))];
                pos(2,:)=[tab.X(swipe(end)),tab.Y(swipe(end))];
                time(1)=tab.Time(swipe(1));
                time(2)=tab.Time(swipe(end));
                sw_speed(j,1) = sqrt((pos(2,1)-pos(1,1))^2+(pos(2,2)-pos(1,2))^2)/(time(2)-time(1));
            end
            speed{f}=sw_speed;
        end
    end
end

save('Data/swipe_speeds.mat','speed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [sub_id,file_loc] = read_subject_file(f,D,file_init)
    currD = D(f).name;              % Get the current subdirectory name
    file_loc =[file_init,currD,'\'];
    sub_id = [currD,'.Sharing.TouchData.typed.csv'];
end