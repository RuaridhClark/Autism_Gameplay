clear all
folder1 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Krysiek_data\subject_data';
addpath(folder1)
% file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';
% file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Krysiek_data\subject_data';
load('swipes_Krysiek_ond.mat')

load('subject_details_Krysiek_ond.mat','subject_details')

% yl = 768;
% obj(1,1:4)=[19, 227, yl-752, yl-525];
% obj(2,1:4)=[234, 774, yl-678, yl-525];
% obj(3,1:4)=[781, 1004, yl-752, yl-525];
% obj(4,1:4)=[106, 101+126, yl-400-86, yl-400]; 
% obj(5,1:4)=[330,325+126, yl-400-86, yl-400];  
% obj(6,1:4)=[556, 551+126, yl-400-86, yl-400]; 
% obj(7,1:4)=[781, 776+126, yl-400-86, yl-400];  
% obj(8,1:4)=[19, 252, yl-190, yl-74];
% obj(9,1:4)=[257, 420, yl-190, yl-27];
% obj(10,1:4)=[426, 589, yl-190, yl-27];
% obj(11,1:4)=[595, 758, yl-190, yl-27];
% obj(12,1:4)=[770, 1005, yl-190, yl-27];
% obj(13,1:4)=[19, 309, yl-523, yl-196];
% obj(14,1:4)=[311, 531, yl-523, yl-196];
% obj(15,1:4)=[533,753, yl-523, yl-196];
% obj(16,1:4)=[755, 1004, yl-523, yl-196];
yl = 768;
obj(1,1:4)=[19, 227, yl-752, yl-510];
obj(2,1:4)=[234, 774, yl-678, yl-510];
obj(3,1:4)=[781, 1004, yl-752, yl-510];
obj(4,1:4)=[106, 101+126, yl-400-96, yl-400];
obj(5,1:4)=[330,325+126, yl-400-96, yl-400];
obj(6,1:4)=[556, 551+126, yl-400-96, yl-400];
obj(7,1:4)=[781, 776+126, yl-400-96, yl-400];
obj(8,1:4)=[19, 252, yl-190, yl-27];
obj(9,1:4)=[257, 420, yl-190, yl-27];
obj(10,1:4)=[426, 589, yl-190, yl-27];
obj(11,1:4)=[595, 758, yl-190, yl-27];
obj(12,1:4)=[770, 1005, yl-190, yl-27];
obj(13,1:4)=[19, 309, yl-508, yl-196];
obj(14,1:4)=[311, 531, yl-508, yl-196];
obj(15,1:4)=[533,753, yl-508, yl-196];
obj(16,1:4)=[755, 1004, yl-508, yl-196];

for i = 46:length(nam_save)  % 393 participants
%     i=393-j;
%     close all
    swipe=swipe_save{i};
    skip=0;
    
    file_id = [nam_save{i},'.mat']; %[nam_save{i},'/',nam_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,'\',file_id])   % check file exists
        
        idc = strfind(subject_details(:,1),nam_save{i});
        idx = find(~cellfun('isempty',idc));
        cat = subject_details{idx,5};
        
        figure
        for ii = 1 : 16
            hold on
%             scatter(obj(i,1),obj(i,3),'+')
            hold on
%             scatter(obj(i,2),obj(i,4),'+')
            hold on
%             scatter(obj(i,1),obj(i,4),'+')
            hold on
%             scatter(obj(i,2),obj(i,3),'+')
            plot([obj(ii,1),obj(ii,2)],[obj(ii,3),obj(ii,3)],'k')
            plot([obj(ii,1),obj(ii,2)],[obj(ii,4),obj(ii,4)],'k')
            plot([obj(ii,1),obj(ii,1)],[obj(ii,3),obj(ii,4)],'k')
            plot([obj(ii,2),obj(ii,2)],[obj(ii,3),obj(ii,4)],'k')
        end
        axis([0 1100 0 800])
%         if strcmp(cat,'TD')
%             colr = 'b';
%             sname = ['TD ',nam_save{i}];
%             title(sname)
%         else
        if strcmp(cat,'ASD')
            colr = 'c';
            sname = [cat,' ',nam_save{i}];
            title(sname)
        elseif strcmp(cat,'TD')
            colr = 'g';
            sname = [cat,' ',nam_save{i}];
            title(sname)
        end
            
        titlename = ['ID ',nam_save{i}];
        filename = [file_loc,'/',file_id];
        savename = ['subject_',nam_save{i}];

        load(filename,'subject_table');
        tab=sortrows(subject_table,8);
        
        if iscell(tab.x(1))                         % convert strings to doubles for X and Y
            tab.x=str2double(strrep(tab.x,',','.'));
            tab.y=str2double(strrep(tab.y,',','.'));
        end
        adj = zeros(17,17);     % 17 defined zones/nodes

        prev_tP =0;
        %% Plot swipes
        for m = 1 : length(swipe)   % process each swipe
            %% start from the back
            nn=[];
            n=1;

            k = swipe{m};

            n=n+1;

            hold on
            plot(tab.x(k),tab.y(k),'r',LineWidth=2)
            plot(tab.x(k),tab.y(k),colr,LineWidth=2)
            pause(.10)
        end
        for m = 1 : length(swipe)
            %% start from the back
            nn=[];
            n=1;

            k = swipe{m};

            n=n+1;

            hold on
            scatter(tab.x(k(1)),tab.y(k(1)),'k.')
            scatter(tab.x(k(end)),tab.y(k(end)),'ko')
        end
        axis off
        axis equal
        saveas(gcf,['Accurate_pngs/',sname,'.png'])
%         end
    end
    
end