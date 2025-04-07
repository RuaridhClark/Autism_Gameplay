clear all
folder1 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Data';
addpath(folder1)
% file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';
load('..\Create_adj\swipes_trial.mat')

load('subject_details.mat','subject_details_776')

zone = defineZonePositions();

for i = 428%1:704   % 704 participants
%     i=705-j;
    close all
    swipe=swipe_save{i};
    skip=0;
    
    file_id = [name_save{i},'/',name_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,file_id])   % check file exists
        
        idc = strfind(subject_details_776(:,1),name_save{i});
        idx = find(~cellfun('isempty',idc));
        cat = subject_details_776{idx,5};
        
        figure
        for ii = 1 : 16
            hold on
%             scatter(zone(i,1),zone(i,3),'+')
            hold on
%             scatter(zone(i,2),zone(i,4),'+')
            hold on
%             scatter(zone(i,1),zone(i,4),'+')
            hold on
%             scatter(zone(i,2),zone(i,3),'+')
            plot([zone(ii,1),zone(ii,2)],[zone(ii,3),zone(ii,3)],'k')
            plot([zone(ii,1),zone(ii,2)],[zone(ii,4),zone(ii,4)],'k')
            plot([zone(ii,1),zone(ii,1)],[zone(ii,3),zone(ii,4)],'k')
            plot([zone(ii,2),zone(ii,2)],[zone(ii,3),zone(ii,4)],'k')
        end
        axis([0 1100 0 800])
%         if strcmp(cat,'TD')
%             colr = 'b';
%             sname = ['TD ',name_save{i}];
%             title(sname)
%         else
        if strcmp(cat,'ASD')
            colr = 'c';
            sname = [cat,' ',name_save{i}];
            title(sname)
        elseif strcmp(cat,'TD')
            colr = 'g';
            sname = [cat,' ',name_save{i}];
            title(sname)
        end
            
        titlename = ['ID ',name_save{i}];
        filename = [file_loc,file_id];
        savename = ['subject_',name_save{i}];

        tab=sortrows(readtable(filename),'Time');   % sort table according to time
        if iscell(tab.X(1))                         % convert strings to doubles for X and Y
            tab.X=str2double(strrep(tab.X,',','.'));
            tab.Y=str2double(strrep(tab.Y,',','.'));
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
            % plot(tab.X(k),tab.Y(k),'r',LineWidth=2)
            p = plot(tab.X(k),tab.Y(k),colr,LineWidth=2);
            delete(p)
            pause(.10)
        end
        for m = 48 : length(swipe)
            nn=[];
            n=1;

            k = swipe{m};

            n=n+1;

            hold on
            scatter(tab.X(k(1)),tab.Y(k(1)),'k.')
            scatter(tab.X(k(end)),tab.Y(k(end)),'ko')
        end
        axis off
        axis equal
        saveas(gcf,['Accurate_pngs/',sname,'.png'])
%         end
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