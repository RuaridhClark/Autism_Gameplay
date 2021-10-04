clear all
folder1 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj_110721';
addpath(folder1)
file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';
load('swipes_all704.mat')

for i = 1:704   % 704 participants
    close all
    swipe=swipe_save{i};
    skip=0;
    
    file_id = [nam_save{i},'/',nam_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,file_id])   % check file exists
        colr = 'b';
        titlename = ['ID ',nam_save{i}];
        filename = [file_loc,file_id];
        savename = ['subject_',nam_save{i}];
        
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
            plot(tab.X(k),tab.Y(k),colr)
        end
    end
    

end