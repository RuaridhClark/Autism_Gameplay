clear all
folder1='H:\My Documents\MATLAB\Autism_MAIN\Create_adj';
addpath(folder1)
% file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';
load('swipes_all704.mat')

for i = 1:704
    close all
    swipe=swipe_save{i};
    skip=0;
    
    file_id = [nam_save{i},'/',nam_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,file_id])
        colr = 'b';
        titlename = ['ID ',nam_save{i}];
        filename = [file_loc,file_id];
        savename = ['subject_',nam_save{i}];
    else
        skip = 1;
    end
    
    if skip == 0
        tab=sortrows(readtable(filename),4); % sort table according to time
        if iscell(tab.X(1)) %% Convert strings to doubles for X and Y
            tab=sortrows(readtable(filename),14); % sort table according to time
            tab.X=strrep(tab.X,',','.');
            tab.Y=strrep(tab.Y,',','.');
            tab.X=str2double(tab.X);
            tab.Y=str2double(tab.Y);
        end
        adj = zeros(17,17);
              
        prev_tP =0;
        for m = 1 : length(swipe)
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