%% Identify each individual swipe

clear all
file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';

load('tCs_all704.mat','tC_save','nam_save') % load each swipe identified 

m = 0; n = m;
swipe_save = cell(1,704);
for i = 1:704
%     close all
    tC=tC_save{i};
    skip=0;
    
    file_id = [nam_save{i},'/',nam_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,file_id])
        n=n+1;
        filename = [file_loc,file_id];
        savename = ['subject_',nam_save{i}];

        swipe = cell(1,max(tC));
        for ii = 1 : max(tC)
            swipe{ii}=find(tC==ii);
        end
    end
    swipe_save{i}=swipe;
end