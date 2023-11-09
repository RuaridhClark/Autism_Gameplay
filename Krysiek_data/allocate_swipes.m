%% Identify each individual swipe

clear all

load('tC_save_Krysiek_ond.mat','tC_save','nam_save') % load each swipe identified 

swipe_save = cell(1,length(nam_save));
for i = 1:length(nam_save)
    tC=tC_save{i};
    skip=0;
    
    file_id = [nam_save{i}];

    swipe = cell(1,max(tC));
    for ii = 1 : max(tC)
        swipe{ii}=find(tC==ii);
    end

    swipe_save{i}=swipe;
end