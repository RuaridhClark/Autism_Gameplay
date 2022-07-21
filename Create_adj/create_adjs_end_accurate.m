% close all
clear all
file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';
load('swipes_all704.mat')

file_adj = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\';

yl = 768;
obj(1,1:4)=[19, 227, yl-752, yl-525];
obj(2,1:4)=[234, 774, yl-678, yl-525];
obj(3,1:4)=[781, 1004, yl-752, yl-525];
obj(4,1:4)=[101, 101+131, yl-400-96, yl-400]; 
obj(5,1:4)=[325,325+131, yl-400-96, yl-400];  
obj(6,1:4)=[551, 551+131, yl-400-96, yl-400]; 
obj(7,1:4)=[776, 776+131, yl-400-96, yl-400];   
% obj(9,1:4)=[113, 113+112, yl-283-84, yl-283];
% obj(10,1:4)=[340, 340+112, yl-283-84, yl-283];
% obj(11,1:4)=[566, 566+112, yl-283-84, yl-283];
% obj(12,1:4)=[790, 790+112, yl-283-84, yl-283];
obj(8,1:4)=[19, 252, yl-190, yl-74];
obj(9,1:4)=[257, 420, yl-190, yl-27];
obj(10,1:4)=[426, 589, yl-190, yl-27];
obj(11,1:4)=[595, 758, yl-190, yl-27];
obj(12,1:4)=[770, 1005, yl-260, yl-52];
obj(13,1:4)=[19, 309, yl-523, yl-196];
obj(14,1:4)=[311, 531, yl-523, yl-196];
obj(15,1:4)=[533,753, yl-523, yl-196];
obj(16,1:4)=[755, 1004, yl-523, yl-264];

for i = 1:704
    close all
    swipe=swipe_save{i};
    skip=0;
    
    file_id = [nam_save{i},'\',nam_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,file_id])
        colr = 'b';
        titlename = ['ID ',nam_save{i}];
        filename = [file_loc,file_id];
        savename = ['subject_',nam_save{i}];
    else
        skip = 1;
    end
    
    %% Load tab missing!

    if skip == 0
        tab=sortrows(readtable(filename),4); % sort table according to time
        if iscell(tab.X(1)) %% Convert strings to doubles for X and Y
            tab=sortrows(readtable(filename),14); % sort table according to time
            tab.X=strrep(tab.X,',','.');
            tab.Y=strrep(tab.Y,',','.');
            tab.X=str2double(tab.X);
            tab.Y=str2double(tab.Y);
        end
        adj = zeros(16,16);
              
        prev_tP =0;
        for m = 1 : length(swipe)
            prev_n = [];
            n=1;
            while isempty(prev_n) && n<=length(swipe{m}) 
                k = swipe{m}(n);
                nn=[];
                for jj = 1 : 16
                    ii = 17-jj;
                    if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                    	nn = ii; 
                    end
                end
                    
                if n>1 && ~isempty(prev_n) && ~isempty(nn)
%                     adj(prev_n,nn)=adj(prev_n,nn)+1;
% %                     prev_n = nn;
                elseif isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
                n=n+1;
            end
            %% start from the back
            nn=[];
            n=-1;
            while isempty(nn) && ~isempty(prev_n)
                n=n+1;
                k = swipe{m}(end-n);
                
                for jj = 1 : 16
                    ii = 17-jj;
                    if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                    	nn = ii; 
                    end
                end
                    
                if ~isempty(prev_n) && ~isempty(nn)
                    adj(prev_n,nn)=adj(prev_n,nn)+1;
%                     prev_n = nn;
                elseif isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
            end
        end
    end
    
    %% save adj
    save([file_adj,'adj_obj_end_accurate2\',savename],'adj') 
    
%     %% plot adj
%     figure;spy(adj)
%     title(titlename)
%     saveas(gcf,['spy',savename,'.png'])
end