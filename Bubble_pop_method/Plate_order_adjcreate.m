% close all
clear all
list =[];
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';
load('swipes_all704.mat')

file_adj = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Bubble_pop_method\adjs\';

yl = 768;
obj(1,1:4)=[19, 227, yl-752, yl-525];
obj(2,1:4)=[234, 774, yl-678, yl-525];
obj(3,1:4)=[781, 1004, yl-752, yl-525];
obj(4,1:4)=[106, 101+126, yl-400-86, yl-400]; 
obj(5,1:4)=[330,325+126, yl-400-86, yl-400];  
obj(6,1:4)=[556, 551+126, yl-400-86, yl-400]; 
obj(7,1:4)=[781, 776+126, yl-400-86, yl-400];  
obj(8,1:4)=[19, 252, yl-190, yl-74];
obj(9,1:4)=[257, 420, yl-190, yl-27];
obj(10,1:4)=[426, 589, yl-190, yl-27];
obj(11,1:4)=[595, 758, yl-190, yl-27];
obj(12,1:4)=[770, 1005, yl-190, yl-27];
obj(13,1:4)=[19, 309, yl-523, yl-196];
obj(14,1:4)=[311, 531, yl-523, yl-196];
obj(15,1:4)=[533,753, yl-523, yl-196];
obj(16,1:4)=[755, 1004, yl-523, yl-196];

zone2_swipes = {};

choice = 'plates'; % 'plates' or 'inter'

if strcmp(choice,'plates')
    optns = 2;
else
    optns = [2,4,5,6,7];
end

for i = 1:335 %562
%     close all
    swipe=swipe_save{i};
    skip=0;
    save_n = [];
    
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

        % swipe start
        prev_tP =0;
        swipe_lens = [];
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
                    
                if isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
                n=n+1;
            end
            nn=[];
            if ismember(prev_n,optns) % switch between plates [2] and inter [2,4,5,6,7]
                %% swipe end: start from the back
                nn=[];
                n=-1;
                while isempty(nn)
                    n=n+1;
                    k = swipe{m}(end-n);
                    
                    for jj = 1 : 16
                        ii = 17-jj;
                        if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                	        nn = ii; 
                        end
                    end

                    if ismember(nn,[13,14,15,16])   % snap-to-plate zones converted to their plate zone IDs
                        nn = nn - 9;
                    end
                            
                    if ~isempty(save_n) && ~isempty(nn)
                        if ismember(nn,[4,5,6,7]) && save_n ~= nn && prev_n ~= nn   % stop plate taps (self-loops) being recorded
                            adj(save_n,nn)=adj(save_n,nn)+1;        % add to adj
                        end
                    end
                end
                if ismember(prev_n,2) && ismember(nn,[4,5,6,7])% 
                    save_n = nn;
                elseif ismember(prev_n,[4,5,6,7]) % 
                    save_n = prev_n;
                end
            end
        end
    end
    
    zone2_swipes{i} = swipe_lens;

    if size(adj,1)>4 && sum(adj(:))>0
    
        adj = adj(4:7,4:7);
    
        %% Save the sorted first eigenvector entries
%         [V,D]=eig(adj');
%         [~,I]=sort(diag(real(D)),'desc');
%         [~,II]=sort(abs(V(:,I(1))),'desc');
    
         %std(V(:,I(1)));
    
        DG = digraph(adj);
    
        figure
        LWidths = 5*DG.Edges.Weight/max(DG.Edges.Weight);
        p = plot(DG,'LineWidth',LWidths);

        adj = adj./sum(adj(:));
        pert = max(min(max(adj,[],2)));
        title(strcat('min. edge =  ',num2str(round(pert,3))))

        for z = 1:4
            labelnode(p,z,num2str(z+3))
        end
        axis off
        axis equal
    
        saveas(gcf,['Graphs/',choice,'_',savename,'.png'])
        close gcf
    end

%     %% save adj
%     if strcmp(choice,'plates')
%         save([file_adj,'plates_order\',savename],'adj') 
%     else
%         save([file_adj,'inter_order\',savename],'adj') 
%     end
    
%     %% plot adj
%     figure;spy(adj)
%     title(titlename)
%     saveas(gcf,['spy',savename,'.png'])
end