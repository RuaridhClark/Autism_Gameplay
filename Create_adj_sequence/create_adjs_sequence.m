% close all
clear all
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';

addpath('C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj')

load('swipes_all704.mat')

file_adj = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\';

load('sets.mat')

yl = 768;
obj(1,1:4)=[19, 227, yl-752, yl-525];
obj(2,1:4)=[234, 774, yl-678, yl-525];
obj(3,1:4)=[781, 1004, yl-752, yl-525];
obj(4,1:4)=[101, 101+131, yl-400-96, yl-400]; 
obj(5,1:4)=[325,325+131, yl-400-96, yl-400];  
obj(6,1:4)=[551, 551+131, yl-400-96, yl-400]; 
obj(7,1:4)=[776, 776+131, yl-400-96, yl-400];   
obj(8,1:4)=[19, 252, yl-190, yl-74];
obj(9,1:4)=[257, 420, yl-190, yl-27];
obj(10,1:4)=[426, 589, yl-190, yl-27];
obj(11,1:4)=[595, 758, yl-190, yl-27];
obj(12,1:4)=[770, 1005, yl-260, yl-52];
obj(13,1:4)=[19, 309, yl-523, yl-196];
obj(14,1:4)=[311, 531, yl-523, yl-196];
obj(15,1:4)=[533,753, yl-523, yl-196];
obj(16,1:4)=[755, 1004, yl-523, yl-264];

num = 16;
for i = 1:704
%     close all
%     %%%
%     figure
%     for ii = 1 : 16
%         hold on
% %             scatter(obj(i,1),obj(i,3),'+')
%         hold on
% %             scatter(obj(i,2),obj(i,4),'+')
%         hold on
% %             scatter(obj(i,1),obj(i,4),'+')
%         hold on
% %             scatter(obj(i,2),obj(i,3),'+')
%         plot([obj(ii,1),obj(ii,2)],[obj(ii,3),obj(ii,3)],'k')
%         plot([obj(ii,1),obj(ii,2)],[obj(ii,4),obj(ii,4)],'k')
%         plot([obj(ii,1),obj(ii,1)],[obj(ii,3),obj(ii,4)],'k')
%         plot([obj(ii,2),obj(ii,2)],[obj(ii,3),obj(ii,4)],'k')
%     end
%     axis([0 1100 0 800])
%     %%%
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
        adj = zeros(num,num);

        dest_array = zeros(2,1);
        iter = 1;

        prev_tP =0;
        save_n = [];
        for m = 1 : length(swipe)
            %% Find start of swipe
            prev_n = [];
            n=1;
            while isempty(prev_n) && n<=length(swipe{m}) 
                k = swipe{m}(n);
                nn=[];
                for jj = 1 : num
                    ii = num+1-jj;
                    if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                    	nn = ii; 
                    end
                end
                    
                if isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
                n=n+1;
            end

            %% Find end of swipe - start from the back
            nn=[];
            n=-1;
            while isempty(nn) && ~isempty(prev_n)
                n=n+1;
                k = swipe{m}(end-n);
                z = swipe{m};
                
                for jj = 1 : num
                    ii = num+1-jj;
                    if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                    	nn = ii; 
                    end
                end
                    
                if ~isempty(prev_n) && ~isempty(nn)
                    %% Define destination array
                    if prev_n == 2 && (ismember(nn,4:7) || ismember(nn,13:16)) && prev_n~=nn
%                         adj(prev_n,nn)=adj(prev_n,nn)+1;
                        dest_array(1,iter) = prev_n;
                        if nn>7
                            nn = nn-9;
                        end
                        dest_array(2,iter) = nn;
                        iter=iter+1;
                    elseif ismember(prev_n,4:7) && (ismember(nn,4:7) || ismember(nn,13:16)) && prev_n~=nn
%                         adj(prev_n,nn)=adj(prev_n,nn)+1;
                        dest_array(1,iter) = prev_n;
                        if nn>7
                            nn = nn-9;
                        end
                        dest_array(2,iter) = nn;
                        iter=iter+1;
                    end
%                     %%%
%                     hold on
%                     plot(tab.X(z),tab.Y(z),'r',LineWidth=2)
%                     plot(tab.X(z),tab.Y(z),colr,LineWidth=2)
%                     %%%
                    save_n = nn;
                elseif isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
            end
        end
    end

    %% Use destination array to define adjacency
    prev_n = 2;
    set = [];
    for j = 1 : size(dest_array,2)
        nn = dest_array(2,j);
        if prev_n ~= nn
            adj(prev_n,nn)=adj(prev_n,nn)+1;
            set = [set;nn];
            if min(ismember(4:7,set))==1
                set = [];
                prev_n = 2;
                % Connection back to the food to highlight a cycle complete
                adj(nn,prev_n)=adj(nn,prev_n)+1;
            else
                prev_n = nn;
            end
        end
    end
    
    x = mean(obj(:,1:2),2);
    y = mean(obj(:,3:4),2);
    
    y(5:6) = y(5:6)+200;
    
    G=digraph(adj);
    G = rmnode(G,[1,3,8:16]);
    x = x([2,4:7]);
    y = y([2,4:7]);
    figure;plot(G,'XData',x,'YData',y,'LineWidth',G.Edges.Weight) % 'EdgeLabel',G.Edges.Weight,

    if ismember(i,sets{1})
        titlename = strcat(titlename,' - ','F TD');
        savename = strcat(savename,'_F_TD');
    elseif ismember(i,sets{5})
        titlename = strcat(titlename,' - ','M TD');
        savename = strcat(savename,'_M_TD');
    elseif ismember(i,sets{2})
        titlename = strcat(titlename,' - ','F ASD');
        savename = strcat(savename,'_F_ASD');
    elseif ismember(i,sets{6})
        titlename = strcat(titlename,' - ','M ASD');
        savename = strcat(savename,'_M_ASD');
    end

    title(titlename)

    %% save adj
    save([file_adj,'adj_sequence\',savename],'adj') 
    
%     %% plot adj
%     figure;spy(adj)
%     title(titlename)
    saveas(gcf,['Graph_plots/',savename,'.png'])
    close gcf
end