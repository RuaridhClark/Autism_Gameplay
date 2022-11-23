% check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
folder4 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Set_allocate';
% folder5 = 'H:\My Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj';
folder7 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
addpath(folder4,folder6,folder7)
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_obj_end_accurate\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =16;    % number of ipad zones (nodes)
redirect = 1; % rewire snap-to-target zones
pert_init=-.80;

saved = zeros(num,704);
ranked = zeros(704,1);


for i = 1:704
    file_id = ['subject_',nam_save{i},'.mat'];

    j=0;
    plot_data=[];

    if isfile([file_loc,file_id])
        load([file_loc,file_id])
        titlename = ['ID ',nam_save{i}];
        savename = ['subject_',nam_save{i}];

        if num == 16 && redirect == 0
            L = adj2L(adj,num);
        elseif num == 16 && redirect == 1
%             L = adj2L_snap2zones(adj,num);
            L = adj2L_snap2zones_foodloc(adj,num);
        end
        
        list = [4,5,6,7]; check=1; 
        pert =40*.01; prev=pert_init; tmp_pert=pert;
        for pert = -0.1 : 0.1 : 0.1
            vec=zeros(1,num);
            vec(list)=ones(1,length(list)).*pert;
            P = -(L+diag(vec));

            %% Save the sorted first eigenvector entries
            [V,D]=eig(P');
            [~,I]=sort(diag(real(D)),'desc');
            [~,II]=sort(abs(V(:,I(1))),'desc');

            saved(:,i)=II;
            
            j=j+1;
            plot_data(j,:)=[pert,V(:,I(1))'];
            holddata=plot_data(j,14:17);
            plot_data(j,14:17)=plot_data(j,9:12);
            plot_data(j,9:12)=holddata;

            [check,tmp_pert] = check_topfour(saved,check,i,tmp_pert,list,pert);
            
        end
        
        %%%%% plot eigenvectors
        
        leglist = [];
        figure;
        cmap = colormap(hot);
        step=60;%round(length(cmap)/size(plot_data,1));
        for j = 1 : size(plot_data,1)
        hold on
        h=plot(1:16,abs(plot_data(j,2:end)),'-o','Color',cmap(j*step,:),'MarkerFaceColor',cmap(j*step,:));
%         h.MarkerFaceColor = h.Color;
        leglist = [leglist;string(strcat('{\it p_i} ='," ",num2str(plot_data(j,1))))];
        end
        axis([.5 16.5 .14 .45])
        xticks([1:16])
        legend(leglist)
        ylabel('{\bf v}_1 entries','fontsize', 14)
        xlabel('vertices','fontsize', 14)
    end
end


%%%%%%%%% functions %%%%%%%%%
function [L] = adj2L(adj,num)
    %% remove zn 4-7 incoming except from 2, 4, 5, 6, 7
    allow=[2];%,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            if num == 16
                adj(it,13:16)=adj(it,4:7)+adj(it,13:16);    % reconnect to 13-16
            end
            adj(it,4:7)=zeros(1,4);                     % remove non-food connections
        end
    end
%     adj=adj-diag(diag(adj));      % remove diagonal
    
    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));   % convert to proportional weights
    end

    bweight=.01;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L=-adj + diag(sum(adj,2));

    %% Remove diagonal - convert L into adj (sort of)
    L=L-diag(diag(L)); 
end

% function [L] = adj2L_snap2zones(adj,num)
%     for it = 1:num
%         adj(it,4:7)=adj(it,4:7)+adj(it,13:16);    % reconnect to 13-16
%         adj(4:7,it)=adj(4:7,it)+adj(13:16,it);    % reconnect to 13-16
%         adj(it,13:16)=zeros(1,4);                     % remove non-food connections
%         adj(13:16,it)=zeros(4,1); 
%     end
%     adj = adj(1:12,1:12);
%     %% remove zn 4-7 incoming except from 2
%     allow=[2,4,5,6,7];
%     for it = 1:num
%         if ~ismember(it,allow)
%             adj(it,4:7)=zeros(1,4);                     % remove non-food connections
%         end
%     end
% %     adj=adj-diag(diag(adj));  % remove diagonal
%     
%     if sum(adj(:))>0
%         adj = (adj./sum(adj(:)));%.*(100); %%%%%%%%%%%%%%%%%%%%%%% TEMP ADDITION Normalising
%     end
% 
%     bweight=.01;
%     [adj] = NNR_adj_conns_OBJ2(adj,bweight);
% 
%     L=-adj + diag(sum(adj,2));
% 
%     %% convert L into adj (sort of)
%     L=L-diag(diag(L)); 
% end

function [L] = adj2L_snap2zones_foodloc(adj,num)
    % 2 to 4-7 = 2 to 4-7 13-16
    % 4-7 
    adj(2,4:7)=adj(2,4:7)+adj(2,13:16);    % reconnect 2 to 13-16
%     adj(4:7,4:7)=adj(4:7,4:7)+adj(4:7,13:16);    % reconnect 2 to 13-16
    adj(2,13:16)=zeros(1,4);               % remove re-connected connections
%     adj(4:7,13:16)=zeros(4,4);               % remove re-connected connections

    %% remove zn 4-7 incoming except from 2
    allow=[2];%,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            if num == 16
                adj(it,13:16)=adj(it,4:7)+adj(it,13:16);    % reconnect to 13-16
            end
            adj(it,4:7)=zeros(1,4);                     % remove non-food connections
        end
    end
    
    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));%.*(100); %%%%%%%%%%%%%%%%%%%%%%% TEMP ADDITION Normalising
    end

    bweight=.01;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L=-adj + diag(sum(adj,2));

    %% Remove diagonal - convert L into adj (sort of)
    L=L-diag(diag(L)); 
end

function [check,tmp_pert] = check_topfour(saved,check,i,tmp_pert,list,pert)
    if max(saved(:,i))>0
        if max(find(ismember(saved(:,i),list)==1))>4 
            if check == 10
                check = 0;
            elseif check == -10
                check = -10;
            else
                check = -1;
            end
            tmp_pert = pert;
        else
            if check == -10
                check = 0;
            elseif check == 10
                check = 10;
            else
                check = 1;
            end
        end
    end
end

function [pert,prev,check] = pert_bisect(pert,prev,pert_init,tmp_pert,check)
%     pert = pert-pert_init;
    dif = pert-prev;
%             [saved(1:4)',check]
    prev=pert;
    if check == 1
        pert = pert + abs(dif)/2;
    elseif check == -1
        pert = pert - abs(dif)/2;
    end

    if abs(dif)/2 < 0.5*.01
        if check == 1
            check = 10;
        elseif check == -1
            check = -10;
        end
        pert=tmp_pert;%-pert_init;
    end
end

function [pert] = pert_iter(pert,check,pert_init)
%     pert = pert-pert_init;
    if check == 10
        iter = 0.001;
        pert = pert + iter;
    elseif check == -10
        iter = -0.001;
        pert = pert + iter;
    end
end
