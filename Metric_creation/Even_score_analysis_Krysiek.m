% check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all

folder4 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Set_allocate';
% folder5 = 'H:\My Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj';
folder7 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
addpath(folder4,folder6,folder7)
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_foodpile\'; %\adj_krysiek\'; % should match zone type

load('swipes_Krysiek_ond.mat','nam_save')

%% stack the adjs
num = 16;    % number of ipad zones (nodes)
redirect = 1; % rewire snap-to-target zones (accurate 0 or snap-to 1)
pert_init=-.80;

lenN = length(nam_save);

saved = zeros(num,lenN);
ranked = zeros(lenN,1);

type = 'plates';
label = '';
bweight = 0.01;

for i = 1:lenN
    file_id = ['subject_',nam_save{i},'.mat'];

    if isfile([file_loc,file_id])
        load([file_loc,file_id])
        titlename = ['ID ',nam_save{i}];
        savename = file_id;

        if num == 16 && redirect == 0
            L = adj2L(adj,num);
        elseif num == 16 && redirect == 1
%             L = adj2L_snap2zones(adj,num);
            if strcmp(type,'plates')
                L = adj2L_snap2zones_foodloc(adj,num,label,bweight);  % direct food delivery
            elseif strcmp(type,'inter')
                L = adj2L_interplate(adj,num,bweight);    % inter-plate analysis
            end
        end
        
        list = [4,5,6,7]; check=1; 
        pert =40*.01; prev=pert_init; tmp_pert=pert;
        while check ~= 0
            vec=zeros(1,num);
            P = -(L+diag(vec));

            %% Save the sorted first eigenvector entries
            [V,D]=eig(P');
            [~,I]=sort(diag(real(D)),'desc');
            [~,II]=sort(abs(V(:,I(1))),'desc');

            saved(:,i)=II;

            %% Save the sorted first eigenvector entries
            [V,D]=eig(P');
            [~,I]=sort(diag(real(D)),'desc');
            [~,II]=sort(abs(V(:,I(1))),'desc');
    
            pert=std(abs(V(4:7,I(1))))/sum(abs(V(4:7,I(1))));
            check=0;
%             pert
        end
        if ranked(i)==0
            ranked(i)=pert;%-pert_init;
        end
    else 
        ranked(i)=0.5*.01;
    end
end

if strcmp(type,'plates')
    if strcmp(label,'accurate')
        save('test_SS_accurate_Krysiek.mat','ranked','nam_save')
    else
        save('Even_ext_Krysiek.mat','ranked','nam_save')
    end
elseif strcmp(type,'inter')
    save('Even_ext_inter_Krysiek.mat','ranked','nam_save')
end
% save('temp_save_OBJ_end.mat','ranked','nam_save')

%%%%%%%%% functions %%%%%%%%%
function [L] = adj2L(adj,num)
    %% remove zn 4-7 incoming except from 2
    allow=[2,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            if num == 16
                adjust = adj(it,4:7);
                adj(it,13:16)=adj(it,13:16)+adjust;    % reconnect to 13-16
%                 if ismember(it,4:7)
%                     adj(it,it+9)=adj(it,it+9)-adjust(ismember(4:7,it));
%                 end
            end
            list = 4:7;
            listedit =list(~ismember(4:7,it));
            adj(it,listedit)=zeros(1,length(listedit));                     % remove non-food connections
        end
    end   

      
    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));   % convert to proportional weights
    end

    bweight=.01;
%     bweight = sum(adj(:))/length(adj)^2;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L=-adj + diag(sum(adj,2));

    %% Remove diagonal - convert L into adj (sort of)
    L=L-diag(diag(L)); 
end

function [L] = adj2L_snap2zones_foodloc(adj,num,label,bweight)
    % 2 to 4-7 = 2 to 4-7 13-16
    % 4-7 
    if ~strcmp(label,'accurate')
        adj(2,4:7)=adj(2,4:7)+adj(2,13:16);    % Collect all food delivery swipe
        adj(2,13:16)=zeros(1,4);               % remove re-connected connections
    end

    %% remove zn 4-7 incoming except from 2
    allow=[2];
    for it = 1:num
        if ~ismember(it,allow)
            if num == 16
                adjust = adj(it,4:7);
                adj(it,13:16)=adj(it,13:16)+adjust;    % reconnect to 13-16
%                 if ismember(it,4:7)
%                     adj(it,it+9)=adj(it,it+9)-adjust(ismember(4:7,it)); % remove plate self-loop
%                 end
            end
            list = 4:7;
            listedit =list(~ismember(4:7,it));
            adj(it,listedit)=zeros(1,length(listedit));                     % remove non-food connections
        end
    end
    
    adj = adj-diag(diag(adj));

%     % reroute connections
%     adj(:,2)=adj(:,2)+sum(adj,1)';
%     adj(2,2)=0;

    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));% Normalising
    end

%     bweight=.01;
%     bweight = 1/length(adj)^2;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L = -adj;
end

function [L] = adj2L_interplate(adj,num,bweight)

    %% remove connections to zn 4-7 except from food zone and other plates
    allow=[2,4,5,6,7];
    % Snap-to-plate 
    for i = 1:length(allow)
        adj(allow(i),4:7)=adj(allow(i),4:7)+adj(allow(i),13:16);    % reconnect 2 to 13-16
        dest = 13:16;
        vals = dest(~ismember(dest,allow(i)+9));
        adj(allow(i),vals)=zeros(1,length(vals));               % remove re-connected connections
%         adj(allow(i),13:16)=zeros(1,4);               % remove re-connected connections
    end
    % Remove direct to plate
    for it = 1:num
        if ~ismember(it,allow)
            if num == 16
                % remove self-loop rellocation (e.g. 4,4 adding to 4,13)
                adj(it,13:16)=adj(it,4:7)+adj(it,13:16);    % reconnect to 13-16
            end
            list = 4:7;
            listedit =list(~ismember(4:7,it));
            adj(it,listedit)=zeros(1,length(listedit));                     % remove non-food connections
        end
    end
    
    adj=adj-diag(diag(adj));
    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));% Normalising
    end

%     bweight=.01;
%     bweight = 1/length(adj)^3;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L = -adj;%+diag(diag(adj));% + diag(sum(adj,2));

    %% Remove diagonal - convert L into adj (sort of)
%     L=L-diag(diag(L)); 
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
