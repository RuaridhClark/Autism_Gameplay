% check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
% folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
% folder2 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate\functions';
folder3 = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate';
folder4 = 'H:\My Documents\GitHub\Autism_Gameplay\Set_allocate';
folder5 = 'H:\My Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'H:\My Documents\GitHub\Autism_Gameplay\Create_adj_110721';
folder7 = 'H:\My Documents\GitHub\Autism_Gameplay';
addpath(folder3,folder4,folder5,folder6,folder7)
file_loc = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =16;    % number of ipad zones (nodes)
pert_init=-80;

saved = zeros(num,704);
ranked = zeros(704,1);

for i = 1:704
    file_id = ['subject_',nam_save{i},'.mat'];

    if isfile([file_loc,file_id])
        load(file_id)
        titlename = ['ID ',nam_save{i}];
        savename = ['subject_',nam_save{i}];

        L = adj2L(adj,num);
        
        list = [4,5,6,7]; check=1; 
        pert =40; prev=pert_init; tmp_pert=pert;
        while check ~= 0
            vec=zeros(1,num);
            vec(list)=ones(1,length(list)).*pert;
            P = -(L+diag(vec));

            %% Save the sorted first eigenvector entries
            [V,D]=eig(P');
            [~,I]=sort(diag(real(D)),'desc');
            [~,II]=sort(abs(V(:,I(1))),'desc');

            saved(:,i)=II;

            [check,tmp_pert] = check_topfour(saved,check,i,tmp_pert,list,pert);
            
            if check == 1 || check == -1
                [pert,prev,check] = pert_bisect(pert,prev,pert_init,tmp_pert,check);
            elseif check == 10 || check == -10
                [pert] = pert_iter(pert,check,pert_init);
            end
        end
        if ranked(i)==0
            ranked(i)=pert;%-pert_init;
        end
    else 
        ranked(i)=0.5;
    end
end

save('temp_save_OBJ_end.mat')

%%%%%%%%% functions %%%%%%%%%
function [L] = adj2L(adj,num)
    %% remove zn 4-7 incoming except from 2
    allow=[2,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            adj(it,13:16)=adj(it,4:7)+adj(it,13:16);    % reconnect to 13-16
            adj(it,4:7)=zeros(1,4);                     % remove non-food connections
        end
    end
    adj=adj-diag(diag(adj));
    
    bweight=1;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L=-adj + diag(sum(adj,2));

    %% convert L into adj (sort of)
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
