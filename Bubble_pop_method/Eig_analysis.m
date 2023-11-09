%%% Sharing Score Analysis

clear all
folder4 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Set_allocate';
% folder5 = 'H:\My Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj';
folder7 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
addpath(folder4,folder6,folder7)
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Bubble_pop_method\adjs\inter_noself\'; %\adj_obj_end_accurate\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =16;    % number of ipad zones (nodes)
redirect = 1; % rewire snap-to-target zones (accurate 0 or snap-to 1)
pert_init=-.80;

saved = zeros(num,704);
ranked = zeros(704,1);

type = 'inter';
label = '';

for i = 1:704
    file_id = ['subject_',nam_save{i},'.mat'];

    if isfile([file_loc,file_id])
        load([file_loc,file_id])
        titlename = ['ID ',nam_save{i}];
        savename = ['subject_',nam_save{i}];
        
        adj = adj(4:7,4:7);
%         adj = adj - diag(diag(adj));

%         red_adj = zeros(4,4);
%         [m,Ind]=max(adj,[],2);
%         for j = 1:4
%             red_adj(j,Ind(j)) = m(j);
%         end
        
        list = [4,5,6,7]; check=1; 
        pert =40*.01; prev=pert_init; tmp_pert=pert;

        %% Save the sorted first eigenvector entries
        [V,D]=eig(adj');
        [~,I]=sort(diag(real(D)),'desc');
        [~,II]=sort(abs(V(:,I(1))),'desc');

        pert = std(V(:,I(1)));
%         pert = max(V(:,I(1)))-min(V(:,I(1)));

        if ranked(i)==0
            ranked(i)=pert;%-pert_init;
        end
    else 
        ranked(i)=0.5*.01;
    end
end

if strcmp(type,'plates')
    if strcmp(label,'accurate')
        save('extend_SS_accurate.mat','ranked','nam_save')
    else
        save('std_eigs.mat','ranked','nam_save')
    end
elseif strcmp(type,'inter')
    save('std_eigs_inter_noself.mat','ranked','nam_save')
end
% save('temp_save_OBJ_end.mat','ranked','nam_save')

%%%%%%%%% functions %%%%%%%%%

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
