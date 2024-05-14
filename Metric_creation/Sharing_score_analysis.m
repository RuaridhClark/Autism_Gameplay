%% Sharing Score Analysis
% Identify the sharing score for each user by perturbing the matrix and
% monitoring the eigenvectors of the plate zones.

clear all
folder1 = '..\Create_adj';
folder2 = '..\';
addpath(folder1,folder2)
file_loc = '..\adjs\adj_zones\';
option = "pretrial";

if option == "trial"
    load('swipes_trial.mat','nam_save')
elseif option == "pretrial"
    load('swipes_pretrial.mat','nam_save')
end

%%
num =16;        % number of ipad zones (nodes)
pert_init=-.80; % initial perturbation

saved = zeros(num,length(nam_save));
ranked = zeros(length(nam_save),1);

type = 'plates';    % 'plates' for direct sharing score, 'inter' for indirect sharing score 
bweight = 0.001;

for i = 1:length(nam_save)
    file_id = ['subject_',nam_save{i},'.mat'];

    if isfile([file_loc,file_id])
        load([file_loc,file_id])
        titlename = ['ID ',nam_save{i}];
        savename = ['subject_',nam_save{i}];

        if num == 16 
            if strcmp(type,'plates')
                L = adj2L_snap2zones_foodloc(adj,num,bweight);  % direct food delivery
            elseif strcmp(type,'inter')
                L = adj2L_interplate(adj,num,bweight);    % inter-plate analysis
            end
        end
        
        list = [4,5,6,7]; check=1; 
        pert =40*.01; prev=pert_init; tmp_pert=pert;
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
                [pert,prev,check] = pert_bisect(pert,prev,tmp_pert,check);
            elseif check == 10 || check == -10
                [pert] = pert_iter(pert,check);
            end
        end
        if ranked(i)==0
            ranked(i)=pert;
        end
    else 
        ranked(i)=NaN;
    end
end

[w_string] = save_title(option,bweight);

if strcmp(type,'plates')
    save(append('..\Results_comparison\Data\SS_ext',w_string,'.mat'),'ranked','nam_save')
elseif strcmp(type,'inter')
    save(append('..\Results_comparison\Data\SS_ext_inter',w_string,'.mat'),'ranked','nam_save')
end

%%%%%%%%% functions %%%%%%%%%
function [L] = adj2L_snap2zones_foodloc(adj,num,bweight)
    
    adj(2,4:7)=adj(2,4:7)+adj(2,13:16);    % Collect all food delivery swipes
    adj(2,13:16)=zeros(1,4);               % remove re-connected connections

    %% remove zn 4-7 incoming except from 2
    allow=2;
    for it = 1:num
        if ~ismember(it,allow)
            if num == 16
                adjust = adj(it,4:7);
                adj(it,13:16)=adj(it,13:16)+adjust;    % reconnect to 13-16
            end
            list = 4:7;
            listedit =list(~ismember(4:7,it));
            adj(it,listedit)=zeros(1,length(listedit));                     % remove non-food connections
        end
    end

    adj = adj-diag(diag(adj));

    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));% Normalising
    end

    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L = -adj;
end

function [L] = adj2L_interplate(adj,num,bweight)
    %% All food delivery version

    %% remove zn 4-7 incoming except from 2,4,5,6,7
    allow=[2,4,5,6,7];
    % Snap-to-plate
    for i = 1:length(allow)
        adj(allow(i),4:7)=adj(allow(i),4:7)+adj(allow(i),13:16);    % reconnect 2 to 13-16
        dest = 13:16;
        vals = dest(~ismember(dest,allow(i)+9));
        adj(allow(i),vals)=zeros(1,length(vals));               % remove re-connected connections
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
    
    adj = adj-diag(diag(adj));
    if sum(adj(:))>0
        adj = (adj./sum(adj(:)));% Normalising
    end

    [adj] = NNR_adj_conns_OBJ2(adj,bweight);

    L=-adj;

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

function [pert,prev,check] = pert_bisect(pert,prev,tmp_pert,check)
    dif = pert-prev;
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
        pert=tmp_pert;
    end
end

function [pert] = pert_iter(pert,check)
    if check == 10
        iter = 0.001;
        pert = pert + iter;
    elseif check == -10
        iter = -0.001;
        pert = pert + iter;
    end
end

function [adj] = NNR_adj_conns_OBJ2(adj,bweight)
    %% create all to all connections
    Complete_graph = ones(length(adj),length(adj))*bweight; 
    Complete_graph(isnan(Complete_graph))=0;
    adj=adj+Complete_graph;
end

function [w_string] = save_title(option,bweight)
    % Convert the value to string
    w_string = sprintf('_%.3f', bweight);
    % Replace the decimal point with underscore
    w_string = strrep(w_string, '.', '_');

    if strcmp(option,'pretrial')
        w_string = append('_pretrial',w_string);
    end
end
