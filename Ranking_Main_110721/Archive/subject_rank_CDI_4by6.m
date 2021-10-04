clear all
% folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
% folder2 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate\functions';
folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_4by6';
folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Github_CDI';
folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj';
folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
addpath(folder3,folder4,folder5,folder6,folder7)
file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_4by6\';

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =4*6; % number of zones
saved = zeros(num,704);

pert_init=-30;  % Initial Perturbation value
pert=pert_init;
pert_chng = 1;%1000;%
bweight=1;

ranked = zeros(704,1);
round=0;
while min(ranked)==0
    round = round+1;
    f_num = 0;
    pert=pert+pert_chng;
    for i = 1:704
        skip=1;
        file_id = ['subject_',nam_save{i},'.mat'];

        if isfile([file_loc,file_id])
            f_num = f_num + 1;
            load(file_id)
            titlename = ['ID ',nam_save{i}];
            savename = ['subject_',nam_save{i}];
            R = f_num;
        else 
            skip=0;
        end

        if skip == 1
            eig_select = [1,2,3]; 
            adj=adj-diag(diag(adj));
            
            sm_a = sum(adj,1);
            list = find(sm_a==0);
            
            adj(isnan(adj))=0;

            [adj] = NNR_adj_conns_OBJ2(adj,bweight);
            
            L=-adj + diag(sum(adj,2));

            %% convert L into adj (sort of)
            L=L-diag(diag(L));

            list = [9,10,11,12];
            for j = 1 : length(list)
                jj = list(j);
                L(jj,jj)=L(jj,jj)+pert;
            end
            C=zeros(num,num);
            P = -(L + C);

            %% Save the sorted first eigenvector entries
            [V,D]=eig(P');
            [~,I]=sort(diag(real(D)),'desc');
            [~,II]=sort(abs(V(:,I(1))),'desc');
            if mean(real(V(:,I(2)))==real(V(:,I(3)))) == 1
                eig_select = [1,2,4];
            end
            saved(:,i)=II;
        end
       
    end
    
    list = [9,10,11,12];
    for k = 1 : 704
        if ranked(k)==0
            if max(saved(:,k))>0
                if max(find(ismember(saved(:,k),list)==1))>4 % in the top 4
        %         if sum(ismember([1,2,3,4],find(ismember(saved(:,k),list)==1)))~=4   %are list the first 4 entries of saved(:,k)
                    ranked(k)=round;
                end
            else
                ranked(k)=0.5;
            end
        end
    end
end

save('temp_save_GRID.mat')

load('subject_details.mat')
[sets] = set_allocate(subject_details_694,nam_save);
plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets)
title('Grid')