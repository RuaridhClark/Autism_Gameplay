% check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
% folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
% folder2 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate\functions';
folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end_1x2';
folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj_110721';
folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
addpath(folder3,folder4,folder5,folder6,folder7)
file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end_1x2\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =12;    % number of ipad objects (nodes)
saved = zeros(num,704);

pert_init=-30;%-1000;%-10000;%
pert=pert_init;
pert_chng = 1;%1000;%

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
            adj(13:end,:)=[];
            adj(:,13:end)=[];
%             %%%
%             AD = load(['H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_thru2\',file_id],'adj');
%             adj=adj-diag(diag(adj))+ diag(diag(AD.adj)); % Add no. of swipe origins
%             %%%
            titlename = ['ID ',nam_save{i}];
            savename = ['subject_',nam_save{i}];
            R = f_num;

            eig_select = [1,2,3]; 
            
%             adj=adj-diag(diag(adj)); %%% TEMP REMOVAL %%%
            sm_a = sum(adj,1);
%             list = find(sm_a==0);

%             if sum(adj(:))>0
%                 adj = (adj./sum(adj(:)));%.*(100); %%%%%%%%%%%%%%%%%%%%%%% TEMP ADDITION Normalising
%             end
            
            bweight=1;
            [adj] = NNR_adj_conns_OBJ2(adj,bweight);
%             adj(4,:)=[];
%             adj(:,4)=[];
                        
            L=-adj + diag(sum(adj,2));

            %% convert L into adj (sort of)
            L=L-diag(diag(L)); % L=L-diag(ones(length(L),1)*15); % 
            

            list = [4,5,6,7];
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
%             if mean(real(V(:,I(2)))==real(V(:,I(3)))) == 1
%                 eig_select = [1,2,4];
%             end
            saved(:,i)=II;

            
        end
    end
    
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

save('temp_save_OBJ_end.mat')

%%%% Plotting
load('subject_details.mat')
[sets] = set_allocate(subject_details_776,nam_save,saved);
%%
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%
[prcnt,x] = plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets);
title('OBJ')
% title('2 year 6 months – 3 year 8 months')%('3 year 9 months – 4 year 10 months')%('4 year 11 months – 6 year 0 months')

%% Difference All
figure;
plot(x,zeros(1,length(x)));
hold on
ref=prcnt{1};
for i = 2 : length(prcnt)
    diff{i}=prcnt{i}-ref;
    plot(x,diff{i})
    hold on
end
legend('TD','ASD','OND*','ONDE')
% ax = gca;
% ax.XAxisLocation = 'origin';
box off

xlabel('Perturbation magnitude')
ylabel('Difference with respect to TD %')
% title('4 year 11 months – 6 year 0 months')

% %% Difference Gender
% figure;
% plot(x,zeros(1,length(x)));
% hold on
% ref=prcnt{4};i=1;
% % for i = 2 : length(prcnt)
%     diff{i}=prcnt{8}-ref;
%     plot(x,diff{i})
% %     hold on
% % end
% legend('ONDE_F','ONDE_M')
% % ax = gca;
% % ax.XAxisLocation = 'origin';
% box off
% 
% xlabel('Perturbation magnitude')
% ylabel('Difference with respect to ONDE_F %')
% % title('2 year 6 months – 3 year 8 months')