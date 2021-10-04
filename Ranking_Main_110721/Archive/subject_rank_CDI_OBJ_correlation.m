% % check NNR_adj_conns_OBJ2 and pert changes for velocity case
% clear all
% % folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
% % folder2 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate\functions';
% folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj';
% folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Github_CDI';
% folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
% folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj';
% folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
% addpath(folder3,folder4,folder5,folder6,folder7)
% file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj\'; % should match zone type
% 
% load('swipes_all704.mat','nam_save')
% 
% %% stack the adjs
% num =16;    % number of ipad objects (nodes)
% saved = zeros(num,704);
% 
% pert_init=-40;%-1000;%-10000;%
% pert=pert_init;
% pert_chng = 1;%1000;%
% bweight=1;
% ranked = zeros(704,1);
% round=0;
% while min(ranked)==0
%     round = round+1;
%     f_num = 0;
%     pert=pert+pert_chng;
%     for i = 1:704
%         skip=1;
%         file_id = ['subject_',nam_save{i},'.mat'];
% 
%         if isfile([file_loc,file_id])
%             f_num = f_num + 1;
%             load(file_id)
%             titlename = ['ID ',nam_save{i}];
%             savename = ['subject_',nam_save{i}];
%             R = f_num;
%         else 
%             skip=0; %% MAY NEED TO REMOVE THESE FROM DATASET? Currently filter them out in set_allocate
%         end
% 
%         if skip == 1
%             eig_select = [1,2,3]; 
%     %         adj=zeros(length(adj),length(adj));
%             adj=adj-diag(diag(adj));
%             sm_a = sum(adj,1);
%             list = find(sm_a==0);
% 
%             [adj] = NNR_adj_conns_OBJ2(adj,bweight);
%             adj(4,:)=[];
%             adj(:,4)=[];
%             
%             L=-adj + diag(sum(adj,2));
% 
%             %% convert L into adj (sort of)
%             L=L-diag(diag(L));
% 
%             list = [5,6,7,8];
%             list(list>3)=list(list>3)-1; % minus 1 from ind
%             for j = 1 : length(list)
%                 jj = list(j);
%                 L(jj,jj)=L(jj,jj)+pert;
%             end
% 
%             C=zeros(num,num);
%             P = -(L + C);
% 
%             %% Save the sorted first eigenvector entries
%             [V,D]=eig(P');
%             [~,I]=sort(diag(real(D)),'desc');
%             [~,II]=sort(abs(V(:,I(1))),'desc');
% %             if mean(real(V(:,I(2)))==real(V(:,I(3)))) == 1
% %                 eig_select = [1,2,4];
% %             end
%             saved(:,i)=II;
% 
%             
%         end
%     end
%     
%     for k = 1 : 704
%         if ranked(k)==0
%             if max(saved(:,k))>0
%                 if max(find(ismember(saved(:,k),list)==1))>4 % in the top 4
%         %         if sum(ismember([1,2,3,4],find(ismember(saved(:,k),list)==1)))~=4   %are list the first 4 entries of saved(:,k)
%                     ranked(k)=round;
%                 end
%             else
%                 ranked(k)=0.5;
%             end
%         end
%     end
% end
% 
% save('temp_save_OBJ.mat')

load('subject_details.mat')
[sets] = set_allocate(subject_details_776,nam_save,saved);

num = 1;

%%%%
% set=cell(1,3);
% for age = 1 : 3
%     [sets] = set_allocate_AGE(subject_details_776,nam_save,saved,age);
%     set{age}=sets;
% end
% sets=cell(1,3);
% for i = 1 :3
%     sets{i}=set{i}{num};
% end

%%%%%
% [prcnt,x] = plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets);
% title('OBJ')
% % title('2 year 6 months – 3 year 8 months')%('3 year 9 months – 4 year 10 months')%('4 year 11 months – 6 year 0 months')

%%%%%%%%%
months=zeros(704,1);
for i = 1 : 704
    months(i,1) = subject_details_776{i,2}*12+subject_details_776{i,3};
end
months(months>=59 & months<73)=65;
months(months>=45 & months<59)=51;
months(months<45)=37;

% months(months>=63 & months<73)=68;
% months(months>=52 & months<63)=57;
% months(months>=41 & months<52)=46;
% months(months<41)=35;

% months(months>=63 & months<73)=68;
% months(months>=52 & months<63)=57;
% months(months>=41 & months<52)=46;
% months(months<39)=35;

months(months>=63 & months<73)=68;
months(months>=52 & months<63)=57;
months(months>=41 & months<52)=46;
months(months<41)=35;

temp=sets{1};
temp(temp==120)=[];
sets{1}=temp;

figure;
f=fit(months(sets{num}),ranked(sets{num}),'poly1');
plot(f,months(sets{num}),ranked(sets{num}),'o')
[rho,pval] = corr(months(sets{num}),ranked(sets{num}),'Type','Kendall');
text(55,mean(ranked(sets{num})),num2str(pval))
title([num2str(num)])