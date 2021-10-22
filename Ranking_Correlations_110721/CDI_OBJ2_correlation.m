% check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
option = 2; % 1 == proportional, 2 == proportion + swipe volume
if option == 1
    load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate_proport.mat')
%     load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')
elseif option == 2
    load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate.mat')
%     load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')
end

% folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
% folder2 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate\functions';
folder3 = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate';
folder4 = 'H:\My Documents\GitHub\Autism_Gameplay\Set_allocate';
folder5 = 'H:\My Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'H:\My Documents\GitHub\Autism_Gameplay\Create_adj_110721';
folder7 = 'H:\My Documents\GitHub\Autism_Gameplay';
addpath(folder3,folder4,folder5,folder6,folder7)
file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =16;    % number of ipad objects (nodes)
saved = zeros(num,704);
save_V = zeros(num,704);

bweight=1;
ranked = zeros(704,1);
round=0;
% while min(ranked)==0
    round = round+1;
    f_num = 0;
    pert=0%
    for i = 1:704
        skip=1;
        file_id = ['subject_',nam_save{i},'.mat'];

        if isfile([file_loc,file_id])
            f_num = f_num + 1;
            load(file_id)
            tadj=adj;
            for xi = 1:num
                if xi>7 && xi<13
                    nxi=xi+4;
                elseif xi>12
                    nxi=xi-5;
                else
                    nxi=xi;
                end
                for yi = 1:num
                    if yi>7 && yi<13
                        nyi=yi+4;
                    elseif yi>12
                        nyi=yi-5;
                    else
                        nyi=yi;
                    end
                    adj(nxi,nyi)=tadj(xi,yi);
                end
            end
            adj=adj(1:num,1:num);
            titlename = ['ID ',nam_save{i}];
            savename = ['subject_',nam_save{i}];
            R = f_num;
%         else 
%             skip=0;
%         end
% 
%         if skip == 1
            eig_select = [1,2,3]; 
    %         adj=zeros(length(adj),length(adj));
            adj=adj-diag(diag(adj));
            sm_a = sum(adj,1);
            list = find(sm_a==0);

            [adj] = NNR_adj_conns_OBJ2(adj,bweight);
            
%             adj=adj'; %% REVERSE ADJ
            
            L=-adj + diag(sum(adj,2));

            %% convert L into adj (sort of)
            L=L-diag(diag(L));

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
            if mean(real(V(:,I(2)))==real(V(:,I(3)))) == 1
                eig_select = [1,2,4];
            end
            saved(:,i)=II;

            save_V(:,i)=abs(V(:,I(1)))-min(abs(V(:,I(1))));
            save_V(:,i)=save_V(:,i)./sum(save_V(:,i));
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
%     ranked= ones(704,1); % exit loop
% end

save('temp_save_OBJ.mat')

load('subject_details.mat')
[sets] = set_allocate(subject_details_776,nam_save,saved);
%%  Remove ADHD
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%
[months] = list_AGE(subject_details_776,nam_save,saved);
% months=months';
% %%% Doesn't remove old subject or put subjects in bins
% num = 1;
% for j = 1:16
%     V_all = save_V(j,:)';
%     figure;
% %     scatter(months(sets{1}),V_all(sets{1}))
%     
%     f=fit(months(sets{num}),V_all(sets{num}),'poly1');
%     plot(f,months(sets{num}),V_all(sets{num}),'x')
% %     [rho,pval] = corr(months(sets{num}),V_all(sets{num}),'Type','Kendall');
% %     text(55,mean(V_all(sets{num})),['p = ',num2str(pval)])
%     [R,p] = corr(months(sets{num}),V_all(sets{num}),'Type','Kendall');
%     text(55,mean(V_all(sets{num})),['R = ',num2str(R),', p = ',num2str(p)])
% %     [R,p] = corrcoef(months(sets{num}),V_all(sets{num}));
% %     text(55,mean(V_all(sets{num}))/2,['R = ',num2str(R(2,1)),', p = ',num2str(p(2,1))])
%     title([num2str(num),' zone ',num2str(j)])
% end

% %%%% Individual trend plot
% for j = 1:16
%     figure;
%     for num = 1 : 4
%         V_all = save_V(j,:)';
% 
%         %%%%%
%         x= months(sets{num});
%         y=V_all(sets{num});
%         [p,S] = polyfit(x,y,1);
%         % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
%         [y_fit,delta] = polyval(p,x,S);
%         % Plot the original data, linear fit, and 95% prediction interval y±2?.
%         [~,Ind]=sort(x,'asc');
%     %     plot(x(Ind),y(Ind),'bx')
%         hold on
%         
%         %%%%%
%         [R,p] = corr(months(sets{num}),V_all(sets{num}),'Type','Kendall');
%         if p<0.05
%             linetype='-';
%         else
%             linetype='--';
%         end
%         plot(x(Ind),y_fit(Ind),linetype)
% %         text(x(Ind(end))+2,y_fit(Ind(end)),['p = ',num2str(p)]) %'R = ',num2str(R),', 
%         
%     end
% %     title([' zone ',num2str(j)])
%     legend('TD','ASD','OND','ONDE','Location','NorthWest')
%     xlabel('Age (months)')
%     ylabel('Adjusted eigenvector entry')
%     axis([30 75 0 0.3])
%     saveas(gcf,[' zone ',num2str(j),'.pdf'])
% end

%% Group trend plot
f=figure;
list = [1,2,3,4,5,6,7,8,9,10,11];%,12,13,14,15,16];
% list=1:16;
xtick_list=[];
for i = 1:length(list)
    j=list(i);
    
    for num = 1 : 4
        V_all = save_V(j,:)';

        %%%%%
        x= months(sets{num})';
        x=(x-30)+((i-1)*80)+10;
        
        y=V_all(sets{num});
        [pf,S] = polyfit(x,y,2);
        % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
        [y_fit,delta] = polyval(pf,x,S);
        % Plot the original data, linear fit, and 95% prediction interval y±2?.
        [~,Ind]=sort(x,'asc');
    %     plot(x(Ind),y(Ind),'bx')
        hold on
        
        %%%%%
        [R,p] = corr(x,y,'Type','Kendall');
        if p<0.001
            linetype='-';
        elseif p<0.01
            linetype='-.';
        elseif p<0.05
            linetype='--';
        else
            linetype=':';
        end
        if num == 1
            clr = [0, 0.4470, 0.7410];
        elseif num == 2
            clr = [0.8500, 0.3250, 0.0980];
        elseif num == 3
            clr = [0.9290, 0.6940, 0.1250];
        elseif num == 4
            clr = [0.4940, 0.1840, 0.5560];
        end
        plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)
    end
    
    legend('TD','ASD','OND*','ONDE','Location','NorthEast')
    xlabel('Age (months)')
    ylabel('Swipe destination popularity (Adjusted eigenvector entry)')
    axis([0 max(x)+10 0 0.4])%axis([0 max(x)+10 0 0.3])

xtick_list = [xtick_list,min(x),(max(x)-min(x))/2+min(x),max(x)];
end
   set(gca,'XTick',xtick_list);

tick_nums=repmat([30,51,72],1,length(list));
xticklabels(tick_nums)

% axis([0 620 0 .3])
grid on
f.Position = [3,335,1360,318];