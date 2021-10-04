% clear all
% fileloc='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\IQ_severity';
% tab_sev = readtable([fileloc,'\eCRF.csv']);
% 
% folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end';
% folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
% folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
% folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj';
% folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
% addpath(folder3,folder4,folder5,folder6,folder7)
% file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end\'; % should match zone type
% 
% load('swipes_all704.mat','nam_save')
% 
% %% stack the adjs
% num =16;    % number of ipad objects (nodes)
% saved = zeros(num,704);
% % save_V = zeros(num,704);
% 
% pert_init=-20;%-1000;%-10000;%
% pert=pert_init;
% pert_chng = 1;%1000;%
% bweight=1;
% ranked = zeros(704,1);
% round=0;
% 
% asdname_save = cell(704,1);
% all_severity = zeros(704,1);
% sets = cell(3,1);
% check = 0;
% 
% while min(ranked)==0
%     round = round+1;
%     m = 0;
%     pert=pert+pert_chng;
%     for i = 1:length(nam_save)
%         skip=1;
%         file_id = ['subject_',nam_save{i},'.mat'];
%         
%         jj = find(strcmp(tab_sev.id_study_id,nam_save{i}));
%         score = tab_sev.additional_patient_data__mood(jj);%tab_sev.recording_day_data__distractibility(i);%tab_sev.recording_day_data__interest_in_tablet_games(i);%tab_sev.additional_patient_data__exposure_to_tablets(i);%tab_sev.additional_patient_data__cooperativity(i);%tab_sev.additional_patient_data__arousal(i);%
%         if isfile([file_loc,file_id]) && ~isempty(score)
%             m = m + 1;
%             load(file_id)
%             if check == 0
%                 [score,sets] = data_score(score,m,sets); 
%                 asdname_save{m} = nam_save{jj};
%             end
% %             titlename = ['ID ',nam_save{i}];
% %             savename = ['subject_',nam_save{i}];
%             
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
% 
%             saved(:,m)=II;
% 
%         end
%     end
%     
%     
%     if check == 0 
%     	ranked = zeros(m,1);
%         saved(:,m+1:end) = [];
%         asdname_save(m+1:end) = [];
%     end
%     for k = 1 : m
%         if ranked(k)==0
%             if max(saved(:,k))>0
%                 if max(find(ismember(saved(:,k),list)==1))>4 % in the top 4
%                     ranked(k)=round;
%                 end
%             else
%                 ranked(k)=0.5;
%             end
%         end
%     end
%     check = 1;
% end
% 
% save('temp_save_OBJ.mat')

clear all
for kk = 1 : 6
    load('temp_save_OBJ.mat')

    %% load start
    for a = 1 : 2
        sets = cell(5,1);
        m=0;
        if a == 1    
            subj_grp = 'TD';
        elseif a == 2
            subj_grp = 'ASD';
        end
        for i = 1:height(tab_sev)
            file_id = ['subject_',tab_sev.id_study_id{i},'.mat'];

            [score,type] = score_choice(kk,i,tab_sev);
    %         score = tab_sev.recording_day_data__distractibility(i);%tab_sev.recording_day_data__interest_in_tablet_games(i);%tab_sev.additional_patient_data__exposure_to_tablets(i);%tab_sev.additional_patient_data__cooperativity(i);%tab_sev.additional_patient_data__arousal(i);%tab_sev.additional_patient_data__mood(i);%
            if isfile([file_loc,file_id]) && ~isempty(score) && strcmp(tab_sev.diagnosis_category{i},subj_grp)
                m = m + 1;
                [score,sets] = data_score(score,m,sets); 
            end
        end

        %% Correlation
        % load('diags_incoming.mat')
        % 
        load('subject_details.mat')
        [months] = list_AGE(subject_details_776,asdname_save,saved);
        if size(months,1)<size(months,2)
            months=months';
        end

        exclude = find(ranked==0.5);
        exclude = [exclude;find(months>75)];
        for i = 1 : length(sets)
            rmv=find(ismember(sets{i},exclude')==1);
            if ~isempty(sets{i})
                sets{i}(rmv)=[];
            end
        end

    %     %% plot_cmprsn4 changes ranked internally to match pert
    %     [prcnt,x] = plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets);
    %     legend('1','2','3','4','5','Location','SouthWest')

        %% Boxplot age - Food as an Origin
        if min(ranked)>=0   % change ranked to match pert
            ranked = (ranked-1)*pert_chng+pert_init;
        end
        if a == 1
            sets_save=sets;
        end
    end

    sets = [sets;sets_save];
    
    load('subject_details.mat')
    saved=ones(1,704); saved(list)=zeros(1,length(list));

    iter=0;all_sets=[];grps=[];
    zones=1:17;zones(4)=[];

    for num = 1:length(sets)
    %         h = kstest(diagsA(zones(2),sets{num}))
        iter=iter+1;
        vals = ranked(sets{num})';
        if isempty(vals)
            vals=0;
        end
        all_sets = [all_sets,vals];
        grps = [grps,ones(1,length(vals))*iter];
    end

%     f=figure;
%     boxplot(all_sets,grps,'Notch','on')
%     h = findobj(gca,'Tag','Box');
%     colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880];
%     for m=1:length(h)
%         temp_m = length(h)-m+1;
%         mm=rem(m,6);mm(mm==0)=1;
%         patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
%     end
% 
%     % text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
%     xticklabels({'1','2','3','4','5'});
%     % set(gca,'xticklabel',entries,'fontsize',10)
% 
%     ylabel('Perturbation threshold')
%     box off
% 
%     num_sets=[];
%     for i = 1 : length(sets)
%         num_sets = ['n = ',num2str(length(sets{i}))];
%         text(0.055+(i-1)/5,1,num_sets,'Units','normalized')
%     end
%     title([type,' ',subj_grp])
%     % text(0.05,1,num_sets,'Units','normalized')

    combos=nchoosek([1,2,3,4,5],2);
    save_p = zeros(1,size(combos,1));
    for j = 1 : size(combos,1)
        num = combos(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds,'off');
        save_p(1,j) = pval;
    end
    type
    save_p
%     if min(save_p)>0.05
%         close gcf
%     end
end

% f.Position = [403,340,300,313];

% function [] = stars_line(n_stars,height,strt,nd,age)
%     drp = 2; drp2 = 3.5;
%     hold on
%     if n_stars == 3
%         scatter((age-1)*5+(nd-strt)/2+(strt-0.18),height,'pk','filled')
%         scatter((age-1)*5+(nd-strt)/2+(strt),height,'pk','filled')
%         scatter((age-1)*5+(nd-strt)/2+(strt+0.18),height,'pk','filled')
%     elseif n_stars == 2
%         scatter((age-1)*5+(nd-strt)/2+(strt-0.09),height,'pk','filled')
%         scatter((age-1)*5+(nd-strt)/2+(strt+0.09),height,'pk','filled')
%     elseif n_stars == 1
%         scatter((age-1)*5+(nd-strt)/2+(strt),height,'pk','filled')
%     end
%     hold on
%     plot([(age-1)*5+strt,(age-1)*5+nd],[height-drp,height-drp],'k')
%     plot([(age-1)*5+strt,(age-1)*5+strt],[height-drp,height-drp2],'k')
%     plot([(age-1)*5+nd,(age-1)*5+nd],[height-drp,height-drp2],'k')
% end

%%%%%%%%%%%%%% Function %%%%%%%%%%%%%%%%
function [score,sets] = data_score(score,i,sets)
    if ~isnan(score)
        score(score<2)=2;
        if score>length(sets)
            sets{score}=i;
        else
            sets{score} = [sets{score},i];
        end
    end
end

function [score,type] = score_choice(k,i,tab_sev)
    if k == 1
        score = tab_sev.additional_patient_data__mood(i);
        type = 'mood';
    elseif k ==2
        score = tab_sev.additional_patient_data__arousal(i);
        type = 'arousal';
    elseif k == 3
        score = tab_sev.additional_patient_data__cooperativity(i);
        type = 'cooperativity';
    elseif k == 4
        score = tab_sev.additional_patient_data__exposure_to_tablets(i);
        type = 'exposure to tablet';
    elseif k == 5
        score = tab_sev.recording_day_data__interest_in_tablet_games(i);
        type = 'interest in tablet';
    elseif k == 6
        score = tab_sev.recording_day_data__distractibility(i);
        type = 'distractability';
    end
end