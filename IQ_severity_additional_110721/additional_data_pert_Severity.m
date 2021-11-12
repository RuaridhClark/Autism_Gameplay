clear all
folder1='H:\My Documents\GitHub\Autism_Gameplay\Set_allocate';
folder2 = 'H:\My Documents\GitHub\Autism_Gameplay\';
num=16;
addpath(folder1,folder2)
option = 2;
for kk = 1 : 6
    if option == 1
%         load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')
        load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate.mat')
    elseif option == 2
%         load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')
        load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate_proport_bi.mat')
    end
    file_loc = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate\'; % should match zone type
    
    fileloc='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\IQ_severity';
    tab_sev = readtable([fileloc,'\eCRF.csv']);
    
    %% load start
    sets = cell(5,1);
    sev_sets = cell(3,1);
    map = [];
    subj_grp = 'ASD';
    for i = 1:height(tab_sev)
        file_id = ['subject_',tab_sev.id_study_id{i},'.mat'];
        [score,type] = score_choice(kk,i,tab_sev);
        
        severity = tab_sev.clinical_diagnosis__asd_severity_level{i};
        if isfile([file_loc,file_id]) && ~isempty(score) && strcmp(tab_sev.diagnosis_category{i},subj_grp)
            if ~strcmp(severity,'')
                [sev_num,~] = severity_score(severity,1,sev_sets); 
                if sev_num == 3
                    map(i)=find(strcmp(nam_save,tab_sev.id_study_id{i}));
                    load([file_loc,file_id])
                    adj = adj(1:num,1:num);
                    n_swipes(i) = sum(adj(2,[4,5,6,7]));
                    [score,sets] = data_score(score,map(i),sets);

                    asdname_save{i} = tab_sev.id_study_id{i};
                    titlename = ['ID ',tab_sev.id_study_id{i}];
                    savename = ['subject_',tab_sev.id_study_id{i}];
                end
            end
        end
    end

    load('subject_details.mat')
    [tmp_months] = list_AGE(subject_details_776,asdname_save,saved);
    months=zeros(1,694);
    months(1:length(tmp_months))=tmp_months;
%     n_swipes=[n_swipes,0];
    tmp_swipes = n_swipes;
	n_swipes = zeros(1,694);
    n_swipes(1:length(tmp_swipes))=tmp_swipes;

    rmv = find(map==0);
    months(rmv)=[];
    map(rmv)=[];
    n_swipes(rmv)=[];

    saved=saved(map);

    rmv = find(sum(saved,1)==0);
    map(rmv)=[];
    n_swipes(rmv)=[];
    ranked=ranked(map);
    ranked(ranked==0.5)=[];

    if size(months,1)<size(months,2)
        months=months';
    end

    exclude = find(months>75);
    for i = 1 : length(sets)
       	for jj = 1 : length(sets{i})
            sets{i}(jj)=find(map==sets{i}(jj));
        end
        rmv=find(ismember(sets{i},exclude')==1);
        if ~isempty(sets{i})
            sets{i}(rmv)=[];
        end
    end

    %% Boxplot age - Food as an Origin
    if min(ranked)>=0   % change ranked to match pert
        ranked = (ranked-1)*pert_chng+pert_init;
    end

    load('subject_details.mat')
    saved=ones(1,704); saved(list)=zeros(1,length(list));

    iter=0;all_sets=[];grps=[];

    for num = 1:length(sets)
    %         h = kstest(diagsA(zones(2),sets{num}))
        iter=iter+1;
        vals = ranked(sets{num})';%n_swipes(sets{num});%
        if isempty(vals)
            vals=0;
        end
        all_sets = [all_sets,vals];
        grps = [grps,ones(1,length(vals))*iter];
    end

    f=figure;
    boxplot(all_sets,grps,'Notch','on')
    h = findobj(gca,'Tag','Box');
%     colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880];
%     colors = [.93,.69,.13;.85,.33,.1;0,.45,.74;.64,.08,.18;.3,.75,.93];
    colors = [.64,.08,.18;.85,.33,.1;.93,.69,.13;.3,.75,.93;0,.45,.74];
    for m=1:length(h)
        temp_m = length(h)-m+1;
        mm=rem(m,6);mm(mm==0)=1;
        patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
    end

    % text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
    xticklabels({'1','2','3','4','5'});
    % set(gca,'xticklabel',entries,'fontsize',10)

    if option == 1
        ylabel('Perturbation threshold (proportion + swipe volume)')
    elseif option == 2
        ylabel('Perturbation threshold (proportion)')
    end
    box off

    num_sets=[];
    for i = 1 : length(sets)
        num_sets = ['n = ',num2str(length(sets{i}))];
        text(0.055+(i-1)/5,.98,num_sets,'Units','normalized')
    end
%     title([type,' ',subj_grp])
    % text(0.05,1,num_sets,'Units','normalized')

    len_rankeds = [ones(length(sets{1}),1);2*ones(length(sets{2}),1);...
                    3*ones(length(sets{3}),1);4*ones(length(sets{4}),1);...
                    5*ones(length(sets{5}),1)];
    [pval,tbl,stats] = anova1([ranked(sets{1});ranked(sets{2});ranked(sets{3});ranked(sets{4});ranked(sets{5})],len_rankeds,'off');
    save_p = pval;
    
    title(num2str(sev_choice))
    
    if min(save_p)>0.05
        close gcf
    else
        type
        save_p
        if pval<0.001
            text(0.02,1.03,['p_{ANOVA} = ',sprintf('%1.1e', pval)],'Units','normalized')
        elseif pval<0.01
            text(0.02,1.03,['p_{ANOVA} = ',sprintf('%1.4f', pval)],'Units','normalized')
        elseif pval<10%0.1
            text(0.02,1.03,['p_{ANOVA} = ',sprintf('%1.4f', pval)],'Units','normalized')
        end
            xlabel([type,' score'])
%             title(subj_grp)
    end
    
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
%         score(score<2)=2;
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

function [sev_num,sets] = severity_score(severity,m,sets)
    if strcmp(severity,'1. Level 1 "Requiring support"')
        sev_num = 1;
        sets{1} = [sets{1},m];
    elseif strcmp(severity,' 2. Level 2 "Requiring substantial support"')
        sev_num = 2;
        sets{2} = [sets{2},m];
    elseif strcmp(severity,'3. Level 3 "Requiring very substantial support"')
        sev_num = 3;
        sets{3} = [sets{3},m];
    end
end