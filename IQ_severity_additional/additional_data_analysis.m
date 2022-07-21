clear all

subj_grp = 'ASD';
sev_choice = 2;        % set as empty string or 1,2,3
num=16;                 % accurate = 16, snap-to = 12
option = 1;             % 1 = n_swipes, 2 = sharing score
destination = 'plates';   % n_swipes for 'plates' or 'food' destinations
measure = 'additional'; % measure 'additional'
gender = 'Male';

[folder_loc,alt_folder_loc,file_loc,floc,n_meas] = setup(measure);

for kk = 1 : n_meas      % different scores (e.g. mood)
    [nam_save,saved,ranked,list] = load_dataset(option,num,folder_loc);
    tab_sev = readtable([floc,'\eCRF.csv']);

    [n_swipes,score,sets,asdname_save,titlename,savename,map,type] = analyse_subjects(kk,tab_sev,file_loc,subj_grp,nam_save,num,destination,measure,sev_choice,gender);

    [months] = load_months(folder_loc,asdname_save,saved);

    [ranked,n_swipes,map,months] = clean_results(ranked,n_swipes,saved,map,months);

    months = check_months_format(months);

    exclude = find(months>75);
    [sets] = create_sets(sets,map,exclude);

    load('subject_details.mat')

    if option == 1
        results = n_swipes;
    elseif option == 2
        results = ranked';
    end

    [all_sets,grps] = create_grps_allsets(results,sets);	
		 
    plot_results(subj_grp,sev_choice,all_sets,sets,grps,option,destination,num,gender)

    [pval,tbl,stats] = test_significance(results,sets);

    check_significance(type,pval)
    
end

%%%%%%%%%%%%%% Function %%%%%%%%%%%%%%%%
function [folder_loc,alt_folder_loc,file_loc,floc,n_meas] = setup(measure)
    folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
    alt_folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data';
    file_loc = [folder_loc,'\adjs_110721\adj_obj_end_accurate\']; % should match zone type
    floc=[alt_folder_loc,'\IQ_severity'];

    folder1=[folder_loc,'\Set_allocate'];
    folder2 = folder_loc;
    addpath(folder1,folder2)

    if strcmp(measure,'additional')
        n_meas = 6;
    elseif strcmp(measure,'IQ')
        n_meas = 12;
    end
end

function [nam_save,saved,ranked,list] = load_dataset(option,num,folder_loc)
    if option == 1      % volume
        if num == 16
            load([folder_loc,'\Ranking_Correlations_110721\Data\OBJ_end_accurate_bi.mat'],'list','nam_save','saved','ranked')
        elseif num == 12
            load([folder_loc,'\Ranking_Correlations_110721\Data\OBJ_end_12zones.mat'],'list','nam_save','saved','ranked')
        end
    elseif option == 2  % proportion
        if num == 16
            load([folder_loc,'\Ranking_Correlations_110721\Data\OBJ_end_accurate_proport_bi.mat'],'list','nam_save','saved','ranked')
        elseif num == 12
            load([folder_loc,'\Ranking_Correlations_110721\Data\OBJ_end_12zones_proport.mat'],'list','nam_save','saved','ranked')
        end
    end
end

function [score,sets] = reshape_score(score,i,sets)
    if ~isnan(score)
%         score(score<2)=2;
        if score>length(sets) % create first entry
            sets{score}=i;
        else
            sets{score} = [sets{score},i];
        end
    end
end

function [score,type] = score_choice(k,i,tab_sev,measure)
    if strcmp(measure,'additional')
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
end


function [adj] = adj_snap2zones(adj,num)
    % rewire adjacency from 16 to 12 zones
    for it = 1:num
        adj(it,4:7)=adj(it,4:7)+adj(it,13:16);    % reconnect to 4-7
        adj(4:7,it)=adj(4:7,it)+adj(13:16,it);    % reconnect to 13-16
        adj(it,13:16)=zeros(1,4);                     % remove non-food connections
        adj(13:16,it)=zeros(4,1); 
    end
    adj = adj(1:12,1:12);
    %% remove zn 4-7 incoming except from 2
    allow=[2,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            adj(it,4:7)=zeros(1,4);                     % remove non-food connections
        end
    end
    adj=adj-diag(diag(adj));
    
    bweight=1;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);
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

function [n_swipes,score,sets,asdname_save,titlename,savename,map,type] = analyse_subjects(kk,tab_sev,file_loc,subj_grp,nam_save,num,desination,measure,sev_choice,gender)
    sets = cell(5,1);sev_sets = cell(3,1);
    n_swipes = zeros(1,height(tab_sev));
    asdname_save = cell(1,height(tab_sev));
    map = [];
    for i = 1:height(tab_sev)
        file_id = ['subject_',tab_sev.id_study_id{i},'.mat'];
        [score,type] = score_choice(kk,i,tab_sev,measure);
        severity = tab_sev.clinical_diagnosis__asd_severity_level{i};
        if isfile([file_loc,file_id]) && ~isempty(score) && strcmp(tab_sev.diagnosis_category{i},subj_grp)
            if strcmp(tab_sev.patient_data__sex{i},gender) || strcmp(gender,'')  % gender either 'Male','Female',or ''
                if strcmp(subj_grp,'ASD') && isnumeric(sev_choice)
                    if ~strcmp(severity,'')
                        [sev_num,~] = severity_score(severity,1,sev_sets); 
                        if sev_num == sev_choice
                            [n_swipes,score,sets,asdname_save,titlename,savename,map] = swipe_analysis(i,map,score,sets,n_swipes,tab_sev,nam_save,num,file_loc,file_id,asdname_save,desination);
                        end
                    end
                else
                    [n_swipes,score,sets,asdname_save,titlename,savename,map] = swipe_analysis(i,map,score,sets,n_swipes,tab_sev,nam_save,num,file_loc,file_id,asdname_save,desination);
                end
            end
        end
    end
end

function [n_swipes,score,sets,asdname_save,titlename,savename,map] = swipe_analysis(i,map,score,sets,n_swipes,tab_sev,nam_save,num,file_loc,file_id,asdname_save,destination)
    map(i)=find(strcmp(nam_save,tab_sev.id_study_id{i}));
    load([file_loc,file_id],'adj')
    if num == 12
        [adj] = adj_snap2zones(adj,num);
    end
    adj = adj(1:num,1:num);
    if strcmp(destination,'plates')
        n_swipes(map(i)) = sum(adj(2,[4,5,6,7]));
    elseif strcmp(destination,'food')
%         tmp = sum(diag(adj));
%         tmp(2)=0;
        n_swipes(map(i)) = um(adj(2,2));
    end
    [score,sets] = reshape_score(score,map(i),sets);
    asdname_save{i} = tab_sev.id_study_id{i};
    titlename = ['ID ',tab_sev.id_study_id{i}];
    savename = ['subject_',tab_sev.id_study_id{i}];
end

function [months] = load_months(folder_loc,asdname_save,saved) 
    load([folder_loc,'\subject_details.mat'])
    [tmp_months] = list_AGE(subject_details_776,asdname_save,saved);
    months=zeros(1,694);
    months(1:length(tmp_months))=tmp_months;
end

function [ranked,n_swipes,map,months] = clean_results(ranked,n_swipes,saved,map,months)
    % map, n_swipes have the same index, ranked 
    rmv = find(map==0);
    months(rmv)=[];
    map(rmv)=[];
%     n_swipes(rmv)=[];
% 
%     saved=saved(map);
% 
%     rmv = find(sum(saved,1)==0);
%     map(rmv)=[];
%     n_swipes(rmv)=[];
    n_swipes = n_swipes(map);
    ranked=ranked(map);
    ranked(ranked==0.5)=[];
end

function months = check_months_format(months)
    if size(months,1)<size(months,2)
        months=months';
    end
end

% function [ranked] = convert_ranked(ranked,pert_chng,pert_init);
%     if min(ranked)>=0   % change ranked to match pert
%         ranked = (ranked-1)*pert_chng+pert_init;
%     end
% end

function [sets] = create_sets(sets,map,exclude)
    for i = 1 : length(sets)
       	for jj = 1 : length(sets{i})
            sets{i}(jj)=find(map==sets{i}(jj));
        end
        rmv=find(ismember(sets{i},exclude')==1);
        if ~isempty(sets{i})
            sets{i}(rmv)=[];
        end
    end
end

function [all_sets,grps] = create_grps_allsets(results,sets)
    iter=0;all_sets=[];grps=[];
    for num = 1:length(sets)
    %         h = kstest(diagsA(zones(2),sets{num}))
        iter=iter+1;
        vals = results(sets{num});
        if isempty(vals)
            vals=0;
        end
        all_sets = [all_sets,vals];
        grps = [grps,ones(1,length(vals))*iter];
    end	
end

function [] = plot_results(subj_grp,sev_choice,all_sets,sets,grps,option,destination,num,gender)

    f=figure;
    b=boxplot(all_sets,grps,'Notch','on','Color',[.5,.5,.25]);
    set(b,'LineWidth',1.5)
    h = findobj(gca,'Tag','Box');

    %% Add scatter points
    hold on
    [C,~,ic]=unique([grps],'stable');

    if isnumeric(sev_choice)
        if sev_choice==1
            clr = [.93,.69,.13];
        elseif sev_choice==2
            clr = [.85,.33,.1];
        elseif sev_choice==3
            clr=[.64,.08,.18];
        end
    else
        if strcmp(subj_grp,'TD')
            clr = [0, 0.4470, 0.7410];
        elseif strcmp(subj_grp,'ASD')
            clr = [0.8500, 0.3250, 0.0980];
        end
    end

    scatter(ic,all_sets,[],clr,'filled','MarkerFaceAlpha',0.5,'jitter','on','jitterAmount',0.15);

    xticklabels({'1','2','3','4','5'});

    if option == 1
        if strcmp(destination,'plates')
            if num == 16
                ylabel('No. of swipes - food to plates','fontsize',14)
            elseif num == 12
                ylabel('No. of swipes - food to snap-to-target','fontsize',14)
            end
        elseif strcmp(destination,'food')
            ylabel('No. of swipes - zone 2 only','fontsize',14)
        end
    elseif option == 2
        if num == 16
            ylabel('Sharing score (plates)','fontsize',14)
        elseif num == 12
            ylabel('Sharing score (snap-to-target)','fontsize',14)
        end
    end

    box off

    for i = 1 : length(sets)
        num_sets = ['n = ',num2str(length(sets{i}))];
        text(0.055+(i-1)/5,.99,num_sets,'Units','normalized','fontsize',14)
    end

    if strcmp(gender,'Male')
        title('Male')
    elseif strcmp(gender,'Female')
        title('Female')
    end
end

function [pval,tbl,stats] = test_significance(results,sets)
    len_rankeds = [ones(length(sets{1}),1);2*ones(length(sets{2}),1);...
                    3*ones(length(sets{3}),1);4*ones(length(sets{4}),1);...
                    5*ones(length(sets{5}),1)];
    [pval,tbl,stats] = anova1([results(sets{1})';results(sets{2})';results(sets{3})';results(sets{4})';results(sets{5})'],len_rankeds,'off');
end

function [] = check_significance(type,pval)
    if min(pval)>0.05
        close gcf
    else
        type
        pval
        if pval<0.001
            text(0.02,1.05,['p_{ANOVA} = ',sprintf('%1.1e', pval)],'Units','normalized','fontsize',14)
        elseif pval<0.01
            text(0.02,1.05,['p_{ANOVA} = ',sprintf('%1.4f', pval)],'Units','normalized','fontsize',14)
        elseif pval<10%0.1
            text(0.02,1.05,['p_{ANOVA} = ',sprintf('%1.4f', pval)],'Units','normalized','fontsize',14)
        end
            xlabel([type,' score'],'fontsize',14)
%             title(subj_grp)
    end
end