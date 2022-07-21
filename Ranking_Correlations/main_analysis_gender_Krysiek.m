% 
clear all
num = 16;               % accurate = 16, snap-to = 12
option = 2;             % 1 = n_swipes, 2 = sharing score
destination = 'plates'; % n_swipes for 'plates' or 'food' destinations
gender = 'Female';            % '' or 'Male' or 'Female' or 'compare'
severity = '';        % 'on' or ''

[folder_loc,alt_folder_loc,file_loc,floc] = setup();

[nam_save,~,ranked,~] = load_dataset(option,num,folder_loc);

% load('diags_incoming.mat')
load('subject_details_Krysiek.mat')

[n_swipes,list] = swipe_analysis(num,file_loc,nam_save,destination);

saved=ones(1,386); saved(list)=zeros(1,length(list)); % create list of 1s delete those without adjs
[sets,months] = create_sets_months(subject_details,nam_save,saved,gender,[],severity);

sets = rmv_frm_sets(sets,months);

if ~strcmp(severity,'on')
    sets = rmv_ADHD(sets,subject_details,nam_save,saved);
end

if option == 1
    results = n_swipes;
elseif option == 2
    results = ranked';
end

if strcmp(gender,'compare') % compare male and female
    [all_sets,grps] = create_grps_allsets_gender(results,sets);
    boxplot_gender_cmpr(all_sets,grps,option,destination,num)
    [save_p] = significance_gender(results',sets);
else                        % display all or male/female genders
    for id = 1 : 4          % TD, ASD, OND, ONDE
        if strcmp(severity,'on')
            sev_choice = id-1;
        else
            sev_choice = 0;
        end
        if ~isempty(sets{id})
            plot_results(results',months,sets,id,option,destination,num,gender,sev_choice)
        end
    end
    
    [all_sets,grps] = create_grps_allsets(results,sets);
    
    Plot_boxplots(all_sets,grps,option,destination,num,gender,severity)
    
    [save_p,combos] = significance_check(results,sets);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [folder_loc,alt_folder_loc,file_loc,floc] = setup()
    folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
    alt_folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data';
    file_loc = [folder_loc,'\adjs\adj_krysiek\']; % should match zone type
    floc=[alt_folder_loc,'\IQ_severity'];

    folder1 = [folder_loc,'\Set_allocate'];
%     folder2 = [folder_loc,'\Plots'];
    folder3 = [folder_loc,'\Create_adj'];
    folder4 = [folder_loc,'\adjs\adj_krysiek'];
    folder5 = [folder_loc,'\Data'];
    addpath(folder1,folder3,folder4,folder5)
end

function [nam_save,saved,ranked,list] = load_dataset(option,num,folder_loc)
    if option == 1      % volume
        if num == 16
%             load([folder_loc,'\Ranking_Correlations\Data\OBJ_end_accurate_bi.mat'],'list','nam_save','saved','ranked')
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_accurate_proport_Krysiek.mat'],'list','nam_save','saved')
            ranked=[];
        elseif num == 12
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_end_12zones.mat'],'list','nam_save','saved','ranked')
        end
    elseif option == 2  % proportion
        if num == 16
%             load([folder_loc,'\Ranking_Correlations\Data\OBJ_end_accurate_proport_bi.mat'],'list','nam_save','saved','ranked')
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_accurate_proport_Krysiek.mat'],'list','nam_save','saved','ranked')
        elseif num == 12
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_end_12zones_proport.mat'],'list','nam_save','saved','ranked')
        end
    end
end

function [sets,months] = create_sets_months(subject_details,nam_save,saved,gender,tab_sev,severity)
    [months] = list_AGE(subject_details,nam_save,saved);
    if size(months,1)<size(months,2)
        months=months';
    end
    if strcmp(severity,'on')
        [sets] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,gender);
    elseif strcmp(gender,'Male')
        [tmp_sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
        sets = tmp_sets(5:8);
    elseif strcmp(gender,'Female')
        [tmp_sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
        sets = tmp_sets(1:4);
    elseif strcmp(gender,'compare')
        [sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
    else
        [sets] = set_allocate(subject_details,nam_save,saved);
    end
end


function [sets] = rmv_frm_sets(sets,months)
    exclude = find(months>105); %% 75 months threshold
    % exclude = find(ranked==0.5);
    % exclude = [exclude;489]; % no food-to-plate swipes recorded
    exclude = [exclude;find(months==0)]; %% 45 Months threshold
    for i = 1 : length(sets)
        rmv=find(ismember(sets{i},exclude')==1);
        sets{i}(rmv)=[];
    end
end

function [sets] = rmv_ADHD(sets,subject_details_776,nam_save,saved)
    %%%  Remove ADHD
    load('OND_details.mat','OND_details')
    [sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
    curr_set = sets{3};
    keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
    sets{3}=curr_set(ismember(curr_set,keep));
end

function [all_sets,grps] = create_grps_allsets_gender(results,sets)
    iter=0; all_sets=[]; grps=[];
    order = [1,5,2,6,3,7,4,8];
    
    for i = 1 : length(order)    
        iter = iter +1;
        add_to = [results(sets{order(i)})];
        all_sets = [all_sets,add_to];
        grps = [grps;ones(length(add_to),1)*iter];
        if rem(i,2)==0
            iter=iter+1;
            grps = [grps;iter];
            all_sets = [all_sets,-100];
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

function [n_swipes,list] = swipe_analysis(num,file_loc,nam_save,destination) 
    n_swipes = zeros(1,386);
    f_num = 0;
    
%     diagsA=zeros(12,704);
    list=[];
    for jj = 1:386
%         skip=1;
        file_id = ['subject_',nam_save{jj},'.mat'];
    
        if isfile([file_loc,file_id])
            f_num = f_num + 1;
            load(file_id,'adj')
            if num == 12
                [adj] = adj_snap2zones(adj,num);
            end
            adj=adj(1:num,1:num);
%             titlename = ['ID ',nam_save{jj}];
%             savename = ['subject_',nam_save{jj}];
%             R = f_num;
%         else 
%             skip=0;
        end
    
%         if max(adj(:))==0
%             list=[list,jj];
%         end
    
        if strcmp(destination,'plates')
            n_swipes(jj) = sum(adj(2,[4,5,6,7]));
        elseif strcmp(destination,'food')
            n_swipes(jj) = sum(adj(2,2));
%             n_swipes(jj) = sum(adj(2,:));
%             temp=diag(adj);
% %             sv_2 = temp(2);
%             temp(2)=0;
%             n_swipes(jj) = sum(temp);
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
%     adj=adj-diag(diag(adj));
    
    bweight=1;
    [adj] = NNR_adj_conns_OBJ2(adj,bweight);
end

function [] = plot_results(results,months,sets,id,option,destination,num,gender,sev_choice)
    figure;
    
    x=months(sets{id});
    y=results(sets{id});
    [pf,S] = polyfit(x,y,2);
    % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
    [y_fit,~] = polyval(pf,x,S);
    % Plot the original data, linear fit, and 95% prediction interval y±2?.
    [~,Ind]=sort(x,'asc');
    [~,p] = corr(x,y,'Type','Spearman');
    
    if p<0.01
        linetype='-';
    elseif p<0.05
        linetype='--';
    else
        linetype=':';
    end

    if sev_choice>0
        if sev_choice==1
            clr = [.93,.69,.13];
        elseif sev_choice==2
            clr = [.85,.33,.1];
        elseif sev_choice==3
            clr=[.64,.08,.18];
        end
    else
        if id == 1
            clr = [0, 0.4470, 0.7410];
        elseif id == 2
            clr = [0.8500, 0.3250, 0.0980];
        elseif id == 3
            clr = [0.9290, 0.6940, 0.1250];
        elseif id == 4
            clr = [0.4940, 0.1840, 0.5560];
        end
    end

%     tempclr=[.5 .5 .5];
    scatplot=scatter(x(Ind),y(Ind),55,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatplot.MarkerFaceAlpha = .3;
    hold on
    if p<0.05
        plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)
    end

    if p<0.001
        text(0.05,.95,['p = ',sprintf('%1.1e', p)],'Units','normalized','fontsize', 11)
    elseif p<0.01
        text(0.05,.95,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    elseif p<10%0.1
        text(0.05,.95,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    end
    
    if option == 2
        axis([29 75 -1 .25])
%         axis([min(x) max(x) min(y) max(y)])
    elseif option == 1
        axis([29 75 0 110])
%        axis([min(x) max(x) min(y) max(y)])
    end

    xlabel('Age (months)','fontsize', 11)
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
    legend('subject','2nd order fit','Location','SouthEast')
    
    if id == 1
        titletxt = 'TD';
    elseif sev_choice == 1
        titletxt = 'Level 1';
    elseif sev_choice == 2
        titletxt = 'Level 2';
    elseif sev_choice == 3
        titletxt = 'Level 3';
    elseif id == 2
        titletxt = 'ASD';
    elseif id == 3
        titletxt = 'OND';
    elseif id == 4
        titletxt = 'ONDE';
    end

    if strcmp(gender,'')
        title(titletxt)
    else
        title([titletxt,' - ',gender])
    end

    f=gcf;
    f.Position = [403,340,330,313];
end

function [] = stars_line(n_stars,height,strt,nd,drp)
    drp2 = 3.5*drp/2;
    hold on
    if n_stars == 3
        scatter((nd-strt)/2+(strt-0.18),height,'pk','filled')
        scatter((nd-strt)/2+(strt),height,'pk','filled')
        scatter((nd-strt)/2+(strt+0.18),height,'pk','filled')
    elseif n_stars == 2
        scatter((nd-strt)/2+(strt-0.09),height,'pk','filled')
        scatter((nd-strt)/2+(strt+0.09),height,'pk','filled')
    elseif n_stars == 1
        scatter((nd-strt)/2+(strt),height,'pk','filled')
    end
    hold on
    plot([strt,nd],[height-drp,height-drp],'k')
    plot([strt,strt],[height-drp,height-drp2],'k')
    plot([nd,nd],[height-drp,height-drp2],'k')
end

function [] = Plot_boxplots(all_sets,grps,option,destination,num,gender,severity)
    f=figure;
    b=boxplot(all_sets,grps,'Notch','on','Color',[.5,.5,.25]);
    set(b,'LineWidth',1.5)
%     h = findobj(gca,'Tag','Box');
    
    %% Add scatter points
    hold on
    [~,~,ic]=unique(grps,'stable');
    
    if strcmp(severity,'on')
        colors = [0, 0.4470, 0.7410;.93,.69,.13; .85,.33,.1; .64,.08,.18];
    else
        colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
    end
        
    for i = 1:4
        Ind=find(ic==i);
        scatter(ic(Ind),all_sets(Ind),[],colors(i,:),'filled','MarkerFaceAlpha',0.15,'jitter','on','jitterAmount',0.15);
    end
    
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

    if option == 2
        axis([0.5 4.5 -1 .39])
    elseif option == 1
        axis([0.5 4.5 -3 150])
    end

    %% Plot significance stars
    if option == 2
    %     heights = [.295,.345,.41,.295];
        heights = [.41,.345,.45,.31];
        drp = .0225;
    elseif option == 1
        heights = [140,128,147,122];
        drp = 2;
    end
    
    n_stars = 3; 
%     stars_line(3,heights(1),1,2,drp) % 3 stars,h,1,2,1
%     stars_line(3,heights(2),1,3,drp) % 3 stars,h,1,2,1
%     stars_line(3,heights(3),2,4,drp) % 2 stars,h,1,3,2
%     stars_line(3,heights(4),3,4,drp) % 2 stars,h,3,4,1
       
    % legend('TD','ASD','OND*','ONDE','Orientation','horizontal')
    
    f.Position = [403,340,300,313];

    title(gender)
end

function [save_p,combos] = significance_check(results,sets)
    % pairwise significance check
    combos=nchoosek([1,2,3,4],2);
    save_p = zeros(1,size(combos,1));
    for j = 1 : size(combos,1)
        num = combos(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval = kruskalwallis([results(sets{num(1)})';results(sets{num(2)})'],len_rankeds,'off');
        save_p(1,j) = pval;
    end
end

function [save_p] = significance_gender(results,sets)
    save_p = zeros(1,4);
    numopt = [1,5;2,6;3,7;4,8];
    for j = 1 : 4
        num=numopt(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval = kruskalwallis([results(sets{num(1)});results(sets{num(2)})],len_rankeds,'off');
        save_p(1,j) = pval;
        if pval > .001
            text(0.05+(j-1)*.95/4,.95,['p = ',sprintf('%1.4f', pval)],'Units','normalized','fontsize', 11)
        else
            text(0.05+(j-1)*.95/4,.95,['p = ',sprintf('%1.1e', pval)],'Units','normalized','fontsize', 11)
        end
    end
end

function [ ] = boxplot_gender_cmpr(all_sets,grps,option,destination,num)

    f=figure;
    b=boxplot(all_sets,grps,'Notch','on','Color',[.5,.5,.25]);
    set(b,'LineWidth',1.5)
%     h = findobj(gca,'Tag','Box');

     %% Add scatter points
    hold on
    [~,~,ic]=unique(grps,'stable');
    
    colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
    for j = 1:12
        i=ceil(j/3);
        Ind=find(ic==j);
        scatter(ic(Ind),all_sets(Ind),[],colors(i,:),'filled','MarkerFaceAlpha',0.15,'jitter','on','jitterAmount',0.15);
    end
    
    hAx=gca; 
    xt=hAx.XTick;                   % retrieve ticks
    xt(3:3:end)=[];
    hAx.XTick=xt;                   % clear some
    
    entries ={'F','M','F','M','F','M','F','M'};
    xticklabels(entries)

    if option == 1
        axis([0 12 -5 120])
    elseif option == 2
        axis([0 12 -.7 .4])
    end
    
    %% Plot significance stars
    % height=30; n_stars = 2; drp = 1.75;
    % stars_line(n_stars,height,1,2,drp) % 3 stars,h,1,2,1
    % height=30; n_stars = 1; drp = 1.75;
    % stars_line(n_stars,height,10,11,drp) % 3 stars,h,1,2,1
    box off
    ylabel('Perturbation threshold (proportion + swipe volume)')
    f.Position = [403,340,574,313];
    xlabel('Gender')

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
end