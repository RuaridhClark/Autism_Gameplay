% 
clear all
num = 12;               % accurate = 16, snap-to = 12
option = 2;             % 1 = n_swipes, 2 = sharing score, 3 = swipe accuracy ratio
destination = 'plates';   % n_swipes for 'plates', 'food' or 'inter' (inter-plates) destinations
gender = '';     % '' or 'Male' or 'Female' or 'compare'
severity = '';        % 'on' or ''
combine = 1;

[folder_loc,alt_folder_loc,file_loc,floc] = setup();
tab_sev = readtable([floc,'\eCRF.csv']);

[nam_save,~,ranked,~] = load_dataset(option,num,folder_loc,destination);

if combine == 0
    load('subject_details.mat')
    subject_details = subject_details_776;
elseif combine == 1
    addpath 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_krysiek'
    load('subject_details_combine_ond.mat')
    subject_details = subject_details_combine;
    if num == 16
        if strcmp(destination,'inter')
            extra = load([folder_loc,'\Ranking_Correlations\Data\OBJ_accurate_Krysiek.mat'],'nam_save','ranked'); % inter-plate sharing score
        else
            extra = load([folder_loc,'\Ranking_Correlations\Data\accurate_2only_Krysiek.mat'],'nam_save','ranked');
        end
    elseif num == 12
        if strcmp(destination,'inter')
            extra = load([folder_loc,'\Ranking_Correlations\Data\OBJ_snapto_redirect_Krysiek2.mat'],'nam_save','ranked'); % inter-plate sharing score
        else
            extra = load([folder_loc,'\Ranking_Correlations\Data\snapto_2only_Krysiek.mat'],'nam_save','ranked');
        end
    end
    ranked = [ranked;extra.ranked];
    nam_save = [nam_save,extra.nam_save];
end

[n_swipes,list] = swipe_analysis(num,file_loc,nam_save,destination);

saved=ones(1,length(subject_details)); saved(list)=zeros(1,length(list)); % create list of 1s delete those without adjs
[sets,months] = create_sets_months(subject_details,nam_save,saved,gender,tab_sev,severity);

sets = rmv_frm_sets(sets,months,ranked);

if ~strcmp(severity,'on')
    sets = rmv_ADHD(sets,subject_details,nam_save,saved);
end

if option == 1
    results = n_swipes;
elseif option == 2
    results = ranked';
%     results(results<-.3)=-0.3;
elseif option == 3
    [n_swipes,list] = swipe_analysis(num,file_loc,nam_save,destination);
    [n_swipes12,list12] = swipe_analysis(12,file_loc,nam_save,destination);
    results = n_swipes./n_swipes12;
end

if strcmp(gender,'compare') % compare male and female
    [all_sets,grps] = create_grps_allsets_gender(results,sets);
    boxplot_gender_cmpr(all_sets,grps,option,destination,num)
    [save_p] = significance_gender(results',sets);
    saveas(gcf,['Figures/boxplot_option',num2str(option),'_',num2str(num),'.png'])
%     %% save
%     figHandles = get(groot, 'Children');
%     for j = 1 : 4
%         set(0, 'currentfigure', figHandles(5-j));
%         saveas(gcf,['Figures/boxplot_option',num2str(option),'_',num2str(num),'_',num2str(j),'.png'])
%     end
else                        % display all or male/female genders
    for id = 1 : 3          % TD, ASD, OND
        if strcmp(severity,'on')
            sev_choice = id-1;
        else
            sev_choice = 0;
        end
        plot_results(results',months,sets,id,option,destination,num,gender,sev_choice)
    end
    
    [all_sets,grps] = create_grps_allsets(results,sets);
    
    Plot_boxplots(all_sets,grps,option,destination,num,gender,severity)
    
    [save_p,combos] = significance_check(results,sets);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [folder_loc,alt_folder_loc,file_loc,floc] = setup()
    folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
    alt_folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data';
    file_loc = [folder_loc,'\adjs\adj_obj_end_accurate\']; % should match zone type
    floc=[alt_folder_loc,'\IQ_severity'];

    folder1 = [folder_loc,'\Set_allocate'];
    folder2 = [folder_loc,'\Plots'];
    folder3 = [folder_loc,'\Create_adj'];
    folder4 = [folder_loc,'\adjs\adj_obj_end_accurate'];
    folder5 = [folder_loc,'\Data'];
    folder6 = folder_loc;
    addpath(folder1,folder2,folder3,folder4,folder5,folder6)
end

function [nam_save,saved,ranked,list] = load_dataset(option,num,folder_loc,destination)
    saved = []; list = [];
%     if option == 1 || option == 3     % volume
%         if num == 16
%             load([folder_loc,'\Ranking_Correlations\Data\OBJ_end_accurate_bi.mat'],'list','nam_save','saved','ranked')
%         elseif num == 12
%             load([folder_loc,'\Ranking_Correlations\Data\OBJ_end_12zones.mat'],'list','nam_save','saved','ranked')
%         end
%     elseif option == 2  % proportion
    if num == 16
        if strcmp(destination,'inter')
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_accurate.mat'],'nam_save','ranked') % inter-plate sharing score
        else
            load([folder_loc,'\Ranking_Correlations\Data\accurate_2only.mat'],'nam_save','ranked')
        end
    elseif num == 12
        if strcmp(destination,'inter')
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_snapto_redirect2.mat'],'nam_save','ranked') % inter-plate sharing score
        else
            load([folder_loc,'\Ranking_Correlations\Data\snapto_2only.mat'],'nam_save','ranked')
        end
    end
%     end
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
%         [sets_sev] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,gender);
%         sets{2} = setdiff(sets{2},sets_sev{4});
    end
end


function [sets] = rmv_frm_sets(sets,months,ranked)
    exclude = find(months>72);%72); %% 75 months threshold
%     exclude = [exclude;find(ranked<-0.5)];
    % exclude = [exclude;489]; % no food-to-plate swipes recorded
    exclude = [exclude;find(months<30)]; %% 45 Months threshold
    for i = 1 : length(sets)
        rmv=find(ismember(sets{i},exclude')==1);
        sets{i}(rmv)=[];
    end
end

function [sets] = rmv_ADHD(sets,subject_details,nam_save,saved)
    %%%  Remove ADHD
    load('OND_details_+Krysiek.mat','OND_details')
    if length(sets)<8
        [sets_OND] = set_allocate_TYPE_OND(subject_details,OND_details,nam_save,saved);
        curr_set = sets{3};
        keep = [sets_OND{6}];%[sets_OND{2},sets_OND{3},sets_OND{4}];
        sets{3}=curr_set(ismember(curr_set,keep));
%         rmv = sets_OND{1};
%         sets{3}=curr_set(~ismember(curr_set,rmv));
    else
        [sets_OND] = set_allocate_TYPE_OND(subject_details,OND_details,nam_save,saved);
        curr_set = sets{3};
        keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
        sets{3}=curr_set(ismember(curr_set,keep));
        curr_set = sets{7};
        sets{7}=curr_set(ismember(curr_set,keep));
    end
end

function [all_sets,grps] = create_grps_allsets_gender(results,sets)
    iter=0; all_sets=[]; grps=[];
    order = [1,5,2,6,3,7];
    
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
    n_swipes = zeros(1,length(nam_save));
    f_num = 0;
    
    list=[];
    for jj = 1:length(nam_save)
        adj=zeros(num,num);
        file_id = ['subject_',nam_save{jj},'.mat'];
    
        if isfile([file_loc,file_id]) || isfile(['C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_krysiek\',file_id])
            f_num = f_num + 1;
            load(file_id,'adj')
            if num == 12
                [adj] = adj_snap2zones(adj,num);
            end
%             adj=adj(1:num,1:num);
        end
    
        if max(adj(:))==0
            list=[list,jj];
        end
    
        if strcmp(destination,'plates')
            n_swipes(jj) = sum(adj(2,[4,5,6,7]));
        elseif strcmp(destination,'food')
            n_swipes(jj) = sum(adj(2,2));
%             n_swipes(jj) = sum(adj(2,:));
%             temp=diag(adj);
% %             sv_2 = temp(2);
%             temp(2)=0;
%             n_swipes(jj) = sum(temp);
        elseif strcmp(destination,'inter')
            n_swipes(jj) = sum(adj(4,[5,6,7]))+sum(adj(5,[4,6,7]))+sum(adj(6,[4,5,7]))+sum(adj(7,[4,5,6]));
        end
    end
end

function [adj] = adj_snap2zones(adj,num)
    % rewire adjacency from 16 to 12 zones
    init=adj;
%     allow=[2,4,5,6,7];
    for it = 2
        adj(it,4:7)=adj(it,4:7)+adj(it,13:16);    % reconnect to 4-7
        adj(it,13:16)=zeros(1,4);                     % remove non-food connections
    end
%     adj = adj(1:12,1:12);
    %% remove zn 4-7 incoming except from 2
    allow=[2];%,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            adj(it,13:16) = adj(it,4:7)+adj(it,13:16);
            adj(it,4:7)=zeros(1,4);                     % remove non-food connections
        end
    end
%     adj=adj-diag(diag(adj));
    
    bweight=1;
%     [adj] = NNR_adj_conns_OBJ2(adj,bweight);
end

function [] = plot_results(results,months,sets,id,option,destination,num,gender,sev_choice)
    figure;
    x=months(sets{id});
    y=results(sets{id});
    [pf,S] = polyfit(x,y,2);
    % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
    [y_fit,delta] = polyval(pf,x,S);
    mean(delta)
    % Plot the original data, linear fit, and 95% prediction interval y�2?.
    [~,Ind]=sort(x,'asc');
    [~,p] = corr(x,y,'Type','Spearman');
    
    if p<10.05
        linetype='-';
    end
%     if p<0.01
%         linetype='-';
%     elseif p<0.05
%         linetype='--';
%     else 
%         linetype=':';
%     end

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
%             clr = [.4 .4 .4];
        elseif id == 2
            clr = [0.8500, 0.3250, 0.0980];
        elseif id == 3
%             clr = [0.9290, 0.6940, 0.1250];
            clr = [0.47,0.67,0.19];
        elseif id == 4
            clr = [0.4940, 0.1840, 0.5560];
        end
    end

%     tempclr=[.5 .5 .5];
    scatplot=scatter(x(Ind),y(Ind),55,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatplot.MarkerFaceAlpha = .3;
    scatplot.MarkerEdgeAlpha = .5;
    hold on
    if p<10.05
        plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)
%         patch([x(Ind)' fliplr(x(Ind)')], [(y_fit(Ind)+1*delta(Ind))' fliplr((y_fit(Ind)-1*delta(Ind))')], [0, .6, .77],'FaceColor',clr, 'FaceAlpha',0.15, 'EdgeColor','none')
%         plot(x(Ind),y_fit(Ind),'-','color',clr,'LineWidth',3)
        plot(x(Ind),y_fit(Ind)+1*delta(Ind),linetype,'color',clr,'LineWidth',1.5,'LineStyle',':')
        plot(x(Ind),y_fit(Ind)-1*delta(Ind),linetype,'color',clr,'LineWidth',1.5,'LineStyle',':')
    end

    if p<0.001
        text(0.05,.97,['p = ',sprintf('%1.1e', p)],'Units','normalized','fontsize', 11)
        stars_only(3)
    elseif p<0.01
        text(0.05,.97,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
        stars_only(2)
    elseif p<0.05
        text(0.05,.97,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
        stars_only(1)
    elseif p<10%0.1
        text(0.05,.97,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    end
    
    if option == 2
%         axis([28 74 -.5 .25])
        axis([28 74 -.3 .25])
    elseif option == 1
%        axis([26 83 min(y) max(y)])
%         axis([28 74 0 120])
        axis([28 74 0 141])
%         axis([28 74 0 60])
        if strcmp(destination,'inter') 
        	axis([28 74 0 30])
        end
    elseif option == 3
        axis([28 74 0 1])
    end

    xlabel('Age (months)','fontsize', 11)
    if option == 1
        if strcmp(destination,'plates')
            if num == 16
                ylabel('No. of swipes (plates)','fontsize',14)
%                 ylabel('No. of swipes (all)','fontsize',14)
                nam = 'n_swipes_16_';
            elseif num == 12
                ylabel('No. of swipes','fontsize',14)
%                 ylabel('No. of swipes (inter-plates)','fontsize',14)
                nam = 'n_swipes_12_';
            end
        elseif strcmp(destination,'food')
            ylabel('No. of swipes (zone 2 only)','fontsize',14)
            nam = 'n_swipes_food_';
        elseif strcmp(destination,'inter')
            ylabel('No. of swipes (+ inter-plate)','fontsize',14)
            nam = 'n_swipes_inter_';
        end
    elseif option == 2
        if num == 16
            ylabel('Sharing score (plates)','fontsize',14)
%             ylabel('Sharing score (plates)','fontsize',14)
            nam = 'sharing_16_';
        elseif num == 12
            ylabel('Sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Sharing score (+ inter-plate)','fontsize',14)
            end
            nam = 'sharing_12_';
        end
    elseif option == 3
        ylabel('Swipe accuracy ratio','fontsize',14)
        nam = 'ratio_';
    end
%     legend('subject','2nd order fit','Location','SouthEast','Orientation','horizontal','fontsize',14)
    
    if id == 1
        titletxt = 'TD';
    elseif sev_choice == 1
        titletxt = 'ASD - Level 1';
    elseif sev_choice == 2
        titletxt = 'ASD - Level 2';
    elseif sev_choice == 3
        titletxt = 'ASD - Level 3';
    elseif id == 2
        titletxt = 'ASD';
    elseif id == 3
        titletxt = 'OND';
    end

    if strcmp(gender,'')
        title(titletxt)
    else
        title([titletxt,' - ',gender])
    end

    f=gcf;
    f.Position = [403,340,330,313];

    if strcmp(destination,'food')
        saveas(gcf,['Figures/food_',nam,titletxt,'_',gender,'.png'])
    else
        saveas(gcf,['Figures/',nam,titletxt,'_',gender,'.png'])
    end
end

function [] = stars_line(n_stars,height,strt,nd,drp)
    drp2 = 3.5*drp/2;
    shft = 0.105;
    hold on
    if n_stars == 3
%         scatter((nd-strt)/2+(strt-0.18),height,'pk','filled')
%         scatter((nd-strt)/2+(strt),height,'pk','filled')
%         scatter((nd-strt)/2+(strt+0.18),height,'pk','filled')
        text((nd-strt)/2+(strt-0.18)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt+0.18)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
    elseif n_stars == 2
        text((nd-strt)/2+(strt+0.09)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt-0.09)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
%         scatter((nd-strt)/2+(strt+0.09),height,'pk','filled')
    elseif n_stars == 1
%         scatter((nd-strt)/2+(strt),height,'pk','filled')
        text((nd-strt)/2+(strt)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
    end
    hold on
    plot([strt,nd],[height-drp,height-drp],'k')
    plot([strt,strt],[height-drp,height-drp2],'k')
    plot([nd,nd],[height-drp,height-drp2],'k')
end

function [] = stars_line_cmpr(n_stars,height,strt,nd,drp,drp2)
%     drp2 = 1*drp/2;%3.5*
    scale = 3;%3
    chng=.305;%chng = .45;%chng = .2;%
    shft=.105;%shft = .2;%shft = .1;%
    hold on
    if n_stars == 3
%         scatter((nd-strt)/2+(strt-0.18),height,'pk','filled')
%         scatter((nd-strt)/2+(strt),height,'pk','filled')
%         scatter((nd-strt)/2+(strt+0.18),height,'pk','filled')
        text((nd-strt)/2+(strt-chng-shft),height+drp2*scale,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt-shft),height+drp2*scale,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt+chng-shft),height+drp2*scale,'$$\star$$','Interpreter', 'latex','FontSize',18)
    elseif n_stars == 2
        text((nd-strt)/2+(strt+chng/2-shft),height+drp2*scale,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt-chng/2-shft),height+drp2*scale,'$$\star$$','Interpreter', 'latex','FontSize',18)
%         scatter((nd-strt)/2+(strt+0.09),height,'pk','filled')
    elseif n_stars == 1
%         scatter((nd-strt)/2+(strt),height,'pk','filled')
        text((nd-strt)/2+(strt-shft),height+drp2*scale,'$$\star$$','Interpreter', 'latex','FontSize',18)
    end
    hold on
    plot([strt,nd],[height-drp,height-drp],'k')
    plot([strt,strt],[height-drp,height-drp-drp2],'k')
    plot([nd,nd],[height-drp,height-drp-drp2],'k')
end

function [] = stars_only(n_stars)
    hold on
    gap =0.055;
    x=0.18; y=1.03;
    if n_stars == 3
        text(x-gap,y,'$$\star$$','Interpreter', 'latex','Units','normalized','FontSize',18)
        text(x,y,'$$\star$$','Interpreter', 'latex','Units','normalized','FontSize',18)
        text(x+gap,y,'$$\star$$','Interpreter', 'latex','Units','normalized','FontSize',18)
    elseif n_stars == 2
        text(x-gap/2,y,'$$\star$$','Interpreter', 'latex','Units','normalized','FontSize',18)
        text(x+gap/2,y,'$$\star$$','Interpreter', 'latex','Units','normalized','FontSize',18)
    elseif n_stars == 1
        text(x,y,'$$\star$$','Interpreter', 'latex','Units','normalized','FontSize',18)
    end
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
        colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.47,0.67,0.19;0.4940, 0.1840, 0.5560];
%         colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
    end
        
    for i = 1:4
        Ind=find(ic==i);
        scatter(ic(Ind),all_sets(Ind),[],colors(i,:),'filled','MarkerFaceAlpha',0.15,'jitter','on','jitterAmount',0.15);
    end
    
    xticklabels({'TD','ASD','OND'});
%     xticklabels({'TD','ASD 1','ASD 2','ASD 3'});
    
    if option == 1 
        if strcmp(destination,'plates')
            if num == 16
                ylabel('No. of swipes (plates)','fontsize',14)
%                 ylabel('No. of swipes (all)','fontsize',14)
            elseif num == 12
                ylabel('No. of swipes','fontsize',14)
%                 ylabel('No. of swipes (inter-plates)','fontsize',14)
            end
        elseif strcmp(destination,'food')
            ylabel('No. of swipes (zone 2 only)','fontsize',14)
        elseif strcmp(destination,'inter')
            ylabel('No. of swipes (inter-plate)','fontsize',14)
        end
    elseif option == 2
        if num == 16
            ylabel('Sharing score (plates)','fontsize',14)
        elseif num == 12
            ylabel('Sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Sharing score (+ inter-plate)','fontsize',14)
            end
        end
    elseif option == 3
        ylabel('Swipe accuracy ratio','fontsize',14)
    end
    
    box off

    if strcmp(destination,'food')
        axis([0.5 3.5 -3 164])
    elseif strcmp(destination,'inter') && option == 1
        axis([0.5 3.5 0 30])
    elseif option == 2
        axis([0.5 3.5 -.3 .25])%.39])
    elseif option == 1
        if num == 16
%             axis([0.5 4.5 -3 145])
            axis([0.5 3.5 0 141])
        elseif num == 12
%             axis([0.5 4.5 -3 172])
%             axis([0.5 4.5 -3 60])
            axis([0.5 3.5 0 141])
        end
    elseif option == 3
        axis([0.5 3.5 0 1.2])
    end

    %% Plot significance stars
    if strcmp(destination,'food')
        heights = [135,143,151,130];
        drp = 3;
    elseif strcmp(severity,'on')
        if option == 2
            heights = [.4,.35,.3];
            drp = .0225;
        elseif option == 1
%             heights = [125,120,115];
            heights = [155,147.5,140];
            drp = 3;
        end
    elseif option == 2
        heights = [.28,.33,.39,.28];%,.25];
%         heights = [.41,.345,.45,.31];
        drp = .02;
    elseif option == 1
        if num == 16
            heights = [122,130,140,122];
        elseif num == 12
            heights = [155,162,171,155];
        end
        drp = 3;
    elseif option == 3
        heights = [1.05,1.105,1.17,1.05];
        drp = 0.025;
    end
 
%     if strcmp(destination,'food')
%         stars_line(3,heights(1),1,2,drp) % 3 stars,h,1,2,1
%         stars_line(2,heights(2),1,3,drp) % 3 stars,h,1,2,1
%         stars_line(1,heights(3),1,4,drp) % 2 stars,h,1,3,2
%         stars_line(3,heights(4),2,4,drp) % 2 stars,h,3,4,1
%     elseif strcmp(severity,'on')
%         stars_line_cmpr(2,heights(1),1,2,drp) % 
%         stars_line_cmpr(1,heights(2),2,3,drp) % 
%         stars_line_cmpr(1,heights(3),3,4,drp) % 
%     else
%         stars_line(3,heights(1),1,2,drp) % 3 stars,h,1,2,1
%         stars_line(3,heights(2),1,3,drp) % 3 stars,h,1,2,1
%         stars_line(3,heights(3),2,4,drp) % 2 stars,h,1,3,2
%         stars_line(3,heights(4),3,4,drp) % 2 stars,h,3,4,1
% %         stars_line(1,heights(5),2,3,drp) 
%     end
       
    % legend('TD','ASD','OND*','ONDE','Orientation','horizontal')
    
    f.Position = [403,340,300,313];

    title(gender)

    if strcmp(destination,'food')
        saveas(gcf,['Figures/boxplot_food_',gender,'.png'])
    elseif strcmp(destination,'inter')
        saveas(gcf,['Figures/boxplot_option',num2str(option),'_inter_',num2str(num),'_',gender,'.png'])       
    else
        saveas(gcf,['Figures/boxplot_option',num2str(option),'_',num2str(num),'_',gender,'.png'])
    end
end

function [save_p,combos] = significance_check(results,sets)
    % pairwise significance check
    combos=nchoosek([1,2,3],2);
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
    numopt = [1,5;2,6;3,7];
    for j = 1 : 3
        num=numopt(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval = kruskalwallis([results(sets{num(1)});results(sets{num(2)})],len_rankeds,'off');
        save_p(1,j) = pval;
%         if pval > .001
%             text(0.05+(j-1)*1/4,.92,['p = ',sprintf('%1.4f', pval)],'Units','normalized','fontsize', 10)
%         else
%             text(0.05+(j-1)*1/4,.92,['p = ',sprintf('%1.1e', pval)],'Units','normalized','fontsize', 10)
%         end
    end
end

% function [save_p] = significance_gender_separate(results,sets)
%     save_p = zeros(1,4);
%     numopt = [1,5;2,6;3,7;4,8];
%     figHandles = get(groot, 'Children');
%     for j = 1 : 4
%         set(0, 'currentfigure', figHandles(5-j));
%         num=numopt(j,:);
%         len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
%         pval = kruskalwallis([results(sets{num(1)});results(sets{num(2)})],len_rankeds,'off');
%         save_p(1,j) = pval;
%         if pval > .001
%             text(0.375,.92,['p = ',sprintf('%1.4f', pval)],'Units','normalized','fontsize', 10)
%         else
%             text(0.375,.92,['p = ',sprintf('%1.1e', pval)],'Units','normalized','fontsize', 10)
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = boxplot_gender_cmpr(all_sets,grps,option,destination,num)

    f=figure;
    b=boxplot(all_sets,grps,'Notch','on','Color',[.5,.5,.25]);
    set(b,'LineWidth',1.5)
%     h = findobj(gca,'Tag','Box');

     %% Add scatter points
    hold on
    [~,~,ic]=unique(grps,'stable');
    
    colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.47,0.67,0.19;0.4940, 0.1840, 0.5560];
    for j = 1:12
        i=ceil(j/3);
        Ind=find(ic==j);
        scatter(ic(Ind),all_sets(Ind),[],colors(i,:),'filled','MarkerFaceAlpha',0.08,'jitter','on','jitterAmount',0.15);
    end
    
    hAx=gca; 
    xt=hAx.XTick;                   % retrieve ticks
    xt(3:3:end)=[];
    hAx.XTick=xt;                   % clear some
    
    entries ={'F','M','F','M','F','M'};
    xticklabels(entries)

    if strcmp(destination,'food')
        axis([0 9 -5 150])
    elseif option == 1
        if num == 16
            axis([0 9 -5 141])
        elseif num == 12
            axis([0 9 -5 141])
        end
    elseif option == 2
        axis([0 9 -.3 .25])
    elseif option == 3
        axis([0 9 0 1])
    end
    
    %% Plot significance stars
    % height=30; n_stars = 2; drp = 1.75;
    % stars_line(n_stars,height,1,2,drp) % 3 stars,h,1,2,1
    % height=30; n_stars = 1; drp = 1.75;
    % stars_line(n_stars,height,10,11,drp) % 3 stars,h,1,2,1
    box off
    f.Position = [403,340,574,313];
%     f.Position = [403,340,400,313];
    xlabel('Gender')

    if option == 1
        if strcmp(destination,'plates')
            if num == 16
                ylabel('No. of swipes (plates)','fontsize',14)
%                 ylabel('No. of swipes (all)','fontsize',14)
            elseif num == 12
                ylabel('No. of swipes','fontsize',14)
            end
        elseif strcmp(destination,'food')
            ylabel('No. of swipes (zone 2 only)','fontsize',14)
        end
    elseif option == 2
        if num == 16
            ylabel('Sharing score (plates)','fontsize',14)
        elseif num == 12
            ylabel('Sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Sharing score (+ inter-plate)','fontsize',14)
            end
        end
    elseif option == 3
        ylabel('Swipe accuracy ratio','fontsize',14)
    end

    %% Plot significance stars
    if strcmp(destination,'food')
        heights = 135.*ones(1,4);
        drp = 3;
    elseif option == 2
    %     heights = [.295,.345,.41,.295];
        heights = .25.*ones(1,4);
        drp = .01;
        drp2=3*drp/2;
    elseif option == 1
        if num == 16
            heights = 122.*ones(1,4);
        elseif num == 12
            heights = 139.*ones(1,4);
        end
        drp = 2;
        drp2=4.65*drp/2;
    elseif option == 3
        heights = 1.05.*ones(1,4);
        drp = 0.025;
    end
 
%     if strcmp(destination,'food')
%         stars_line_cmpr(3,heights(1),1,2,drp) % 3 stars,h,1,2,1
%         stars_line_cmpr(2,heights(2),4,5,drp) % 3 stars,h,1,2,1
%         stars_line_cmpr(1,heights(3),7,8,drp) % 2 stars,h,1,3,2
%     else
%         stars_line_cmpr(2,heights(1),1,2,drp,drp2) % 3 stars,h,1,2,1
% %         stars_line_cmpr(2,heights(2),4,5,drp,drp2) % 3 stars,h,1,2,1
% %         stars_line_cmpr(3,heights(3),7,8,drp) % 2 stars,h,1,3,2
% %         stars_line_cmpr(2,heights(4),10,11,drp,drp2) % 2 stars,h,3,4,1
%     end
    if option == 1
        textadd = 'swipes';
    else
        textadd = 'sharing';
    end

    saveas(gcf,['Figures/compare_',num2str(num),'_',textadd,'.png'])
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
