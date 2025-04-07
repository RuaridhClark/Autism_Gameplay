%% Main analysis script for sharing gameplay
clear all

option = 1;             % 1 = n_swipes, 2 = sharing score, 3 = sharing score difference
destination = 'plates'; % n_swipes for food delivery swipes ('plates') or inter-plate swipes ('inter')
sex = '';               % '' or 'Male' or 'Female' or 'compare' (switch for analysing sex groupings)
severity = '';          % 'on' or '' (switch for analysing severity levels in ASD)
combine = 1;            % Combine both pre-trial and trial data when combine = 1
bweight = 0.01;         % Edge weight for complete graph

[num,alt_folder_loc,file_loc] = setup(option,destination);

%%% Load data %%%
[name_save,~,ranked,~] = load_trial_data(num,destination,bweight); % load trial data
[name_save,ranked,subject_details] = load_pretrial_data(name_save,ranked,combine,num,destination,bweight); % load pretrial data

%%% Count swipes %%%
[n_swipes,saved] = swipe_analysis(num,file_loc,name_save,destination,length(subject_details)); % list records subjects with adjacency matrix

%%% Create group sets %%%
[sets,months,n] = create_sets_months(subject_details,name_save,saved,sex,severity);

%%% Results variable %%%
[results] = results_variable(option,n_swipes,ranked,destination);

%%% Plot results %%%
[p_values] = plot_options(results,sets,option,destination,num,months,sex,severity,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num,alt_folder_loc,file_loc] = setup(option,destination)
    if option == 2 || strcmp(destination,'plates')
        num = 12;
    else
        num = 16;
    end

    alt_folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data';

    file_loc = '..\adjs\adj_zones\';
    folder1 = '..\Data';
    folder2 = 'Set_allocate';
    addpath(file_loc,folder1,folder2)
end

function [name_save,saved,ranked,list] = load_trial_data(num,destination,bweight)
    [w_string] = save_title('trial',bweight);

    saved = []; list = [];
    if num == 12 || num == 16
        if strcmp(destination,'inter')
            load(['Data\SS_ext_inter',w_string,'.mat'],'name_save','ranked')
        else
            load(['Data\SS_ext',w_string,'.mat'],'name_save','ranked')
        end
    end
end

function [name_save,ranked,subject_details] = load_pretrial_data(name_save,ranked,combine,num,destination,bweight)
    if combine == 0
        load('subject_details_trial.mat')
    elseif combine == 1
        addpath '..\adjs\adj_zones'
        load('subject_details_combine.mat')
        subject_details = subject_details_combine;
        [w_string] = save_title('pretrial',bweight);
        if num == 12 || num == 16
            if strcmp(destination,'inter')
                extra = load(['Data\SS_ext_inter',w_string,'.mat'],'name_save','ranked');
            else
                extra = load(['Data\SS_ext',w_string,'.mat'],'name_save','ranked');
            end
        end
        ranked = [ranked;extra.ranked];
        name_save = [name_save,extra.name_save];
    end
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

function [sets,months,n] = create_sets_months(subject_details,name_save,saved,sex,severity)
    [months] = list_AGE(subject_details,name_save,saved);
    if size(months,1)<size(months,2)
        months=months';
    end
    if strcmp(severity,'on')
        [sets] = set_allocate_severity(subject_details,name_save,saved,sex);
    elseif strcmp(sex,'Male')
        [tmp_sets] = set_allocate_sex_TYPE(subject_details,name_save,saved);
        sets = tmp_sets(5:8);
    elseif strcmp(sex,'Female')
        [tmp_sets] = set_allocate_sex_TYPE(subject_details,name_save,saved);
        sets = tmp_sets(1:4);
    elseif strcmp(sex,'compare')
        [sets] = set_allocate_sex_TYPE(subject_details,name_save,saved);
    else
        [sets] = set_allocate(subject_details,name_save,saved);
    end

    sets = rmv_frm_sets(sets,months); % age range 30 to 72 months

    %%% Adjust sets for ASD severity %%%
    if ~strcmp(severity,'on') % adjust sets to include only ASD severity levels and WP
        sets = rmv_ADHD(sets,subject_details,name_save,saved);
        n=3;
    else 
        n=4;
    end
end

function [sets] = rmv_frm_sets(sets,months)
    exclude = find(months>72); % 72 month threshold
    exclude = [exclude;find(months<30)]; %% 30 Months threshold
    for i = 1 : length(sets)
        rmv=find(ismember(sets{i},exclude')==1);
        sets{i}(rmv)=[];
    end
end

function [sets] = rmv_ADHD(sets,subject_details,name_save,saved)
    %%%  Remove ADHD
    if length(sets)<8
        [sets_OND] = set_allocate_TYPE_OND(subject_details,name_save,saved);
        curr_set = sets{3};
        keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
        sets{3}=curr_set(ismember(curr_set,keep));
    else
        [sets_OND] = set_allocate_TYPE_OND(subject_details,name_save,saved);
        curr_set = sets{3};
        keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
        sets{3}=curr_set(ismember(curr_set,keep));
        curr_set = sets{7};
        sets{7}=curr_set(ismember(curr_set,keep));
    end
end

function [all_sets,grps] = create_grps_allsets_sex(results,sets)
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
        iter=iter+1;
        vals = results(sets{num});
        if isempty(vals)
            vals=0;
        end
        all_sets = [all_sets,vals];
        grps = [grps,ones(1,length(vals))*iter];
    end	
end

function [n_swipes,saved] = swipe_analysis(num,file_loc,name_save,destination,len_subjects) 
    n_swipes = zeros(1,length(name_save));
    f_num = 0;
    
    list=[];
    for jj = 1:length(name_save)
        adj=zeros(num,num);
        file_id = ['subject_',name_save{jj},'.mat'];
    
        if isfile([file_loc,file_id])
            f_num = f_num + 1;
            load([file_loc,file_id],'adj')
            if num == 12
                [adj] = adj_snap2zones(adj,num);
            end
        end
    
        if max(adj(:))==0
            list=[list,jj];
        end
    
        if strcmp(destination,'plates')
            n_swipes(jj) = sum(adj(2,[4,5,6,7]));
        elseif strcmp(destination,'food')
            n_swipes(jj) = sum(adj(2,2));
        elseif strcmp(destination,'inter')
            if num == 12
                n_swipes(jj) = sum(adj(4,[5,6,7]))+sum(adj(5,[4,6,7]))+sum(adj(6,[4,5,7]))+sum(adj(7,[4,5,6]));
            elseif num == 16
                n_swipes(jj) = sum(adj(4,[5,6,7,14,15,16]))+sum(adj(5,[4,6,7,13,15,16]))+sum(adj(6,[4,5,7,13,14,16]))+sum(adj(7,[4,5,6,13,14,15]));
            end
        elseif strcmp(destination,'total')
            n_swipes(jj) = sum(adj(:));
        elseif strcmp(destination,'evenness_indirect')
            V1 = sum(adj(5,[4,13]))+sum(adj(6,[4,13]))+sum(adj(7,[4,13])) + sum(adj(2,[4,13]));
            V2 = sum(adj(4,[5,14]))+sum(adj(6,[5,14]))+sum(adj(7,[5,14])) + sum(adj(2,[5,14]));
            V3 = sum(adj(4,[6,15]))+sum(adj(5,[6,15]))+sum(adj(7,[6,15])) + sum(adj(2,[6,15]));
            V4 = sum(adj(4,[7,16]))+sum(adj(5,[7,16]))+sum(adj(6,[7,16])) + sum(adj(2,[7,16]));
            n_swipes(jj) = std([V1,V2,V3,V4])/sum([V1,V2,V3,V4]);
        elseif strcmp(destination,'evenness_direct')
            V1 = sum(adj(2,[4,13]));
            V2 = sum(adj(2,[5,14]));
            V3 = sum(adj(2,[6,15]));
            V4 = sum(adj(2,[7,16]));
            n_swipes(jj) = std([V1,V2,V3,V4])/sum([V1,V2,V3,V4]); 
        elseif strcmp(destination,'attentive_indirect')
            adj = adj - diag(diag(adj));
            V1 = sum(adj(5,[4,13]))+sum(adj(6,[4,13]))+sum(adj(7,[4,13])) + sum(adj(2,[4,13]));
            V2 = sum(adj(4,[5,14]))+sum(adj(6,[5,14]))+sum(adj(7,[5,14])) + sum(adj(2,[5,14]));
            V3 = sum(adj(4,[6,15]))+sum(adj(5,[6,15]))+sum(adj(7,[6,15])) + sum(adj(2,[6,15]));
            V4 = sum(adj(4,[7,16]))+sum(adj(5,[7,16]))+sum(adj(6,[7,16])) + sum(adj(2,[7,16]));
            n_swipes(jj) = sum([V1,V2,V3,V4])/sum(adj(:));
        elseif strcmp(destination,'attentive_direct')
            adj = adj - diag(diag(adj));
            n_swipes(jj) = sum(adj(2,[4,5,6,7]))/sum(adj(:));
        end
    end

    saved=ones(1,len_subjects); 
    saved(list)=zeros(1,length(list)); % create list of 1s, zero those without adjs
end

function [results] = results_variable(option,n_swipes,ranked,destination)
    if option == 1
        results = n_swipes;
    elseif option == 2
        results = ranked';
        if strcmp(destination,'plates')
            save('results_plates.mat','results')
        elseif strcmp(destination,'inter')
            save('results_inter.mat','results')
        end
    elseif option == 3
        % Sharing score difference
        load('results_inter.mat')
        results_i = results;
        load('results_plates.mat')
        results_p = results;
        results=results_i-results_p;
    end
end

function [adj] = adj_snap2zones(adj,num)
    % rewire adjacency from 16 to 12 zones
    for it = 2
        adj(it,4:7)=adj(it,4:7)+adj(it,13:16);    % reconnect to 4-7
        adj(it,13:16)=zeros(1,4);                     % remove non-food connections
    end

    %% remove zn 4-7 incoming except from 2
    allow=2;
    for it = 1:num
        if ~ismember(it,allow)
            adjust = adj(it,4:7);
            adj(it,13:16)=adj(it,13:16)+adjust;    % reconnect to 13-16
            if ismember(it,4:7)
                adj(it,it+9)=adj(it,it+9)-adjust(ismember(4:7,it));
            end

            % Remove connections (not self-loops)
            list = 4:7;
            listedit =list(~ismember(4:7,it));
            adj(it,listedit)=zeros(1,length(listedit));                     % remove non-food connections
        end
    end
end

function [] = plot_results(results,months,sets,id,option,destination,num,sex,severity)

    [sev_choice] = severity_switch(severity,id);

    x=months(sets{id}); 
    y=results(sets{id});

    x=x(~isnan(y)); y=y(~isnan(y));

    [pf,S] = polyfit(x,y,2);
    [y_fit,delta] = polyval(pf,x,S); % Evaluate the first-degree polynomial fit in p at the points in x
    
    [r,p] = corr(x,y,'Type','Spearman'); % Spearman correlation
    
    plot_and_save(x,y,sev_choice,id,y_fit,delta,r,p,option,destination,num,sex);
    
end

function [sev_choice] = severity_switch(severity,id)
    if strcmp(severity,'on')
        sev_choice = id-1;
    else
        sev_choice = 0;
    end
end

function [] = plot_and_save(x,y,sev_choice,id,y_fit,delta,r,p,option,destination,num,sex)

    [clr] = set_colour(sev_choice,id);
    linetype='-';
    [~,Ind]=sort(x,'asc');

    figure;
    scatplot=scatter(x(Ind),y(Ind),55,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatplot.MarkerFaceAlpha = .3;
    scatplot.MarkerEdgeAlpha = .5;
    hold on
    plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)
    plot(x(Ind),y_fit(Ind)+1*delta(Ind),linetype,'color',clr,'LineWidth',1.5,'LineStyle',':')
    plot(x(Ind),y_fit(Ind)-1*delta(Ind),linetype,'color',clr,'LineWidth',1.5,'LineStyle',':')

    p_r_stars(p,r)

    xlabel('Age (months)','fontsize', 11)

    [name] = add_ylabel(option,destination,num);
    
    [titletxt] = add_title(sev_choice,id,sex);

    f=gcf;
    f.Position = [403,340,330,313];

    set_axis(option,destination)

    save_figure(destination,name,titletxt,sex,sev_choice)
end

function [] = Plot_boxplots(all_sets,grps,option,destination,num,sex,severity,p_values,combos)

    f=figure;
    b=boxplot(all_sets,grps,'Notch','on','Color',[.5,.5,.25]);
    set(b,'LineWidth',1.5)
    
    %% Add scatter points
    hold on
    [~,~,ic]=unique(grps,'stable');
    
    if strcmp(severity,'on')
        colors = [0, 0.4470, 0.7410;.93,.69,.13; .85,.33,.1; .64,.08,.18];
    else
        colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.47,0.67,0.19;0.4940, 0.1840, 0.5560];
    end
        
    for i = 1:4
        Ind=find(ic==i);
        scatter(ic(Ind),all_sets(Ind),[],colors(i,:),'filled','MarkerFaceAlpha',0.15,'jitter','on','jitterAmount',0.15);
    end
    
    if strcmp(severity,'on')
        xticklabels({'WP','ASD 1','ASD 2','ASD 3'});
    else
        xticklabels({'WP','ASD','OND'});
    end
    
    [name] = add_ylabel(option,destination,num);
    
    box off
    
    f.Position = [403,340,300,313];

    set_axis_boxplot(option,destination);

    title(sex)

    %% Plot significance stars    
    auto_star_plot(ic,all_sets,p_values,combos,option,severity)

    save_figure(destination,['boxplot_',name],'',sex,'')
end

function [] = p_r_stars(p,r)
    if p<0.001
        text(0.05,.97+.025,['p = ',sprintf('%1.1e', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,strrep(['r = ',sprintf('%1.2f', r)],'-','−'),'Units','normalized','fontsize', 11)
        stars_only(3)
    elseif p<0.01
        text(0.05,.97+.025,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,strrep(['r = ',sprintf('%1.2f', r)],'-','−'),'Units','normalized','fontsize', 11)
        stars_only(2)
    elseif p<0.05
        text(0.05,.97+.025,['p = ',sprintf('%1.3f', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,strrep(['r = ',sprintf('%1.2f', r)],'-','−'),'Units','normalized','fontsize', 11)
        stars_only(1)
    elseif p<10
        text(0.05,.97+.025,['p = ',sprintf('%1.3f', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,strrep(['r = ',sprintf('%1.2f', r)],'-','−'),'Units','normalized','fontsize', 11)
    end
end

function [] = stars_line(n_stars,height,strt,nd,drp2)
    drp = 2*drp2/3.5;
    shft = 0.105;
    hold on
    if n_stars == 3
        text((nd-strt)/2+(strt-0.18)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt+0.18)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
    elseif n_stars == 2
        text((nd-strt)/2+(strt+0.09)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
        text((nd-strt)/2+(strt-0.09)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
    elseif n_stars == 1
        text((nd-strt)/2+(strt)-shft,height,'$$\star$$','Interpreter', 'latex','FontSize',18)
    end
    hold on
    plot([strt,nd],[height-drp,height-drp],'k')
    plot([strt,strt],[height-drp,height-drp2],'k')
    plot([nd,nd],[height-drp,height-drp2],'k')
end

function [] = stars_only(n_stars)
    hold on
    gap =0.055;
    x=0.18; y=1.03+.025;
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

function [name] = add_ylabel(option,destination,num)
    if option == 1
        if strcmp(destination,'plates')
            if num == 12
                ylabel('No. of food delivery swipes','fontsize',14)
                name = 'n_swipes_food_';
            end
        elseif strcmp(destination,'food')
            ylabel('No. of zone 2 only swipes ','fontsize',14)
            name = 'n_swipes_2_';
        elseif strcmp(destination,'inter')
            ylabel('No. of inter-plate swipes','fontsize',14)
            name = 'n_swipes_inter_';
        elseif strcmp(destination,'evenness_direct')
            ylabel('Direct standard deviation','fontsize',14)
            name = 'direct_std_';
        elseif strcmp(destination,'evenness_indirect')
            ylabel('Indirect standard deviation','fontsize',14)
            name = 'indirect_std_';
        elseif strcmp(destination,'total')
            ylabel('Total no. of swipes','fontsize',14)
            name = 'n_swipes_total_';
        elseif strcmp(destination,'attentive_direct')
            ylabel('Direct attentiveness ratio','fontsize',14)
            name = 'direct_attentive_';
        elseif strcmp(destination,'attentive_indirect')
            ylabel('Indirect attentiveness ratio','fontsize',14)
            name = 'indirect_attentive_';
        end
    elseif option == 2
        if num == 12
            ylabel('Direct sharing score','fontsize',14)
            name = 'sharing_direct_';
            if strcmp(destination,'inter')
                ylabel('Indirect sharing score','fontsize',14)
                name = 'sharing_indirect_';
            end
        end
    elseif option == 3
        ylabel('Sharing score difference','fontsize',14)
        name = 'sharing_diff_';
    end
end

function [titletxt] = add_title(sev_choice,id,sex)
    if id == 1
        titletxt = 'WP';
    elseif sev_choice == 1
        titletxt = '              ASD - Level 1';
    elseif sev_choice == 2
        titletxt = '              ASD - Level 2';
    elseif sev_choice == 3
        titletxt = '              ASD - Level 3';
    elseif id == 2
        titletxt = 'ASD';
    elseif id == 3
        titletxt = 'OND';
    end

    if strcmp(sex,'')
        title(titletxt)
    else
        title([titletxt,' - ',sex])
    end
end

function [clr] = set_colour(sev_choice,id)
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
            clr = [0.47,0.67,0.19];
        elseif id == 4
            clr = [0.4940, 0.1840, 0.5560];
        end
    end
end

function [] = set_axis(option,destination)
    if option == 2
        axis([28 74 -.65 .3])
    elseif option == 1
        if strcmp(destination,'inter') 
        	axis([28 74 0 60])
        elseif strcmp(destination,'food') 
            axis([28 74 0 164])
        elseif strcmp(destination,'total')
            axis([28 74 0 650])
        elseif contains(destination,'attentive')
            axis([28 74 0 1.15])
        elseif contains(destination,'evenness')
            axis([28 74 0 0.5])
        else
            axis([28 74 0 147])
        end
    elseif option == 3
        axis([28 74 -.15 0.5])
    end
end

function [] = set_axis_boxplot(option,destination)
    if option == 2
        axis([0.5 3.5 -.65 .3])
    elseif option == 1
        if strcmp(destination,'inter') 
        	axis([0.5 3.5 0 60])
        elseif strcmp(destination,'food') 
            axis([0.5 3.5 0 164])
        elseif strcmp(destination,'total')
            axis([0.5 3.5 0 650])
        elseif contains(destination,'attentive')
            axis([0.5 3.5 0 1.15])
        elseif contains(destination,'evenness')
            axis([.5 3.5 0 0.5])
        else
            axis([0.5 3.5 0 147])
        end
    elseif option == 3
        axis([0.5 3.5 -.15 0.5])
    end
end

function [] = save_figure(destination,name,titletxt,sex,sev_choice)
    if strcmp(destination,'food')
        saveas(gcf,['Figures/food_',name,titletxt,'_',sex,'.png'])
    elseif sev_choice>0
        saveas(gcf,['Figures/severity_',name,titletxt,'_',sex,'.png'])
    elseif strcmp(destination,'inter')
        saveas(gcf,['Figures/',name,titletxt,'_',sex,'_inter.png'])
    else
        saveas(gcf,['Figures/',name,titletxt,'_',sex,'_plates.png'])
    end
end

function [] = auto_star_plot(ic,all_sets,p_values,combos,option,severity)
    if strcmp(severity,'on')
        n=4;
    else
        n=length(p_values);
    end

    mv = zeros(1,n);
    for i = 1:n
        Ind=find(ic==i);
        mv(i) = max(all_sets(Ind));
    end

    if option == 2 || option == 3
        drp = max(mv)*0.175;
    elseif option == 1
        drp = max(mv)*0.05;
    end
    addh=drp*1.5;
    h_save=0;
    for i = 1 : length(p_values)
        j = length(p_values)+1-i;
        id1 = combos(j,1);
        id2 = combos(j,2);
        height = max([mv(id1:id2)])+addh;
        height = max([h_save+addh,height]);
        if p_values(j) < 0.001
            stars_line(3,height,id1,id2,drp)
        elseif p_values(j) < 0.01
            stars_line(2,height,id1,id2,drp)
        elseif p_values(j) < 0.05
            stars_line(1,height,id1,id2,drp)
        end
        if p_values(j) < 0.05
            h_save = height;
        end
    end
end

function [p_values,combos] = significance_check(results,sets,n)
    % Pairwise significance check and effect size calculation
    combos = nchoosek(1:n, 2); % All pairwise combinations
    p_values = zeros(1, size(combos, 1)); % Preallocate for p-values
    effect_sizes = zeros(1, size(combos, 1)); % Preallocate for effect sizes (e.g., Cohen's d)
    
    for j = 1:size(combos, 1)
        num = combos(j, :);
        
        % Group labels for Kruskal-Wallis test
        len_rankeds = [ones(length(sets{num(1)}), 1); 2 * ones(length(sets{num(2)}), 1)];
        pval_kw = kruskalwallis([results(sets{num(1)})'; results(sets{num(2)})'], len_rankeds, 'off');
    
        % Extract the data for the two groups
        group1 = results(sets{num(1)});
        group2 = results(sets{num(2)});
        
        % Perform the two-tailed unpaired Wilcoxon test
        [pval, ~] = ranksum(group1, group2);
        p_values(1, j) = pval;
        
        % Calculate effect size (Cohen's d) if the comparison is significant
        if pval < 0.05
            % Descriptive statistics
            mean1 = mean(group1);
            mean2 = mean(group2);
            std1 = std(group1);
            std2 = std(group2);
            n1 = length(group1);
            n2 = length(group2);
            
            % Pooled standard deviation
            pooled_std = sqrt(((n1 - 1) * std1^2 + (n2 - 1) * std2^2) / (n1 + n2 - 2));
            
            % Cohen's d
            effect_sizes(1, j) = (mean2 - mean1) / pooled_std;
        else
            effect_sizes(1, j) = NaN; % Not significant, no effect size calculated
        end
    end

    % Display significant comparisons and effect sizes
    disp('Pairwise comparisons:');
    for j = 1:size(combos, 1)
        if p_values(1, j) < 0.05
            fprintf('Comparison %d-%d: p = %.4f, Cohen''s d = %.4f\n', ...
                    combos(j, 1), combos(j, 2), p_values(1, j), effect_sizes(1, j));
        else
            fprintf('Comparison %d-%d: p = %.4f (not significant)\n', ...
                    combos(j, 1), combos(j, 2), p_values(1, j));
        end
    end
end

function [p_values] = significance_sex(results,sets)
    p_values = zeros(1,4);
    numopt = nchoosek([1,5,2,6,3,7],2);
    for j = 1 : length(numopt)
        num=numopt(j,:);

        % Extract the data for the two groups
        group1 = results(sets{num(1)});
        group2 = results(sets{num(2)});
        
        % Perform the two-tailed unpaired Wilcoxon test
        [pval, ~] = ranksum(group1, group2);

        p_values(1,j) = pval;
    end
end

function [ ] = boxplot_sex_cmpr(all_sets,grps,option,destination,num)

    f=figure;
    b=boxplot(all_sets,grps,'Notch','on','Color',[.5,.5,.25]);
    set(b,'LineWidth',1.5)

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

    box off
    f.Position = [403,340,574,313];
    xlabel('Sex')

    if option == 1
        if strcmp(destination,'plates')
            if num == 12
                ylabel('No. of food delivery swipes','fontsize',14)
            end
        elseif strcmp(destination,'food')
            ylabel('No. of swipes (zone 2 only)','fontsize',14)
        elseif strcmp(destination,'evenness')
            ylabel('Standard deviation','fontsize',14)
        end
    elseif option == 2
        if num == 12
            ylabel('Direct sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Indirect sharing score','fontsize',14)
            end
        end
    end
end

function [p_values] = plot_options(results,sets,option,destination,num,months,sex,severity,n)

    if strcmp(sex,'compare') % compare male and female
        [all_sets,grps] = create_grps_allsets_sex(results,sets);
        boxplot_sex_cmpr(all_sets,grps,option,destination,num)
        [p_values] = significance_sex(results',sets);
        axis([0.5 8.5 0 150])
        saveas(gcf,['Figures/Compare_',destination,'.png'])
    else                        % display all or male/female sex
        for id = 1 : n        % TD, ASD, OND
            plot_results(results',months,sets,id,option,destination,num,sex,severity)
        end
        
        [all_sets,grps] = create_grps_allsets(results,sets);
    
        [p_values,combos] = significance_check(results,sets,n);
    
        if strcmp(severity,'on')
            p_values = p_values([1,4,6]);
            combos = combos([1,4,6],:);
        end
        
        Plot_boxplots(all_sets,grps,option,destination,num,sex,severity,p_values,combos)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
