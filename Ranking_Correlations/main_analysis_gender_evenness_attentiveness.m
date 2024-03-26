% 
clear all
num = 16;               % accurate = 16 (only for no. of swipes), snap-to = 12
option = 1;             % 1 = n_swipes, 2 = sharing score, 3 = swipe accuracy ratio
destination = 'inter';   % n_swipes for 'plates', 'food' (food/table-zone) or 'inter' (inter-plates) destinations
sex = 'Female';     % '' or 'Male' or 'Female' or 'compare'
severity = '';        % 'on' or ''
combine = 1;
bweight = '_0_01';

[folder_loc,alt_folder_loc,file_loc,floc] = setup();
tab_sev = readtable([floc,'\eCRF.csv']);

[nam_save,~,ranked,~] = load_dataset(option,num,folder_loc,destination,bweight);
[nam_save,ranked,subject_details] = load_extra_data(nam_save,ranked,combine,option,num,folder_loc,destination,bweight);


[n_swipes,list] = swipe_analysis(num,file_loc,nam_save,destination);

saved=ones(1,length(subject_details)); saved(list)=zeros(1,length(list)); % create list of 1s delete those without adjs
[sets,months] = create_sets_months(subject_details,nam_save,saved,sex,tab_sev,severity);

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

n=3;
if strcmp(severity,'on')
    n=4;
end

if strcmp(sex,'compare') % compare male and female
    [all_sets,grps] = create_grps_allsets_sex(results,sets);
    boxplot_sex_cmpr(all_sets,grps,option,destination,num)
    [save_p] = significance_sex(results',sets);
    saveas(gcf,['Figures/boxplot_option',num2str(option),'_',num2str(num),'.png'])
%     %% save
%     figHandles = get(groot, 'Children');
%     for j = 1 : 4
%         set(0, 'currentfigure', figHandles(5-j));
%         saveas(gcf,['Figures/boxplot_option',num2str(option),'_',num2str(num),'_',num2str(j),'.png'])
%     end
else                        % display all or male/female sexs
    for id = 1 : n        % TD, ASD, OND
        if strcmp(severity,'on')
            sev_choice = id-1;
        else
            sev_choice = 0;
        end
        plot_results(results',months,sets,id,option,destination,num,sex,sev_choice,[])
    end
    
    [all_sets,grps] = create_grps_allsets(results,sets);

    [save_p,combos] = significance_check(results,sets,n);

    if strcmp(severity,'on')
        save_p = save_p([1,4,6]);
        combos = combos([1,4,6],:);
    end
    
    Plot_boxplots(all_sets,grps,option,destination,num,sex,severity,save_p,combos)

end

%%% ADDED - semi-partial correlation input

% % inter only
% [nam_save2,~,~,~] = load_dataset(option,16,folder_loc,'inter');
% if num == 16
%     if strcmp(destination,'inter')
%         extra = load([folder_loc,'\Ranking_Correlations\Data\extend_attentive_Krysiek.mat'],'nam_save','ranked'); % inter-plate sharing score
%     end
% elseif num == 12
%     if strcmp(destination,'inter')
%         extra = load([folder_loc,'\Ranking_Correlations\Data\extend_SS_inter_Krysiek.mat'],'nam_save','ranked');
%     end
% end
% nam_save2 = [nam_save2,extra.nam_save];
% [saved_swipes,~] = swipe_analysis(16,file_loc,nam_save2,'inter');
% % 
% % [saved_swipes,~] = swipe_analysis(num,file_loc,nam_save,'plates');
% % saved_swipes = saved_swipes + saved_swipes2;
% %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [folder_loc,alt_folder_loc,file_loc,floc] = setup()
    folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
    alt_folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data';
    file_loc = [folder_loc,'\adjs\adj_foodpile\']; % should match zone type
    floc=[alt_folder_loc,'\IQ_severity'];

    folder1 = [folder_loc,'\Set_allocate'];
    folder2 = [folder_loc,'\Plots'];
    folder3 = [folder_loc,'\Create_adj'];
    folder4 = [folder_loc,'\adjs\adj_foodpile'];
    folder5 = [folder_loc,'\Data'];
    folder6 = folder_loc;
    addpath(folder1,folder2,folder3,folder4,folder5,folder6)
end

function [nam_save,saved,ranked,list] = load_dataset(option,num,folder_loc,destination,bweight)
    saved = []; list = [];
%     if num == 16
%         % if strcmp(destination,'inter')
%         %     load([folder_loc,'\Ranking_Correlations\Data\extend_even_inter.mat'],'nam_save','ranked') % evenness
%         % else
%         %     load([folder_loc,'\Ranking_Correlations\Data\extend_SS_accurate.mat'],'nam_save','ranked')
%         % end
%     else
    if num == 12 || num == 16
        if strcmp(destination,'inter')
            load([folder_loc,'\Ranking_Correlations\Data\SS_ext_inter',bweight,'.mat'],'nam_save','ranked')
        else
            load([folder_loc,'\Ranking_Correlations\Data\SS_ext',bweight,'.mat'],'nam_save','ranked')
        end
%         if strcmp(destination,'inter')
%             load([folder_loc,'\Ranking_Correlations\Archive\extend_SS_inter.mat'],'nam_save','ranked') % inter-plate sharing score
%         else
%             load([folder_loc,'\Ranking_Correlations\Archive\extend_SS.mat'],'nam_save','ranked')
%         end
    end
    
%     end
end

function [nam_save,ranked,subject_details] = load_extra_data(nam_save,ranked,combine,option,num,folder_loc,destination,bweight)
    if combine == 0
        load('subject_details.mat')
        subject_details = subject_details_776;
    elseif combine == 1
        addpath 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_foodpile'
        load('subject_details_combine_ond.mat')
        subject_details = subject_details_combine;
    %     if num == 16
    %         % if strcmp(destination,'inter')
    %         %     extra = load([folder_loc,'\Ranking_Correlations\Data\extend_even_inter_Krysiek.mat'],'nam_save','ranked'); % evenness
    %         % else
    %         %     extra = load([folder_loc,'\Ranking_Correlations\Data\extend_SS_accurate_Krysiek.mat'],'nam_save','ranked');
    %         % end
    %     else
        if num == 12 || num == 16
            if strcmp(destination,'inter')
                extra = load([folder_loc,'\Ranking_Correlations\Data\SS_ext_inter_Krysiek',bweight,'.mat'],'nam_save','ranked');
            else
                extra = load([folder_loc,'\Ranking_Correlations\Data\SS_ext_Krysiek',bweight,'.mat'],'nam_save','ranked');
            end
%             if strcmp(destination,'inter')
%                 extra = load([folder_loc,'\Ranking_Correlations\Archive\extend_SS_inter_Krysiek.mat'],'nam_save','ranked'); % inter-plate sharing score
%             else
%                 extra = load([folder_loc,'\Ranking_Correlations\Archive\extend_SS_Krysiek.mat'],'nam_save','ranked');
%             end
        end
        ranked = [ranked;extra.ranked];
        nam_save = [nam_save,extra.nam_save];
    end
end

function [sets,months] = create_sets_months(subject_details,nam_save,saved,sex,tab_sev,severity)
    [months] = list_AGE(subject_details,nam_save,saved);
    if size(months,1)<size(months,2)
        months=months';
    end
    if strcmp(severity,'on')
        [sets] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,sex);
    elseif strcmp(sex,'Male')
        [tmp_sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
        sets = tmp_sets(5:8);
    elseif strcmp(sex,'Female')
        [tmp_sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
        sets = tmp_sets(1:4);
    elseif strcmp(sex,'compare')
        [sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
    else
        [sets] = set_allocate(subject_details,nam_save,saved);
%         [sets_sev] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,sex);
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
        keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
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
    
        if isfile([file_loc,file_id]) || isfile(['C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_foodpile\\',file_id])
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
            if num == 12
                n_swipes(jj) = sum(adj(4,[5,6,7]))+sum(adj(5,[4,6,7]))+sum(adj(6,[4,5,7]))+sum(adj(7,[4,5,6]));
            elseif num == 16
                n_swipes(jj) = sum(adj(4,[5,6,7,14,15,16]))+sum(adj(5,[4,6,7,13,15,16]))+sum(adj(6,[4,5,7,13,14,16]))+sum(adj(7,[4,5,6,13,14,15]));
            end
        elseif strcmp(destination,'total')
            n_swipes(jj) = sum(adj(:));
        elseif strcmp(destination,'evenness')
%             V1 = sum(adj(2,[4,13]));
%             V2 = sum(adj(2,[5,14]));
%             V3 = sum(adj(2,[6,15]));
%             V4 = sum(adj(2,[7,16]));
%             n_swipes(jj) = std([V1,V2,V3,V4])/sum([V1,V2,V3,V4]); 
            
            n_swipes(jj) = std(adj(2,[4,5,6,7]))/sum(adj(2,[4,5,6,7]));
% 
%             V1 = sum(adj(5,[4,13]))+sum(adj(6,[4,13]))+sum(adj(7,[4,13])) + sum(adj(2,[4,13]));
%             V2 = sum(adj(4,[5,14]))+sum(adj(6,[5,14]))+sum(adj(7,[5,14])) + sum(adj(2,[5,14]));
%             V3 = sum(adj(4,[6,15]))+sum(adj(5,[6,15]))+sum(adj(7,[6,15])) + sum(adj(2,[6,15]));
%             V4 = sum(adj(4,[7,16]))+sum(adj(5,[7,16]))+sum(adj(6,[7,16])) + sum(adj(2,[7,16]));
% 
%             n_swipes(jj) = std([V1,V2,V3,V4])/sum([V1,V2,V3,V4]);

            if isnan(n_swipes(jj))
                n_swipes(jj) = .11;
            end
        elseif strcmp(destination,'attentive')
            adj = adj - diag(diag(adj));
            n_swipes(jj) = sum(adj(2,[4,5,6,7]))/sum(adj(:));
            
%             V1 = sum(adj(5,[4,13]))+sum(adj(6,[4,13]))+sum(adj(7,[4,13]));% + sum(adj(2,[4,13]));
%             V2 = sum(adj(4,[5,14]))+sum(adj(6,[5,14]))+sum(adj(7,[5,14]));% + sum(adj(2,[5,14]));
%             V3 = sum(adj(4,[6,15]))+sum(adj(5,[6,15]))+sum(adj(7,[6,15]));% + sum(adj(2,[6,15]));
%             V4 = sum(adj(4,[7,16]))+sum(adj(5,[7,16]))+sum(adj(6,[7,16]));% + sum(adj(2,[7,16]));
%             n_swipes(jj) = sum([V1,V2,V3,V4])/sum(adj(:));

            if isnan(n_swipes(jj))
                n_swipes(jj) = 0;
            end
        elseif strcmp(destination,'diag')
            n_swipes(jj) = sum(diag(adj));%/(sum(adj(:)));
%             n_swipes(jj) = sum(adj(4,4))+sum(adj(5,5))+sum(adj(6,6))+sum(adj(7,7));
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

%     %% remove zn 4-7 incoming except from 2
%     allow=[2];%,4,5,6,7];
%     for it = 1:num
%         if ~ismember(it,allow)
%             adj(it,13:16) = adj(it,4:7)+adj(it,13:16);
%             adj(it,4:7)=zeros(1,4);                     % remove non-food connections
%         end
%     end
    %% remove zn 4-7 incoming except from 2
    allow=[2];
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

function [] = plot_results(results,months,sets,id,option,destination,num,sex,sev_choice,n_swipes)
    figure;
    x=months(sets{id});    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y=results(sets{id});
    % z=n_swipes(sets{id})';
%     x=z;
%     x(z==0)=[];
%     y(z==0)=[];
    [pf,S] = polyfit(x,y,2);
    % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
    [y_fit,delta] = polyval(pf,x,S);
    mean(delta)
    % Plot the original data, linear fit, and 95% prediction interval y�2?.
    [~,Ind]=sort(x,'asc');

    [r,p] = corr(x,y,'Type','Spearman');
    
%     %%% Export data and header
%     header = {'age', 'sharing', 'n_swipes'};
%     vars1 = [x y n_swipes(sets{id})'];
%     exportData = [header; num2cell(vars1)];
%     writecell(exportData, ['vars',num2str(id),'.csv']);

%     type = {'TD','ASD','OND'};
%     numopt = nchoosek(1:3,2);
%     for i = 1 : length(numopt)
%         x = results(sets{numopt(i,1)});
%         y = results(sets{numopt(i,2)});
%         header = {type{numopt(i,1)}, type{numopt(i,2)}, 'age'};
%         vars1 = [x y months(sets{id})];
%         exportData = [header; num2cell(vars1)];
%         writecell(exportData, ['vars',num2str(id),'.csv']);
%     end
% 
%     header = {'age', 'sharing', 'n_swipes'};
%     vars1 = [x y n_swipes(sets{id})'];
%     exportData = [header; num2cell(vars1)];
%     writecell(exportData, ['vars',num2str(id),'.csv']);

%     vars = [x y];
%     [r,p] = partialcorr(vars, n_swipes(sets{id})','Type','Spearman');
%     p = p(2);
%     r = r(2);
    
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
        text(0.05,.97+.025,['p = ',sprintf('%1.1e', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,['r = ',sprintf('%1.2f', r)],'Units','normalized','fontsize', 11)
        stars_only(3)
    elseif p<0.01
        text(0.05,.97+.025,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,['r = ',sprintf('%1.2f', r)],'Units','normalized','fontsize', 11)
        stars_only(2)
    elseif p<0.05
        text(0.05,.97+.025,['p = ',sprintf('%1.3f', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,['r = ',sprintf('%1.2f', r)],'Units','normalized','fontsize', 11)
        stars_only(1)
    elseif p<10%0.1
        text(0.05,.97+.025,['p = ',sprintf('%1.3f', p)],'Units','normalized','fontsize', 11)
        text(0.06,.91+.025,['r = ',sprintf('%1.2f', r)],'Units','normalized','fontsize', 11)
    end
    
    if option == 2
%         axis([28 74 -.5 .25])
        axis([28 74 -.65 .3])
    elseif option == 1
        if strcmp(destination,'inter') 
        	axis([28 74 0 60])
        elseif strcmp(destination,'food') 
            axis([28 74 0 164])
        elseif strcmp(destination,'total')
            axis([28 74 0 650])
        else
            axis([28 74 0 147])
        end
    elseif option == 3
        axis([28 74 0 1])
    end

    xlabel('Age (months)','fontsize', 11)
    if option == 1
        if strcmp(destination,'plates')
            if num == 16
                ylabel('No. of food-to-plate swipes','fontsize',14)
%                 ylabel('No. of swipes (all)','fontsize',14)
                nam = 'n_swipes_16_';
            elseif num == 12
                ylabel('No. of food delivery swipes','fontsize',14)
%                 ylabel('No. of swipes (inter-plates)','fontsize',14)
                nam = 'n_swipes_12_';
            end
        elseif strcmp(destination,'food')
            ylabel('No. of zone 2 only swipes ','fontsize',14)
            nam = 'n_swipes_food_';
        elseif strcmp(destination,'inter')
            ylabel('No. of inter-plate swipes','fontsize',14)
            nam = 'n_swipes_inter_';
        elseif strcmp(destination,'evenness')
            ylabel('Standard deviation','fontsize',14)
            nam = 'n_swipes_even_';
        elseif strcmp(destination,'total')
            ylabel('Total no. of swipes','fontsize',14)
            nam = 'n_swipes_total_';
        else
            nam = 'other';
        end
    elseif option == 2
        if num == 16
            ylabel('Direct sharing score','fontsize',14)
%             ylabel('Sharing score (plates)','fontsize',14)
            nam = 'sharing_16_';
        elseif num == 12
            ylabel('Direct sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Indirect sharing score','fontsize',14)
            end
            nam = 'sharing_12_';
        end
    elseif option == 3
        ylabel('Swipe accuracy ratio','fontsize',14)
        nam = 'ratio_';
    end
%     legend('subject','2nd order fit','Location','SouthEast','Orientation','horizontal','fontsize',14)
    
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

    f=gcf;
    f.Position = [403,340,330,313];

    if strcmp(destination,'food')
        saveas(gcf,['Figures/food_',nam,titletxt,'_',sex,'.png'])
    elseif sev_choice>0
        saveas(gcf,['Figures/severity_',nam,titletxt,'_',sex,'.png'])
    elseif strcmp(destination,'inter')
        saveas(gcf,['Figures/',nam,titletxt,'_',sex,'_inter.png'])
    else
        saveas(gcf,['Figures/',nam,titletxt,'_',sex,'_plates.png'])
    end
end

function [] = stars_line(n_stars,height,strt,nd,drp2)
%     drp2 = 3.5*drp/2;
    drp = 2*drp2/3.5;
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

function [] = Plot_boxplots(all_sets,grps,option,destination,num,sex,severity,save_p,combos)
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
    
    if strcmp(severity,'on')
        xticklabels({'WP','ASD 1','ASD 2','ASD 3'});
    else
        xticklabels({'WP','ASD','OND'});
    end
    
    if option == 1 
        if strcmp(destination,'plates')
            if num == 16
                ylabel('No. of swipes (plates)','fontsize',14)
%                 ylabel('No. of swipes (all)','fontsize',14)
            elseif num == 12
                ylabel('No. of food delivery swipes','fontsize',14)
%                 ylabel('No. of swipes (inter-plates)','fontsize',14)
            end
        elseif strcmp(destination,'food')
            ylabel('No. of swipes (zone 2 only)','fontsize',14)
        elseif strcmp(destination,'inter')
            ylabel('No. of swipes (inter-plate)','fontsize',14)
        elseif strcmp(destination,'evenness')
            ylabel('Standard deviation','fontsize',14)
        end
    elseif option == 2
        if num == 16
            ylabel('Direct sharing score','fontsize',14)
        elseif num == 12
            ylabel('Direct sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Indirect sharing score','fontsize',14)
            end
        end
    elseif option == 3
        ylabel('Swipe accuracy ratio','fontsize',14)
    end
    
    box off

    if strcmp(destination,'food')
        axis([0.5 3.5 -3 164])
    elseif strcmp(destination,'inter') && option == 1
        axis([0.5 3.5 0 60])
    elseif strcmp(destination,'total')
        axis([0.5 3.5 0 650])
    elseif option == 2
        if strcmp(severity,'on')
            axis([0.5 4.5 -.64 .3])
        else
            axis([0.5 3.5 -.64 .3])
        end
    elseif option == 1
        if num == 16
            if strcmp(severity,'on')
                axis([0.5 4.5 0 147])
            else
                axis([0.5 3.5 0 147])
            end
        elseif num == 12
            if strcmp(severity,'on')
                axis([0.5 4.5 0 147])
            else
                axis([0.5 3.5 0 147])
            end
        end
    elseif option == 3
        axis([0.5 3.5 0 1.2])
    end

%     %% Plot significance stars
%     if strcmp(destination,'food')
%         heights = [135,143,151,130];
%         drp = 3;
%     elseif strcmp(severity,'on')
%         if option == 2
%             heights = [.4,.35,.3];
%             drp = .0225;
%         elseif option == 1
% %             heights = [125,120,115];
%             heights = [155,147.5,140];
%             drp = 3;
%         end
%     elseif option == 2
%         heights = [.28,.33,.39,.28];%,.25];
% %         heights = [.41,.345,.45,.31];
%         drp = .02;
%     elseif option == 1
%         if num == 16
%             heights = [122,130,140,122];
%         elseif num == 12
%             heights = [155,162,171,155];
%         end
%         drp = 3;
%     elseif option == 3
%         heights = [1.05,1.105,1.17,1.05];
%         drp = 0.025;
%     end
%  
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

    title(sex)

    %% Plot significance stars    
    auto_star_plot(ic,all_sets,save_p,combos,option)

    if strcmp(destination,'food')
        saveas(gcf,['Figures/boxplot_food_',sex,'.png'])
    elseif strcmp(destination,'inter')
        saveas(gcf,['Figures/boxplot_option',num2str(option),'_inter_',num2str(num),'_',sex,'.png'])     
    elseif strcmp(severity,'on')
        saveas(gcf,['Figures/boxplot_severity_option',num2str(option),'_',num2str(num),'_',sex,'.png'])
    else
        saveas(gcf,['Figures/boxplot_option',num2str(option),'_plates_',num2str(num),'_',sex,'.png'])
    end
end

function [] = auto_star_plot(ic,all_sets,save_p,combos,option)
    mv = zeros(1,4);
    for i = 1:4
        Ind=find(ic==i);
        mv(i) = max(all_sets(Ind));
    end

    if option == 2
        drp = max(mv)*0.175;
    elseif option == 1
        drp = max(mv)*0.05;
    end
    addh=drp*1.5;
    h_save=0;
    for i = 1 : length(save_p)
        j = length(save_p)+1-i;
%         j=i;
        id1 = combos(j,1);
        id2 = combos(j,2);
        height = max([mv(id1:id2)])+addh;
        height = max([h_save+addh,height]);
        if save_p(j) < 0.001
            stars_line(3,height,id1,id2,drp)
        elseif save_p(j) < 0.01
            stars_line(2,height,id1,id2,drp)
        elseif save_p(j) < 0.05
            stars_line(1,height,id1,id2,drp)
        end
        if save_p(j) < 0.05
            h_save = height;
        end
    end
end

function [save_p,combos] = significance_check(results,sets,n)
    % pairwise significance check
    combos=nchoosek(1:n,2);
    save_p = zeros(1,size(combos,1));
    for j = 1 : size(combos,1)
        num = combos(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval_kw = kruskalwallis([results(sets{num(1)})';results(sets{num(2)})'],len_rankeds,'off');

        % Extract the data for the two groups
        group1 = results(sets{num(1)});
        group2 = results(sets{num(2)});
        
        % Perform the two-tailed unpaired Wilcoxon test
        [pval, h] = ranksum(group1, group2);
        save_p(1,j) = pval;
    end
end

function [save_p] = significance_sex(results,sets)
    save_p = zeros(1,4);
%     numopt = [1,5;2,6;3,7];
    numopt = nchoosek([1,5,2,6,3,7],2);
    for j = 1 : length(numopt)
        num=numopt(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval_kw = kruskalwallis([results(sets{num(1)});results(sets{num(2)})],len_rankeds,'off');

        % Extract the data for the two groups
        group1 = results(sets{num(1)});
        group2 = results(sets{num(2)});
        
        % Perform the two-tailed unpaired Wilcoxon test
        [pval, h] = ranksum(group1, group2);

        save_p(1,j) = pval;
%         if pval > .001
%             text(0.05+(j-1)*1/4,.92,['p = ',sprintf('%1.4f', pval)],'Units','normalized','fontsize', 10)
%         else
%             text(0.05+(j-1)*1/4,.92,['p = ',sprintf('%1.1e', pval)],'Units','normalized','fontsize', 10)
%         end
    end
end

% function [save_p] = significance_sex_separate(results,sets)
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
function [ ] = boxplot_sex_cmpr(all_sets,grps,option,destination,num)

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
            axis([0 9 -5 147])
        elseif num == 12
            axis([0 9 -5 147])
        end
    elseif option == 2
        axis([0 9 -.64 .3])
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
    xlabel('Sex')

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
        elseif strcmp(destination,'evenness')
            ylabel('Standard deviation','fontsize',14)
        end
    elseif option == 2
        if num == 16
            ylabel('Direct sharing score','fontsize',14)
        elseif num == 12
            ylabel('Direct sharing score','fontsize',14)
            if strcmp(destination,'inter')
                ylabel('Indirect sharing score','fontsize',14)
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

%     if strcmp(destination,'inter')
%         saveas(gcf,['Figures/compare_',num2str(num),'_',textadd,'_inter.png'])
%     elseif strcmp(destination,'food')
%         saveas(gcf,['Figures/compare_',num2str(num),'_',textadd,'_food.png'])
%     else
%         saveas(gcf,['Figures/compare_',num2str(num),'_',textadd,'.png'])
%     end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
