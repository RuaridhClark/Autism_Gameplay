clear all
for s = 1 : 11
    
    load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')
%     load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')

    fileloc='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\IQ_severity';
    tab_sev = readtable([fileloc,'\eCRF.csv']);

    % folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end';
    folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
    % folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
    folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj_110721';
    folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
    addpath(folder4,folder6,folder7)% addpath(folder3,folder4,folder5,folder6,folder7)
    file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end\'; % should match zone type

    load('swipes_all704.mat','nam_save')

    %% stack the adjs
    num =12;    % number of ipad objects (nodes)

    asdname_save = cell(704,1);
    set = [];
    iqs=[];
    val=[];
    rmv_m=[];
    subj_grp = 'ASD';


    [all_scores,titlename] = iq_test(s,tab_sev);
    
    
    m = 0;
    map=zeros(1,height(tab_sev));
    iqs=zeros(1,height(tab_sev));
    for i = 1:height(tab_sev)
        skip=1;
        file_id = ['subject_',nam_save{i},'.mat'];
        
        iq_score = all_scores(i);
        if s==1 && isempty(iq_score) 
            iq_score = tab_sev.GriffithsIQ(i);
            if isnan(iq_score)
                iq_score = tab_sev.MerrillPalmerIQ(i);
            end
        end    
        if ~isempty(iq_score) && ~isnan(iq_score)
            % map table entry to 704 participants
            tmp = find(strcmp(nam_save,tab_sev.id_study_id{i}));
            if ~isempty(tmp)
                map(i)= tmp;
            end
         
                tmp_iq = iq_score;
            
                nanfound=0;
                if ~isnan(tmp_iq)
                    set = [set,i];
                    iqs(i)=tmp_iq;
                    try     % Check for adjacency matrix
                        load([file_loc,file_id])
                        val(i) = sum(adj(2,[4,5,6,7]));
                    catch
                        val(i)=-1;
                    end
                else
                    nanfound=1;
                end
                if sum(iqs>0)~=length(set)
                    set(end)=[];
                    rmv_m = [rmv_m,m];
                elseif nanfound == 0 
                    asdname_save{i} = tab_sev.id_study_id{i}; 
                else
                    rmv_m = [rmv_m,m];
                end
            
            savename = ['subject_',nam_save{i}];
        end
    end
    
    load('H:\My Documents\MATLAB\Autism_MAIN\subject_details.mat')
    [months] = list_AGE_update(subject_details_776,asdname_save);
    months=[months,0]; % hacky
    val=[val,0];
    
    rmv = find(map==0);
    months(rmv)=[];
    iqs(rmv)=[];
    val(rmv)=[];
    map(rmv)=[];
    
    rmv = find(iqs==0);
    months(rmv)=[];
    map(rmv)=[];
    iqs(rmv)=[];
    val(rmv)=[];
    
    saved=saved(map);
    
    rmv = find(sum(saved,1)==0);
    map_pert = map;
    map_pert(rmv)=[];
    iqs_pert = iqs;
    iqs_pert(rmv)=[];
    ranked=ranked(map_pert);

    iqs_pert(val==-1)=[];
    ranked(val==-1)=[];
    val(val==-1)=[];
    
    iqs_pert(ranked==0.5)=[];
    val(ranked==0.5)=[];
    ranked(ranked==0.5)=[];

    if size(months,1)<size(months,2)
        months=months';
    end
    
    pert_init=-0.3;
    pert_chng = 0.01;
    if min(ranked)>=0   % change ranked to match pert ** NEED TO EXCLUDE 0.5s
        ranked = (ranked-1)*pert_chng+pert_init;
    end

%     %% Correlation
% 
%     if min(ranked)>=0   % change ranked to match pert
%         ranked = (ranked-1)*pert_chng+pert_init;
%     end
% 
%     %% Correlation scatter plot
%     for num = 2
%         figure;
%         
%         x=months;%iqs';%
%         y=iqs';%ranked;%
%         [pf,S] = polyfit(x,y,1);
%         % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
%         [y_fit,delta] = polyval(pf,x,S);
%         % Plot the original data, linear fit, and 95% prediction interval y±2?.
%         [~,Ind]=sort(x,'asc');
%         [R,p] = corr(x,y,'Type','Kendall');
%         if p<0.01
%             linetype='-';
%         elseif p<0.05
%             linetype='--';
%         else
%             linetype=':';
%         end
%         if num == 1
%             clr = [0, 0.4470, 0.7410];
%         elseif num == 2
%             clr = [0.8500, 0.3250, 0.0980];
%         elseif num == 3
%             clr = [0.9290, 0.6940, 0.1250];
%         elseif num == 4
%             clr = [0.4940, 0.1840, 0.5560];
%         end
%         tempclr=[.5 .5 .5];
%     %     scatplot=plot(x(Ind),y(Ind),'x','color',tempclr,'linewidth',1);
%         scatplot=scatter(x(Ind),y(Ind),55,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
%         % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
%         scatplot.MarkerFaceAlpha = .3;
%     %     scatplot.MarkerEdgeAlpha = 0;
%         hold on
%         plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)
% 
%     %     [rho,pval] = corr(months(sets{num}),ranked(sets{num}),'Type','Kendall');
%     %     text(55,mean(ranked(sets{num})),['p_{Ken,\tau} = ',num2str(p)])
%         text(0.05,.95,['p_{Ken,\tau} = ',num2str(p)],'Units','normalized')
%     %     axis([27 73 -20 25])
%         xlabel('Age (months)')%xlabel('IQ')%
%         ylabel('IQ')%ylabel('Pertubation threshold')%
%         legend('subject','linear fit','Location','SouthEast')
% 
%         title(titlename)
%         if p>0.05
%             close gcf
%         end
%     end

    for num = 2
        figure;
        
        x=iqs_pert';%months;%
        y=val';%ranked;%iqs';%
        [pf,S] = polyfit(x,y,2);
        % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
        [y_fit,delta] = polyval(pf,x,S);
        % Plot the original data, linear fit, and 95% prediction interval y±2?.
        [~,Ind]=sort(x,'asc');
        [R,p] = corr(x,y,'Type','Kendall');
        if p<0.01
            linetype='-';
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
        tempclr=[.5 .5 .5];
    %     scatplot=plot(x(Ind),y(Ind),'x','color',tempclr,'linewidth',1);
        scatplot=scatter(x(Ind),y(Ind),55,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
        % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
        scatplot.MarkerFaceAlpha = .3;
    %     scatplot.MarkerEdgeAlpha = 0;
        hold on
        plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)

    %     [rho,pval] = corr(months(sets{num}),ranked(sets{num}),'Type','Kendall');
    %     text(55,mean(ranked(sets{num})),['p_{Ken,\tau} = ',num2str(p)])
        text(0.05,.95,['p_{Ken,\tau} = ',num2str(p)],'Units','normalized')
    %     axis([27 73 -20 25])
        xlabel('IQ')%xlabel('Age (months)')%
        ylabel('number of swipes ending in zones 4-7')%ylabel('IQ')%ylabel('Pertubation threshold')%
        legend('subject','2nd order fit','Location','SouthEast')

        title(titlename)
%         if p>0.05
%             close gcf
%         end
    end
end
%%%%%%%%% functions %%%%%%%%%%%%%%
function [all_scores,titlename] = iq_test(s,tab_sev)
    if s == 1
        all_scores = tab_sev.iq_wisc__general;
        all_scores = str2double(all_scores);%cell
        titlename = 'iq_wisc__general';
    elseif s == 2
        all_scores = tab_sev.iq_wisc__verbal_comprehension_index;
        titlename = 'iq_wisc__verbal_comprehension_index';
    elseif s == 3
        all_scores = tab_sev.iq_wisc__visual_spatial_index;
        titlename = 'iq_wisc__visual_spatial_index';
    elseif s == 4
        all_scores = tab_sev.iq_wisc__fluid_reasoning_index;
        titlename = 'iq_wisc__fluid_reasoning_index';
    elseif s == 5
        all_scores = tab_sev.iq_wisc__working_memory_index;
        titlename = 'iq_wisc__working_memory_index';
    elseif s == 6
        all_scores = tab_sev.iq_wisc__processing_speed_index;
        titlename = 'iq_wisc__processing_speed_index';
    elseif s == 7
        all_scores = tab_sev.vabs__overall_score;
        titlename = 'vabs__overall_score';
    elseif s == 8
        all_scores = tab_sev.vabs__communication;
        titlename = 'vabs__communication';
    elseif s == 9
        all_scores = tab_sev.vabs__daily_living_skills;
        titlename = 'vabs__daily_living_skills';
    elseif s == 10
        all_scores = tab_sev.vabs__socialization;
        titlename = 'vabs__socialization';
    elseif s == 11
        all_scores = tab_sev.vabs__motor_skills;
        titlename = 'vabs__motor_skills';
    elseif s == 12
        all_scores = tab_sev.vabs__maladaptive_behavior;
        all_scores = str2double(all_scores);%cell
        titlename = 'vabs__maladaptive_behavior';
    end
        %%%
%%%%%%

end