clear all
% load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')

fileloc='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\IQ_severity';
tab_sev = readtable([fileloc,'\eCRF.csv']);

folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end';
folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj';
folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
addpath(folder3,folder4,folder5,folder6,folder7)
file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =12;    % number of ipad objects (nodes)

asdname_save = cell(704,1);
set = [];
iqs=[];
rmv_m=[];
subj_grp = 'ASD';

    m = 0;
    map=[];
    for i = 1:height(tab_sev)
        skip=1;
        file_id = ['subject_',nam_save{i},'.mat'];
        
        iq_score = tab_sev.iq_wisc__verbal_comprehension_index(i);%tab_sev.iq_wisc__general{i};%tab_sev.vabs__maladaptive_behavior{i};%tab_sev.vabs__motor_skills(i);%tab_sev.vabs__socialization(i);%tab_sev.vabs__daily_living_skills(i);%tab_sev.vabs__communication(i);%tab_sev.vabs__overall_score(i);%tab_sev.iq_wisc__processing_speed_index(i);%tab_sev.iq_wisc__working_memory_index(i);%tab_sev.iq_wisc__fluid_reasoning_index(i);%tab_sev.iq_wisc__visual_spatial_index(i);%
        if ~isempty(iq_score) %&& isfile([file_loc,file_id]) && strcmp(tab_sev.diagnosis_category{i},subj_grp) && max(strcmp(nam_save,tab_sev.id_study_id{i}))
            m = m + 1;
            map(m)=find(strcmp(nam_save,tab_sev.id_study_id{i}));
            load(file_id)
%                 iqs = [iqs,iq_score];%[iqs,str2num(iq_score)];%
         
                tmp_iq = iq_score;
                if ischar(iq_score)
                    tmp_iq = str2num(iq_score);
                end
            
                nanfound=0;
                if ~isnan(tmp_iq)
                    set = [set,i];

                    iqs(m)=tmp_iq;

                else
%                     iqs(end)=[];
                    nanfound=1;
                end
                if sum(iqs>0)~=length(set)
                    set(end)=[];
%                     asdname_save(end) = [];
                    rmv_m = [rmv_m,m];
                elseif nanfound == 0 
                    asdname_save{m} = nam_save{i}; 
                else
                    rmv_m = [rmv_m,m];
                end
            
            
            titlename = ['ID ',nam_save{i}];
            savename = ['subject_',nam_save{i}];

        end
    end
    
    ranked=ranked(map);

% save('temp_save_OBJ.mat')

%% Correlation
% load('diags_incoming.mat')
% 
load('subject_details.mat')
[months] = list_AGE_update(subject_details_776,asdname_save,saved);
iqs(months==0)=[];
ranked(months==0)=[];
months(months==0)=[]; % removing invalid asdname_save entries
if size(months,1)<size(months,2)
    months=months';
end

%% Correlation

if min(ranked)>=0   % change ranked to match pert
    ranked = (ranked-1)*pert_chng+pert_init;
end

% ranked(rmv_m)=[];

%% Correlation scatter plot
for num = 2
    figure;
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=months;%iqs';%
    y=iqs';%ranked;%
    [pf,S] = polyfit(x,y,1);
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
    xlabel('Age (months)')%xlabel('IQ')%
    ylabel('IQ')%ylabel('Pertubation threshold')%
    legend('subject','linear fit','Location','SouthEast')
    
    if num == 2
        title(subj_grp)
    end
end

for num = 2
    figure;
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=iqs';%months;%
    y=ranked;%iqs';%
    [pf,S] = polyfit(x,y,1);
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
    ylabel('Pertubation threshold')%ylabel('IQ')%
    legend('subject','linear fit','Location','SouthEast')
    
    if num == 2
        title(subj_grp)
    end
end