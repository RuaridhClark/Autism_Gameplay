% % % check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
option = 1; % 1 == proportional, 2 == proportion + swipe volume
if option == 1
    load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate_proport_bi.mat')
%     load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')
elseif option == 2
    load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate_bi.mat')
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
% file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end\'; % should match zone type

%%
fileloc='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\IQ_severity';
tab_sev = readtable([fileloc,'\eCRF.csv']);

asdname_save = cell(704,1);
all_severity = zeros(704,1);
sets = cell(3,1);
check = 0;

m = 0;
map=zeros(1,height(tab_sev));
file_loc = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate\';
for i = 1:height(tab_sev)
    skip=1;
    file_id = ['subject_',tab_sev.id_study_id{i},'.mat'];

    severity = tab_sev.clinical_diagnosis__asd_severity_level{i};
    if isfile([file_loc,file_id]) && ~isempty(severity)
%         m=m+1;
        map(i)=find(strcmp(nam_save,tab_sev.id_study_id{i}));
        conf = tab_sev.clinical_diagnosis__asd_severity_level_rating_confidence{i};
%                 if strcmp(conf,'High') || strcmp(conf,'Medium') || strcmp(conf,'Low') 
        [sev_num] = severity_score(severity); 
        sets{sev_num}=[sets{sev_num},map(i)];

    end
end

%% Correlation
load('diags_incoming.mat')

load('subject_details.mat')
[months] = list_AGE(subject_details_776,nam_save,saved);

if size(months,1)<size(months,2)
    months=months';
end

val = -.70;
if option == 1
    ranked(ranked<val)=val;
end

exclude = find(ranked==0.5);
% exclude = [exclude;find(ranked<val)];
exclude = [exclude;find(months>75)]; %% 75 months threshold
% exclude = [exclude;489]; % no food-to-plate swipes recorded
% exclude = [exclude;find(months<32)]; %% 45 Months threshold
for i = 1 : length(sets)
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
end

if min(ranked)>=0   % change ranked to match pert ** NEED TO EXCLUDE 0.5s
    try
        ranked = (ranked-1)*pert_chng+pert_init;
    catch
        ranked = ranked;%+pert_init;
    end
end

for num = 1 : 3
    figure;    
    x=months(sets{num});
    y=ranked(sets{num});
    [pf,S] = polyfit(x,y,2);
    % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
    [y_fit,delta] = polyval(pf,x,S);
    % Plot the original data, linear fit, and 95% prediction interval y±2?.
    [~,Ind]=sort(x,'asc');
    [R,p] = corr(x,y,'Type','Spearman');
%     [R,PValue,H] = corrplot([x,y],'type','Kendall','testR','on');
%     figure
    if p<0.01
        linetype='-';
    elseif p<0.05
        linetype='--';
    else
        linetype=':';
    end
    if num == 1
        clr = [.93,.69,.13];
    elseif num == 2
%         clr = [.85,.33,.1];
        clr = [.84,.35,.15];
    elseif num == 3
        clr = [.64,.08,.18];
    end
    tempclr=[.5 .5 .5];
%     scatplot=plot(x(Ind),y(Ind),'x','color',tempclr,'linewidth',1);
    scatplot=scatter(x(Ind),y(Ind),55,'o','MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatplot.MarkerFaceAlpha = .3;
%     scatplot.MarkerEdgeAlpha = 0;
    hold on
    if p<0.05
        plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)
    end

    if p<0.001
        text(0.05,.05,['p = ',sprintf('%1.1e', p)],'Units','normalized','fontsize', 11)
    elseif p<0.01
        text(0.05,.05,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    elseif p<10%0.1
        text(0.05,.05,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    end
    
    if option == 1
        axis([27 73 -.7 .31])
    elseif option == 2
       axis([27 73 -70 30])
    end

    xlabel('Age (months)')
    if option == 1
        ylabel('Sharing score (plates)')
    elseif option == 2
        ylabel('Sharing score (swipe volume)')
    end
    legend('subject','2nd order fit','Location','SouthEast')
    
    if num == 1
        title('1')
    elseif num == 2
        title('2')
    elseif num == 3
        title('3')
    end
end

%% Boxplot age - Food as an Origin
load('subject_details.mat')
saved=ones(1,704); saved(list)=zeros(1,length(list));

iter=0;all_sets=[];grps=[];
zones=1:17;zones(4)=[];

for num = 1:3
%         h = kstest(diagsA(zones(2),sets{num}))
    iter=iter+1;
    vals = ranked(sets{num})';
    all_sets = [all_sets,vals];
    grps = [grps,ones(1,length(vals))*iter];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=figure;
boxplot(all_sets,grps,'Notch','on')
h = findobj(gca,'Tag','Box');
% colors = [.93,.69,.13;.85,.33,.1;.64,.08,.18;];
colors = [.93,.69,.13;.84,.35,.15;.64,.08,.18;];
for m=1:length(h)
    temp_m = length(h)-m+1;
    mm=rem(m,5);mm(mm==0)=1;
    patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
end

% text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
% xticklabels({'TD','ASD','OND*','ONDE'});
% set(gca,'xticklabel',entries,'fontsize',10)

if option == 1
    ylabel('Sharing score (proportion)')
elseif option == 2
    ylabel('Sharing score (swipe volume)')
end

if option == 1
    axis([0.5 3.5 -.7 .39])
elseif option == 2
    axis([0.5 3.5 -51 39])
end

% legend('TD','ASD','OND','ONDE')
box off

%% Plot significance stars
if option == 1
    heights = [.29,.36];
%     heights = [.315,.38,.45];
    drp = .0275;
elseif option == 2
    heights = [28,32,37];
    drp = 1.75;
end
% 
% n_stars = 3; 
% stars_line(n_stars,heights(1),1,2,drp) % 3 stars,h,1,2,1
% stars_line(2,heights(2),1,3,drp) % 3 stars,h,1,2,1
% stars_line(1,heights(3),2,3,drp) % 2 stars,h,1,3,2

n_stars = 2; 
stars_line(n_stars,heights(1),1,2,drp) % 3 stars,h,1,2,1
stars_line(1,heights(2),1,3,drp) % 3 stars,h,1,2,1

combos=nchoosek([1,2,3],2);
save_p = zeros(1,size(combos,1));
for j = 1 : size(combos,1)
    num = combos(j,:);
    len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
    pval = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds,'off');
    save_p(1,j) = pval;
end


% legend('TD','ASD','OND*','ONDE','Orientation','horizontal')

f.Position = [403,340,300,313];

%%%%%% function %%%%%%
function [sev_num] = severity_score(severity)
    if strcmp(severity,'1. Level 1 "Requiring support"')
        sev_num = 1;
    elseif strcmp(severity,' 2. Level 2 "Requiring substantial support"')
        sev_num = 2;
    elseif strcmp(severity,'3. Level 3 "Requiring very substantial support"')
        sev_num = 3;
    end
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
