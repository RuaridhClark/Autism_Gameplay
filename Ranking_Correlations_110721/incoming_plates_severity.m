% % % check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate.mat')

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

%% stack the adjs
num =16;    % number of ipad objects (nodes)
saved = zeros(num,704);
save_V = zeros(num,704);

n_swipes = zeros(1,704);
f_num = 0;

diagsA=zeros(12,704);
list=[];
for jj = 1:704
    skip=1;
    file_id = ['subject_',nam_save{jj},'.mat'];

    if isfile([file_loc,file_id])
        f_num = f_num + 1;
        load(file_id)
        adj=adj(1:num,1:num);
        titlename = ['ID ',nam_save{jj}];
        savename = ['subject_',nam_save{jj}];
        R = f_num;
    else 
        skip=0;
    end

    if max(adj(:))==0
        list=[list,jj];
    end

    
    n_swipes(jj) = sum(adj(2,[4,5,6,7]));
    adj=adj-diag(diag(adj));
end

%% Correlation
load('diags_incoming.mat')

load('subject_details.mat')
saved=ones(1,704);
saved(list)=zeros(1,length(list));
[months] = list_AGE(subject_details_776,nam_save,saved);

if size(months,1)<size(months,2)
    months=months';
end

exclude = find(ranked==0.5);
exclude = [exclude;find(months>75)]; %% 75 months threshold
% exclude = [exclude;find(months<32)]; %% 45 Months threshold
for i = 1 : length(sets)
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
end

iter=0;all_sets=[];grps=[];
for num = 1:3
%         h = kstest(diagsA(zones(2),sets{num}))
    iter=iter+1;
%         vals = diagsA(zones(zns),sets{num});
    vals = n_swipes(1,sets{num});
    all_sets = [all_sets,vals];%val(sets{num})];
    grps = [grps,ones(1,length(vals))*iter];%val(sets{num})))*iter];
end

for num = 1 : 3
    figure;
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=months(sets{num});
    y=n_swipes(sets{num})';
    [pf,S] = polyfit(x,y,2);
    % Evaluate the first-degree polynomial fit in p at the points in x. Specify the error estimation structure as the third input so that polyval calculates an estimate of the standard error. The standard error estimate is returned in delta.
    [y_fit,delta] = polyval(pf,x,S);
    % Plot the original data, linear fit, and 95% prediction interval y±2?.
    [~,Ind]=sort(x,'asc');
    [R,p] = corr(x,y,'Type','Spearman');
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
        clr = [.85,.33,.1];
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
    plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)

%     [rho,pval] = corr(months(sets{num}),ranked(sets{num}),'Type','Kendall');
%     text(55,mean(ranked(sets{num})),['p_{Ken,\tau} = ',num2str(p)])
    text(0.05,.95,['p_{Ken,\tau} = ',num2str(p)],'Units','normalized')
    
    axis([27 73 0 110])

    xlabel('Age (months)')
    ylabel('No. of swipes - food to plate')
    legend('subject','2nd order fit','Location','NorthWest')
    
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
    vals = n_swipes(sets{num});
    all_sets = [all_sets,vals];
    grps = [grps,ones(1,length(vals))*iter];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=figure;
boxplot(all_sets,grps,'Notch','on')
h = findobj(gca,'Tag','Box');
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
for m=1:length(h)
    temp_m = length(h)-m+1;
    mm=rem(m,5);mm(mm==0)=1;
    patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
end

% text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
% xticklabels({'TD','ASD','OND*','ONDE'});
% set(gca,'xticklabel',entries,'fontsize',10)

ylabel('No. of swipes - food to plate')

axis([0.5 3.5 -3 150])

% legend('TD','ASD','OND','ONDE')
box off

% %% Plot significance stars
% if option == 1
%     heights = [.29,.33,.39,.29];
%     drp = .0175;
% elseif option == 2
%     heights = [28,32,37,28];
%     drp = 1.75;
% end
% % height=28;%0.29; 
% n_stars = 3; 
% stars_line(n_stars,heights(1),1,2,drp) % 3 stars,h,1,2,1
% % height = 32;%.33;
% stars_line(n_stars,heights(2),1,3,drp) % 3 stars,h,1,2,1
% % height = 37;%.39; 
% stars_line(n_stars,heights(3),2,4,drp) % 2 stars,h,1,3,2
% % height = 28;%.29;
% stars_line(n_stars,heights(4),3,4,drp) % 2 stars,h,3,4,1

combos=nchoosek([1,2,3],2);
save_p = zeros(1,size(combos,1));
for j = 1 : size(combos,1)
    num = combos(j,:);
    len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
    pval = kruskalwallis([n_swipes(sets{num(1)})';n_swipes(sets{num(2)})'],len_rankeds,'off');
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
