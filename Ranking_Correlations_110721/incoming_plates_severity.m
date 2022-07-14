% % % check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
load('C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate_bi.mat')

% folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
% folder2 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate\functions';
folder3 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate';
folder4 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Set_allocate';
folder5 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj_110721';
folder7 = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
addpath(folder3,folder4,folder5,folder6,folder7)
% file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\adj_obj_end\'; % should match zone type

%%
fileloc='C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\IQ_severity';
tab_sev = readtable([fileloc,'\eCRF.csv']);

asdname_save = cell(704,1);
all_severity = zeros(704,1);
sets = cell(3,1);
check = 0;

m = 0;
map=zeros(1,height(tab_sev));
file_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate\';
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

        if strcmp(tab_sev.patient_data__sex{i},'Female') %%%%%%%%%%%% TEMPORARY - REMOVE %%%%%%%%%%%%%%%
            sets{sev_num}=[sets{sev_num},map(i)];
        end
    end
end

%% stack the adjs
num =16;    % number of ipad objects (nodes)
saved = zeros(num,704);
save_V = zeros(num,704);

n_swipes = zeros(1,704);
f_num = 0;

% diagsA=zeros(12,704);
list=[];
for jj = 1:704
    skip=1;
    file_id = ['subject_',nam_save{jj},'.mat'];

    if isfile([file_loc,file_id])
        f_num = f_num + 1;
        load(file_id)
        if num == 12
            [adj] = adj_snap2zones(adj,num);
        end
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

% load('n_swipes_snapto.mat') %load('n_swipes_accurate.mat')

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
        text(0.05,.95,['p = ',sprintf('%1.1e', p)],'Units','normalized','fontsize', 11)
    elseif p<0.01
        text(0.05,.95,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    elseif p<10%0.1
        text(0.05,.95,['p = ',sprintf('%1.4f', p)],'Units','normalized','fontsize', 11)
    end
    
    axis([27 73 0 110])

    xlabel('Age (months)','fontsize', 11)
%     ylabel('No. of swipes - food to plates','fontsize', 11)
    ylabel('No. of swipes - food to {\it{snap-to-target}})','fontsize', 11)
%     legend('subject','2nd order fit','Location','NorthWest')
%     legend('subject','2nd order fit','Location','SouthEast')
    
    if num == 1
        title('1')
        legend('subject','2nd order fit','Location','SouthEast')
    elseif num == 2
        title('2')
        legend('subject','2nd order fit','Location','NorthWest')
    elseif num == 3
        title('3')
        legend('subject','2nd order fit','Location','NorthWest')
    end
    
    f=gcf;
    f.Position = [403,340,330,313];
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
colors = [.93,.69,.13;.84,.35,.15;.64,.08,.18;];
for m=1:length(h)
    temp_m = length(h)-m+1;
    mm=rem(m,5);mm(mm==0)=1;
    patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
end

% text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
% xticklabels({'TD','ASD','OND*','ONDE'});
% set(gca,'xticklabel',entries,'fontsize',10)

xlabel('ASD Severity')
ylabel('No. of swipes - food to plates')
% ylabel('No. of swipes - food to {\it snap-to-target} zone')

% axis([0.5 3.5 -3 145])
axis([0.5 3.5 -3 125])

% legend('TD','ASD','OND','ONDE')
box off

%% Plot significance stars
heights = [112.5,120];
% heights = [131,137,145];
drp = 2.5;

stars_line(2,heights(1),1,2,drp) % 3 stars,h,1,2,1
stars_line(3,heights(2),1,3,drp) % 3 stars,h,1,2,1
% stars_line(1,heights(3),2,3,drp) % 2 stars,h,1,3,2

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

function [adj] = adj_snap2zones(adj,num)
    for it = 1:num
        adj(it,4:7)=adj(it,4:7)+adj(it,13:16);    % reconnect to 4-7
        adj(4:7,it)=adj(4:7,it)+adj(13:16,it);    % reconnect to 13-16
        adj(it,13:16)=zeros(1,4);                 % remove non-food connections
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