% % % check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
option = 2; % 1 == proportional, 2 == proportion + swipe volume
if option == 1
    load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate_proport.mat')
%     load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')
elseif option == 2
    load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate.mat')
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

%% Main Plots
% load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\OBJ_end.mat')
% load('temp_save_OBJ_slow.mat')%load('temp_save_OBJ_fast.mat')
load('subject_details.mat')
% [sets] = set_allocate(subject_details_776,nam_save,saved);
[sets] = set_allocate(subject_details_776,nam_save,saved);
%%%  Remove ADHD
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%%
%%% Plot all
[prcnt,x] = plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets);
% title('OBJ')
% title('2 year 6 months – 3 year 8 months')%('3 year 9 months – 4 year 10 months')%('4 year 11 months – 6 year 0 months')

%% Difference All
figure;
plot(x,zeros(1,length(x)),'-o','LineWidth',1.5,'MarkerSize',5);
hold on
ref=prcnt{1};
for i = 2 : length(prcnt)
    diff{i}=prcnt{i}-ref;
    plot(x,diff{i},'-o','LineWidth',1.5,'MarkerSize',5)
    hold on
end
legend('TD','ASD','OND*','ONDE','Location','SouthWest')
% ax = gca;
% ax.XAxisLocation = 'origin';
box off
% 
xlabel('Perturbation Threshold (proportion)')
ylabel('Difference with respect to TD %')
% % title('4 year 11 months – 6 year 0 months')

% % %% Difference Gender
% figure;
% plot(x,zeros(1,length(x)));
% hold on
% ref=prcnt{4};i=1;
% % for i = 2 : length(prcnt)
%     diff{i}=prcnt{8}-ref;
%     plot(x,diff{i})
% %     hold on
% % end
% legend('ONDE_F','ONDE_M')
% % ax = gca;
% % ax.XAxisLocation = 'origin';
% box off
% 
% xlabel('Perturbation magnitude')
% ylabel('Difference with respect to ONDE_F %')
% % title('2 year 6 months – 3 year 8 months')
% 

%% Correlation
load('diags_incoming.mat')

load('subject_details.mat')
[months] = list_AGE(subject_details_776,nam_save,saved);
[sets] = set_allocate(subject_details_776,nam_save,saved);

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

if min(ranked)>=0   % change ranked to match pert ** NEED TO EXCLUDE 0.5s
    ranked = (ranked-1)*pert_chng+pert_init;
end

for num = 1 : 4
    figure;
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=months(sets{num});
    y=ranked(sets{num});
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
    
    if option == 1
        axis([27 73 -.3 .3])
    elseif option == 2
       axis([27 73 -30 30])
    end

    xlabel('Age (months)')
    if option == 1
        ylabel('Pertubation threshold (proportion)')
    elseif option == 2
        ylabel('Pertubation threshold (proportion + swipe volume)')
    end
    legend('subject','2nd order fit','Location','SouthEast')
    
    if num == 1
        title('TD')
    elseif num == 2
        title('ASD')
    elseif num == 3
        title('OND')
    elseif num == 4
        title('ONDE')
    end
end

%% OND Ranked trendline and correlation
folder1='H:\My Documents\MATLAB\Colormaps\Colormaps';
addpath(folder1)
clrs=fake_parula(11);

load('subject_details.mat')
load('OND_details.mat')
[months] = list_AGE(subject_details_776,nam_save,saved);
[sets] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);

if size(months,1)<size(months,2)
    months=months';
end

exclude = find(ranked==0.5);
exclude = [exclude;find(months>75)];
for i = 1 : length(sets)
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
end


figure;
for num = 1 : 4
    
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=months(sets{num});
    y=ranked(sets{num});
    if ~isempty(x)
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
        clr = clrs((num-1)*2+num,:);
        mrkr = 's';
    elseif num == 2
        clr = clrs((num-1)*2+num,:);
        mrkr = 'o';
    elseif num == 3
        clr = clrs((num-1)*2+num,:);
        mrkr = 'p';
    elseif num == 4
        clr = clrs((num-1)*2+num,:);
        mrkr = 'd';
    end
    tempclr=[.5 .5 .5]; 
%     scatplot=plot(x(Ind),y(Ind),mrkr,'color',clr,'linewidth',1);
    scatplot=scatter(x(Ind),y(Ind),55,mrkr,'MarkerFaceColor',clr,'MarkerEdgeColor',clr); 
    % Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
    scatplot.MarkerFaceAlpha = .3;
%     scatplot.MarkerEdgeAlpha = 0;

    hold on
    plot(x(Ind),y_fit(Ind),linetype,'color',clr,'LineWidth',1.5)

%     [rho,pval] = corr(months(sets{num}),ranked(sets{num}),'Type','Kendall');
%     text(55,mean(ranked(sets{num})),['p_{Ken,\tau} = ',num2str(p)])
    if option == 1
        axis([27 73 -.3 .3])
    elseif option == 2
       axis([27 73 -30 30])
    end
    xlabel('Age (months)')
    if option == 1
        ylabel('Pertubation threshold (proportion)')
    elseif option == 2
        ylabel('Pertubation threshold (proportion + swipe volume)')
    end
    legend('subject','linear fit','Location','SouthEast')
  
end
end
% legend('ADHD','Down Syndrome','Language delay/disorder','Other')
legend('ADHD','','Down Syndrome','','Language','delay/disorder','Other','')

%% Boxplot age - Food as an Origin
load('subject_details.mat')
saved=ones(1,704); saved(list)=zeros(1,length(list));

iter=0;all_sets=[];grps=[];
zones=1:17;zones(4)=[];

% [sets] = set_allocate(subject_details_776,nam_save,saved);
%%%
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%%
for num = 1:4
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
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
for m=1:length(h)
    temp_m = length(h)-m+1;
    mm=rem(m,5);mm(mm==0)=1;
    patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
end

% text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
xticklabels({'TD','ASD','OND*','ONDE'});
% set(gca,'xticklabel',entries,'fontsize',10)

if option == 1
    ylabel('Perturbation threshold (proportion)')
elseif option == 2
    ylabel('Perturbation threshold (proportion + swipe volume)')
end

if option == 1
    axis([0.5 4.5 -.31 .39])
elseif option == 2
    axis([0.5 4.5 -31 39])
end

% legend('TD','ASD','OND','ONDE')
box off

%% Plot significance stars
if option == 1
    heights = [.29,.33,.39,.29];
    drp = .0175;
elseif option == 2
    heights = [28,32,37,28];
    drp = 1.75;
end
% height=28;%0.29; 
n_stars = 3; 
stars_line(n_stars,heights(1),1,2,drp) % 3 stars,h,1,2,1
% height = 32;%.33;
stars_line(n_stars,heights(2),1,3,drp) % 3 stars,h,1,2,1
% height = 37;%.39; 
stars_line(n_stars,heights(3),2,4,drp) % 2 stars,h,1,3,2
% height = 28;%.29;
stars_line(n_stars,heights(4),3,4,drp) % 2 stars,h,3,4,1

combos=nchoosek([1,2,3,4],2);
save_p = zeros(1,size(combos,1));
for j = 1 : size(combos,1)
    num = combos(j,:);
    len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
    pval = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds,'off');
    save_p(1,j) = pval;
end


% legend('TD','ASD','OND*','ONDE','Orientation','horizontal')

f.Position = [403,340,300,313];

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
