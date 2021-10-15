% check NNR_adj_conns_OBJ2 and pert changes for velocity case
clear all
load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\OBJ_end_accurate.mat')
% load('H:\My Documents\GitHub\Autism_Gameplay\Ranking_Correlations_110721\Data\save_OBJ_end.mat')

% folder1 = 'H:\My Documents\MATLAB\Autism_MAIN\EEG_eigalign_validate';
folder2 = 'H:\My Documents\GitHub\Autism_Gameplay\Set_allocate';
folder3 = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate';
folder5 = 'H:\My Documents\GitHub\Autism_Gameplay\Plots';
folder6 = 'H:\My Documents\GitHub\Autism_Gameplay\Create_adj_110721';
folder7 = 'H:\My Documents\GitHub\Autism_Gameplay';
addpath(folder2,folder3,folder5,folder6,folder7)
file_loc = 'H:\My Documents\GitHub\Autism_Gameplay\adjs_110721\adj_obj_end_accurate\'; % should match zone type


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

    
    n_swipes(jj) = sum(adj(2,[13,14,15,16]));%[4,5,6,7]));
    adj=adj-diag(diag(adj));
%     diagsA(1:12,jj)=adj(2,:);%sum(adj,1);
end

% val = sum(diagsA(4:7,:));
% % diagsA(4:7,:)=val.*ones(4,704);

%% Setup 
load('subject_details.mat')
saved=ones(1,704);
saved(list)=zeros(1,length(list));
[sets] = set_allocate(subject_details_776,nam_save,saved);
%%
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%
[months] = list_AGE(subject_details_776,nam_save,saved);

%% Boxplot age - Food as an Origin
load('subject_details.mat')
saved=ones(1,704); saved(list)=zeros(1,length(list));

iter=0;all_sets=[];grps=[];
zones=1:12;

[sets] = set_allocate(subject_details_776,nam_save,saved);
%%%
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%%
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

for num = 1:4
%         h = kstest(diagsA(zones(2),sets{num}))
    iter=iter+1;
%         vals = diagsA(zones(zns),sets{num});
    vals = n_swipes(1,sets{num});
    all_sets = [all_sets,vals];%val(sets{num})];
    grps = [grps,ones(1,length(vals))*iter];%val(sets{num})))*iter];
end

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

ylabel('No. of swipes ending in zones 4--7')
% axis([0.5 4.5 0 165])
% % legend('TD','ASD','OND','ONDE')
% box off
% 
% %%%
n_stars = 3;
height=140;%145; 
stars_line(n_stars,height,1,2,1) % 3 stars,h,1,2,1
height = 128;%133;
stars_line(n_stars,height,1,3,1) % 3 stars,h,1,2,1
height = 147;%152; 
stars_line(n_stars,height,2,4,1) % 2 stars,h,1,3,2
n_stars = 2;
height = 122;%127; 
stars_line(n_stars,height,3,4,1) % 2 stars,h,3,4,1

combos=nchoosek([1,2,3,4],2);
save_p = zeros(3,size(combos,1));
% [sets] = set_allocate(subject_details_776,nam_save,saved);
%%%
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
curr_set = sets{3};
keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
sets{3}=curr_set(ismember(curr_set,keep));
%%%
for j = 1 : size(combos,1)
    num = combos(j,:);
    len_rankeds = [ones(1,length(sets{num(1)})),2*ones(1,length(sets{num(2)}))];
    pval = kruskalwallis([n_swipes(sets{num(1)}),n_swipes(sets{num(2)})],len_rankeds,'off');
    save_p(1,j) = pval;
end


% legend('TD','ASD','OND*','ONDE','Orientation','horizontal')
axis([0.5 4.5 -3 150])%axis([0.5 4.5 -3 155])
f.Position = [403,340,300,313];
box off

%% Correlation    
% if size(months,1)<size(months,2)
%     months=months';
% end
% 
% exclude = find(ranked==0.5);
% exclude = [exclude;find(months>75)];
% for i = 1 : length(sets)
%     rmv=find(ismember(sets{i},exclude')==1);
%     sets{i}(rmv)=[];
% end

for num = 1 : 4
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
%     text(55,mean(diagsA(5,sets{num})),['p_{Ken,\tau} = ',num2str(p)])
    text(0.05,.95,['p_{Ken,\tau} = ',num2str(p)],'Units','normalized')
    axis([27 73 0 110])%axis([27 73 0 125])
    xlabel('Age (months)')
    ylabel('No. of swipes - food to plate')
    legend('subject','2nd order fit','Location','SouthEast')
    
    if num == 1
        title('TD')
    elseif num == 2
        title('ASD')
    elseif num == 3
        title('OND*')
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
[sets,other] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);

if size(months,1)<size(months,2)
    months=months';
end

exclude = find(ranked==0.5);
exclude = [exclude;find(months>75)];
for i = 1 : length(sets)
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
    if i == 4
        other(rmv)=[];
    end
end


figure;
for num = 1 : 4
    
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=months(sets{num});
    y=n_swipes(sets{num})';
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
%     axis([27 73 -20 25])
    xlabel('Age (months)')
    ylabel('No. of swipes - food to plate')
    legend('subject','linear fit','Location','SouthEast')
    
%     if num == 1
%         title('TD')
%     elseif num == 2
%         title('ASD')
%     elseif num == 3
%         title('OND')
%     elseif num == 4
%         title('ONDE')
%     end
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
    vals = n_swipes(sets{num});
    all_sets = [all_sets,vals];
    grps = [grps,ones(1,length(vals))*iter];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%
function [] = stars_line(n_stars,height,strt,nd,age)
    drp = 5; drp2 = 8;
    hold on
    if n_stars == 3
        scatter((age-1)*5+(nd-strt)/2+(strt-0.18),height,'pk','filled')
        scatter((age-1)*5+(nd-strt)/2+(strt),height,'pk','filled')
        scatter((age-1)*5+(nd-strt)/2+(strt+0.18),height,'pk','filled')
    elseif n_stars == 2
        scatter((age-1)*5+(nd-strt)/2+(strt-0.09),height,'pk','filled')
        scatter((age-1)*5+(nd-strt)/2+(strt+0.09),height,'pk','filled')
    elseif n_stars == 1
        scatter((age-1)*5+(nd-strt)/2+(strt),height,'pk','filled')
    end
    hold on
    plot([(age-1)*5+strt,(age-1)*5+nd],[height-drp,height-drp],'k')
    plot([(age-1)*5+strt,(age-1)*5+strt],[height-drp,height-drp2],'k')
    plot([(age-1)*5+nd,(age-1)*5+nd],[height-drp,height-drp2],'k')
end