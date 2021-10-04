load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')
% load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')
% load('n_swipes.mat')
% n_swipes = n_swipes';

if min(ranked)>=0   % change ranked to match pert ** NEED TO EXCLUDE 0.5s
    ranked = (ranked-1)*pert_chng+pert_init;
end

[sets] = set_allocate_GENDER_TYPE(subject_details_776,nam_save,saved);
%
load('OND_details.mat')
[sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
ids = [3,7]; % Female,Male
for i = 1:2
    curr_set = sets{ids(i)};
    keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
    sets{ids(i)}=curr_set(ismember(curr_set,keep));
end
%

iter=0;
all_sets=[];
grps=[];
order = [1,5,2,6,3,7,4,8];

for i = 1 : length(order)    
    iter = iter +1;
    add_to = [ranked(sets{order(i)})];%[n_swipes(sets{order(i)})];%
    all_sets = [all_sets;add_to];
    grps = [grps;ones(length(add_to),1)*iter];
    if rem(i,2)==0
        iter=iter+1;
        grps = [grps;iter];
        all_sets = [all_sets;-100];
    end
end
    
f=figure;
boxplot(all_sets,grps,'Notch','on')
h = findobj(gca,'Tag','Box');
colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
for m=1:length(h)
    temp_m = length(h)-m+1;
    mm=ceil(m/3);
    patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
end

hAx=gca; 
xt=hAx.XTick;                   % retrieve ticks
xt(3:3:end)=[];
hAx.XTick=xt;                   % clear some

entries ={'F','M','F','M','F','M','F','M'};
xticklabels(entries)

% text(1.3,-5,'TD');text(4.225,-5,'ASD');text(7.225,-5,'OND');text(10.175,-5,'ONDE')

axis([0 12 -.3 .3])%axis([0 12 0 50])

%% Plot significance stars
height=0.29; n_stars = 3; drp = .0175;
stars_line(n_stars,height,1,2,drp) % 3 stars,h,1,2,1
height = .29; n_stars = 1;
stars_line(n_stars,height,10,11,drp) % 2 stars,h,3,4,1

% legend('TD','ASD','OND*','ONDE','Orientation','horizontal')

xlabel('Gender')
ylabel('Perturbation threshold (proportion)')
f.Position = [403,340,574,313];
box off

save_p = zeros(1,4);

numopt = [1,5;2,6;3,7;4,8];
for j = 1 : 4
    num=numopt(j,:);
    len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
    pval = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds,'off');
%     pval = kruskalwallis([n_swipes(sets{num(1)});n_swipes(sets{num(2)})],len_rankeds,'off');
    save_p(1,j) = pval;
end

%%%%%%%%%% function %%%%%%%%%%%%%
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