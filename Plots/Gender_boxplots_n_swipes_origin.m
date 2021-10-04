% load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\OBJ_end_proport_110721.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')
load('n_swipes_origin.mat')
n_swipes = n_swipes';

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
    add_to = [n_swipes(sets{order(i)})];%[ranked(sets{order(i)})];
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

axis([0 12 0 130])
% hold on
% scatter(1.32,49,'pk','filled')
% scatter(1.5,49,'pk','filled')
% scatter(1.68,49,'pk','filled')
% scatter(10.5,49,'pk','filled')
% hold on
% height=49;
% plot([10,11],[height-1.5,height-1.5],'k')
% plot([10,10],[height-1.5,height-2.5],'k')
% plot([11,11],[height-1.5,height-2.5],'k')
% height=49;
% plot([1,2],[height-1.5,height-1.5],'k')
% plot([1,1],[height-1.5,height-2.5],'k')
% plot([2,2],[height-1.5,height-2.5],'k')
box off
ylabel('No. of swipes originating in zone 2')
xlabel('Gender')
f.Position = [403,340,574,313];
axis([0 12 0 125])
save_p = zeros(1,4);

numopt = [1,5;2,6;3,7;4,8];
for j = 1 : 4
    num=numopt(j,:);
    len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
%     pval = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds,'off');
    pval = kruskalwallis([n_swipes(sets{num(1)});n_swipes(sets{num(2)})],len_rankeds,'off');
    save_p(1,j) = pval;
end