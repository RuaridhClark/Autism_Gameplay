load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\OBJ_end.mat')
load('OND_details.mat')

ranked = ranked-20;

iter=0;
all_sets=[];
grps=[];
order = [1,2,3,4];
combos=nchoosek([1,2,3,4],2);
save_p = zeros(3,size(combos,1));

for age = 1:3
    [sets] = set_allocate_AGE(subject_details_776,nam_save,saved,age);
    %%  Remove ADHD
    [sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);
    ids = 3;
    curr_set = sets{ids};
    keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
    sets{ids}=curr_set(ismember(curr_set,keep));
    %%
    for j = 1 : size(combos,1)
        num = combos(j,:);
        len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
        pval = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds,'off');
        save_p(age,j) = pval;
    end

    for i = 1 : length(order)    
        iter = iter +1;
        add_to = [ranked(sets{order(i)})];
        all_sets = [all_sets;add_to];
        grps = [grps;ones(length(add_to),1)*iter];
        if rem(i,4)==0
            iter=iter+1;
            grps = [grps;iter];
            all_sets = [all_sets;0];
        end
    end
    saved_p{age}=save_p;
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
hAx=gca;
xt=hAx.XTick;                            % retrieve ticks
xt(5:5:end)=[];
hAx.XTick=xt;                   % clear some
entries = cell(1,12);
% text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
text_x = {'            30–44 months','            45–58 months','            59–72 months'};
for ind = 1 : 3
    %         entries = num2str(ind)
    entries(2+(ind-1)*4)=text_x(ind);
end
% xticklabels(entries)
set(gca,'xticklabel',entries,'fontsize',10)
ylabel(['Perturbation threshold'])

% text(1.3,-5,'TD');text(4.225,-5,'ASD');text(7.225,-5,'OND');text(10.175,-5,'ONDE')

axis([0.5 14.5 -20 35])
%%
height=136.5/4.3235; n_stars = 3;
stars_line(n_stars,height,1,2,2) % 3 stars,h,1,2,2
height = 140/4.3235;
stars_line(n_stars,height,1,2,3) % 3 stars,h,1,2,3
height = 120/4.3235; n_stars = 2;
stars_line(n_stars,height,1,3,1) % 2 stars,h,1,3,1
height = 147/4.3235;
stars_line(n_stars,height,1,3,2) % 2 stars,h,1,3,2
height = 130/4.3235;
stars_line(n_stars,height,2,4,2) % 2 stars,h,2,4,2
height = 133/4.3235;
stars_line(n_stars,height,2,4,3) % 2 stars,h,2,4,3
height = 112/4.3235;
stars_line(n_stars,height,3,4,2) % 2 stars,h,3,4,2
box off
ylabel('Perturbation threshold')
f.Position = [403,340,674,313];

% num = [1,5];
% len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
% p = kruskalwallis([ranked(sets{num(1)});ranked(sets{num(2)})],len_rankeds);

function [] = stars_line(n_stars,height,strt,nd,age)
    drp = 5/4.3235; drp2 = 8/4.3235;
    hold on
    if n_stars == 3
        scatter((age-1)*5+(nd-strt)/2+(strt-0.18),height,'pk','filled')
        scatter((age-1)*5+(nd-strt)/2+(strt),height,'pk','filled')
        scatter((age-1)*5+(nd-strt)/2+(strt+0.18),height,'pk','filled')
    elseif n_stars == 2
        scatter((age-1)*5+(nd-strt)/2+(strt-0.09),height,'pk','filled')
        scatter((age-1)*5+(nd-strt)/2+(strt+0.09),height,'pk','filled')
    end
    hold on
    plot([(age-1)*5+strt,(age-1)*5+nd],[height-drp,height-drp],'k')
    plot([(age-1)*5+strt,(age-1)*5+strt],[height-drp,height-drp2],'k')
    plot([(age-1)*5+nd,(age-1)*5+nd],[height-drp,height-drp2],'k')
end