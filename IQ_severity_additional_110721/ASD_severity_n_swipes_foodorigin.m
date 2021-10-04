clear all
load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Correlations_110721\Data\save_OBJ_end.mat')

fileloc='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\IQ_severity';
tab_sev = readtable([fileloc,'\eCRF.csv']);

folder3 = 'H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end';
folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
folder5 = 'H:\My Documents\MATLAB\Autism_MAIN\Plots';
folder6 = 'H:\My Documents\MATLAB\Autism_MAIN\Create_adj_110721';
folder7 = 'H:\My Documents\MATLAB\Autism_MAIN';
addpath(folder3,folder4,folder5,folder6,folder7)
file_loc = 'H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end\'; % should match zone type

load('swipes_all704.mat','nam_save')

%% stack the adjs
num =12;    % number of ipad objects (nodes)

asdname_save = cell(704,1);
all_severity = zeros(704,1);
sets = cell(3,1);
check = 0;

map = [];
    m = 0;
    n_swipes = [];
    for i = 1:height(tab_sev)
        file_id = ['subject_',tab_sev.id_study_id{i},'.mat'];
        
        severity = tab_sev.clinical_diagnosis__asd_severity_level{i};
        if isfile([file_loc,file_id]) && ~isempty(severity)
            map(i)=find(strcmp(nam_save,tab_sev.id_study_id{i}));
            load(file_id)
            adj = adj(1:12,1:12);
            n_swipes(i) = sum(adj(:,2));
            [sev_num] = severity_score(severity); 
            sets{sev_num}=[sets{sev_num},map(i)];
            all_severity(i) = sev_num;
            asdname_save{i} = tab_sev.id_study_id{i};
            titlename = ['ID ',tab_sev.id_study_id{i}];
            savename = ['subject_',tab_sev.id_study_id{i}];
        end
    end
    
%% Correlation
load('subject_details.mat')
[tmp_months] = list_AGE(subject_details_776,asdname_save,saved);
months=zeros(1,694);
months(1:length(tmp_months))=tmp_months;
n_swipes=[n_swipes,0];

rmv = find(map==0);
months(rmv)=[];
map(rmv)=[];
n_swipes(rmv)=[];

saved=saved(map);

rmv = find(sum(saved,1)==0);
map(rmv)=[];
n_swipes(rmv)=[];
ranked=ranked(map);
ranked(ranked==0.5)=[];

if size(months,1)<size(months,2)
    months=months';
end

exclude = find(months>75);
for i = 1 : length(sets)
    for jj = 1 : length(sets{i})
        sets{i}(jj)=find(map==sets{i}(jj));
    end
    
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
end

%% plot_cmprsn4 changes ranked internally to match pert
[prcnt,x] = plot_cmprsn4(n_swipes,pert_init,pert_chng,nam_save,sets);
legend('1','2','3','Location','SouthWest')

%% Boxplot age - Food as an Origin
load('subject_details.mat')
saved=ones(1,704); saved(list)=zeros(1,length(list));

iter=0;all_sets=[];grps=[];
zones=1:12;

for num = 1:3
%         h = kstest(diagsA(zones(2),sets{num}))
    iter=iter+1;
    vals = n_swipes(sets{num});
    all_sets = [all_sets,vals];
    grps = [grps,ones(1,length(vals))*iter];
end

f=figure;
boxplot(all_sets,grps,'Notch','on')
h = findobj(gca,'Tag','Box');
colors = [.93,.69,.13;.85,.33,.1;.64,.08,.18;];%colors = [0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.5560];
for k=1:length(h)
    temp_m = length(h)-k+1;
    mm=rem(k,5);mm(mm==0)=1;
    patch(get(h(temp_m),'XData'),get(h(temp_m),'YData'),colors(mm,:),'FaceAlpha',.5);
end

% text_x = {'2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months'};
xticklabels({'1','2','3'});
% set(gca,'xticklabel',entries,'fontsize',10)

xlabel('ASD severity')
ylabel('No. of swipes originating in zones 2')
box off

% num_sets=[];
% for i = 1 : length(sets)
%     num_sets = ['n = ',num2str(length(sets{i}))];
%     text(0.1+(i-1)/3,1,num_sets,'Units','normalized')
% end

combos=nchoosek([1,2,3],2);
save_p = zeros(1,size(combos,1));
for j = 1 : size(combos,1)
    num = combos(j,:);
    len_rankeds = [ones(length(sets{num(1)}),1);2*ones(length(sets{num(2)}),1)];
    pval = kruskalwallis([n_swipes(sets{num(1)})';n_swipes(sets{num(2)})'],len_rankeds,'off');
    save_p(1,j) = pval;
end

f.Position = [403,340,275,313];
axis([0.5 3.5 -3 140])

%% Plot significance stars
height=133; drp = 2.5;
% stars_line(2,height,1,2,drp) % 3 stars,h,1,2,1
height = 135;
stars_line(3,height,1,3,drp) % 3 stars,h,1,2,1
height = 128;
stars_line(1,height,2,3,drp) % 3 stars,h,1,2,1

%%%%%%%%%%%%%% Function %%%%%%%%%%%%%%%%
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