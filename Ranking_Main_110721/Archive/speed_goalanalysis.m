clear all
folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
addpath(folder4)

load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\save_OBJ_end_KFslow.mat','ranked')
ranked_all_slow = ranked;

load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\save_OBJ_end_slow_mvmean.mat','ranked')
ranked_age_slow = ranked;

load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\save_OBJ_end_KFfast.mat','ranked')
ranked_all_fast = ranked;

load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\save_OBJ_end_fast_mvmean.mat','ranked','saved')
ranked_age_fast = ranked;

load('H:\My Documents\MATLAB\Autism_MAIN\subject_details.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\swipes_all704.mat','nam_save')

[months] = find_ages(subject_details_776,nam_save);
if size(months,1)<size(months,2)
    months=months';
end

[sets] = set_allocate(subject_details_776,nam_save,saved);

diff_all = ranked_all_fast-ranked_all_slow;
diff_age = ranked_age_fast-ranked_age_slow;

exclude = find(ranked==0.5);
exclude = [exclude;find(months>75)];
for i = 1 : length(sets)
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
end

for m = 1 : 2
    if m == 1
        diff = diff_all;
        addtext='all';
    else
        diff = diff_age;
        addtext='age';
    end
    for num = 1 : 4
        figure;
    %     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
    %     plot(f,months(sets{num}),ranked(sets{num}),'x')

        x=months(sets{num});
        y=diff(sets{num});
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
%         text(10,mean(diff(sets{num})),['p_{Ken,\tau} = ',num2str(p)])
        text(0.025,.95,['p_{Ken,\tau} = ',num2str(p)],'Units','normalized')
        axis([30 75 -15 20])
        xlabel('Age (months)')
        ylabel('Pertubation threshold difference (fast - slow)')
        legend('subject','linear fit','Location','SouthEast')

        if num == 1
            title(['TD ',addtext])
        elseif num == 2
            title(['ASD ',addtext])
        elseif num == 3
            title(['OND ',addtext])
        elseif num == 4
            title(['ONDE ',addtext])
        end
    end
end

%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%
function [months] = find_ages(subject_details,nam_save)

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
%         if max(saved(:,I))>0  %% removing subjects with no sharing game data
            yr = subject_details{i,2} + subject_details{i,3}/12;
            mnths = subject_details{i,2}*12 + subject_details{i,3};
            months(I)=mnths;
%         end
    end

end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end