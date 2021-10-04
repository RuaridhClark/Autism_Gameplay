clear all
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\paths.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Track_Finger\Data\swipe_speeds.mat')

mn_plt_speed = zeros(704,1);
rng = 8;
for m = 1 : 704
    trgt_speed = zeros(17,17);
    trgt_swipes = zeros(17,17);
    for i = 1 : length(speed{m})
        path = all_paths{m};
        if speed{m}(i,1)>0 && min(path(i,:))>0 && ~isinf(speed{m}(i,1))
            trgt_speed(path(i,1),path(i,2)) = trgt_speed(path(i,1),path(i,2)) + speed{m}(i,1);
            trgt_swipes(path(i,1),path(i,2))= trgt_swipes(path(i,1),path(i,2)) + 1;
        end
    end
    
    mn_speed = sum(trgt_speed(2,rng))/sum(trgt_swipes(2,rng));
    mn_speed(isnan(mn_speed))=0;
    mn_plt_speed(m) = mn_speed;
end

%% sets
folder4 = 'H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
addpath(folder4)
load('H:\My Documents\MATLAB\Autism_MAIN\subject_details.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\swipes_all704.mat','nam_save')
[months] = find_ages(subject_details_776,nam_save);

load('H:\My Documents\MATLAB\Autism_MAIN\Ranking_Main\Data\save_OBJ_end.mat','ranked','saved')

[sets] = set_allocate(subject_details_776,nam_save,saved);
% exclude = find(ranked==0.5);
% exclude = [exclude;find(months>75)];
% for i = 1 : length(sets)
%     rmv=find(ismember(sets{i},exclude')==1);
%     sets{i}(rmv)=[];
% end



%% plots
for num = 1 : 4
    figure;
%     f=fit(months(sets{num}),ranked(sets{num}),'poly1');
%     plot(f,months(sets{num}),ranked(sets{num}),'x')
    
    x=months(sets{num})';
    y=mn_plt_speed(sets{num});
    x(find(y==0))=[];
    y(y==0)=[];
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
    text(0.05,.95,['p_{Ken,\tau} = ',num2str(p)],'Units','normalized')
%     axis([27 73 -20 25])
    xlabel('Age (months)')
    ylabel('Mean Speed to Plate from Food')
    legend('subject','linear fit','Location','SouthEast')
    
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
    
%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%
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