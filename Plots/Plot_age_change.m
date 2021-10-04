
for N =1:4
    figure;
%     clrs=viridis(4);
    iter=0;
    for age = 1 : 3
        [sets] = set_allocate_AGE(subject_details_776,nam_save,saved,age);
        [sets_OND] = set_allocate_TYPE_OND(subject_details_776,OND_details,nam_save,saved);

        legend_names = [];
        prcnt = {};

        min_len = 1; % prevent small sets
%         for i = 1 : length(sets)
            curr_set = sets{N};
            if N==3     % Remove ADHD
                keep = [sets_OND{2},sets_OND{3},sets_OND{4}];
                curr_set=curr_set(ismember(curr_set,keep));
            end
            if length(curr_set)>min_len
                percentage = [];
                for ii = 1 : max(ranked)
                    minus_len = length(find(ranked(curr_set)==-1));
                    percentage(ii)=100*length(find(ranked(curr_set)>ii))/(length(curr_set)-minus_len);
                end
                prcnt{i}=percentage;
%                 temp_name = nam_save(curr_set(1));
%                 legend_names = [legend_names;temp_name{1}(1:end-4)];
            end
%         end


        x=linspace(pert_init,pert_init+max(ranked)*pert_chng-1,max(ranked));
        % x=linspace(-40,22,63);
        colormap(autumn)
        
        if length(sets{i})> min_len
            iter=iter+1;
            [clr] = pick_colour(N,iter);
            plot(x,prcnt{i},'-s','color',clr,'linewidth',1.5)
            hold on
        end

    %     [prcnt,x] = plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets);
    %     hold on
    end

    box off
    xlabel('Perturbation')
    ylabel('Percentage of subjects')
    grid on
    if N == 1
        title('TD')
    elseif N == 2
        title('ASD')
    elseif N == 3
        title('OND')
    elseif N ==4
        title('ONDE')
    end

    legend('30–44 months','45–58 months','59–72 months','Location','southwest')
%     legend('2 years 6 months - 3 years 8 months','3 years 9 months - 4 years 10 months','4 years 11 months - 6 years 0 months', 'Location','southwest' )
end

function [clr] = pick_colour(num,age)

    if num == 1
        if age == 1
            clr = [0.5765,    0.6706,    1.0000];
        elseif age == 2
            clr = [0, 0.4470, 0.7410];
        elseif age == 3
            clr = [7,53,88]./255;
        end
    elseif num == 2
        if age==1
            clr = [1.0000    0.5804    0.3961];
        elseif age == 2
            clr = [0.8500, 0.3250, 0.0980];
        elseif age == 3
            clr = [104,42,1]./255;
        end
    elseif num == 3
        if age==1
            clr = [0.9290, 0.6940, 0.1250];
        elseif age == 2
            clr = [0.9290, 0.6940, 0.1250];
        elseif age == 3
            clr = [111,85,10]./255;
        end
    elseif num == 4
        if age==1
            clr = [0.7020    0.5373    0.8000];
        elseif age == 2
            clr = [0.5765, 0.2157, 0.6510];
        elseif age == 3
            clr = [62,17,65]./255;
        end
    end
end