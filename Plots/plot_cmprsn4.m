function [prcnt,x] = plot_cmprsn4(ranked,pert_init,pert_chng,nam_save,sets)

% sets = [{1},{2:328},{329:351},{352:495},{496:500},{501:504},{505:554},{555:580},{581:582},{583},{584:601},{602},{603},{604:693},{694:704}];

% legend_names = [];
prcnt = {};

min_len = 1; % prevent small sets
for i = 1 : length(sets)
    if length(sets{i})>min_len
        percentage = [];
        for ii = 1 : max(ranked)
            minus_len = length(find(ranked(sets{i})==-1));
            percentage(ii)=100*length(find(ranked(sets{i})>ii))/(length(sets{i})-minus_len);
        end
        prcnt{i}=percentage;
%         temp_name = nam_save(sets{i}(1));
%         legend_names = [legend_names;temp_name{1}(1:end-4)];
    end
end

figure;
% x=linspace(pert_init,pert_init+max(ranked)*pert_chng-1,max(ranked));
x=linspace(pert_init,pert_init+(max(ranked)-1)*pert_chng,max(ranked));
% x=linspace(-40,22,63);
for i = 1 : length(sets)
    if length(sets{i})> min_len
        plot(x,prcnt{i},'-o','LineWidth',1.5,'MarkerSize',5)
        hold on
    end
end

legend('TD','ASD','OND*','ONDE','Location','SouthWest')
box off
xlabel('Perturbation')
ylabel('Percentage of subjects')
grid on

% saveas(gcf,[savename,'.png'])
% saveas(gcf,'prcnt_4by6.png')
end