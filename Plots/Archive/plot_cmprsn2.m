function [] = plot_cmprsn2(ranked,pert_init,pert_chng)

% [1,2:328,329:351,352:495,496:500,501:504,505:554,555:580,581:582,583,584:601,602,603,604:693,694:704]

ASD_n = [329:351];
TD_n = [352:495];

prcnt_ASD = [];
prcnt_TD = [];

for i = 1 : max(ranked)
    prcnt_ASD(i)=100*length(find(ranked(ASD_n)>i))/length(ASD_n);
    prcnt_TD(i)=100*length(find(ranked(TD_n)>i))/length(TD_n);
end

figure;
x=linspace(pert_init,pert_init+length(prcnt_ASD)*pert_chng-1,length(prcnt_ASD));
% x=linspace(-40,22,63);
plot(x,prcnt_ASD,'-s')
hold on
plot(x,prcnt_TD,'-s')

legend('ASD','TD')
box off
xlabel('Perturbation magnitude')
ylabel('Percentage of subjects')
grid on

% saveas(gcf,[savename,'.png'])
% saveas(gcf,'prcnt_4by6.png')
end