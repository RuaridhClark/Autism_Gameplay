load('C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj\swipes_all704.mat')
swipe=swipe_save{i};
name=nam_save{i};

tab=sortrows(readtable('GLA-SL-084.Sharing.TouchData.typed.csv'),4);

Img=imread('C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Figures\Figures_0821\tablet_segments_swipes_update3.png'); %tablet_bw.png');

figure;
imshow(Img)
truesize([1023.2 768]);

% load("list.mat")
% for i = 1 : length(list)
%     m=list(i);
%     k = swipe{m};
%     hold on
%     plot(tab.X(k),768-tab.Y(k),'color','b','linewidth',1.5)
%     pause(.2)
% end
for m = 1 : length(swipe)
    k = swipe{m};
    hold on
    plot(tab.X(k),768-tab.Y(k),'color','b','linewidth',1.5)
    pause(.2)
end