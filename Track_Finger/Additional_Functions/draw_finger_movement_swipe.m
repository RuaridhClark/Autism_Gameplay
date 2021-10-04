% close all
% clear all
file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\Files\segmented_by_SCL\';
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\swipes_all704.mat')
% load('swipe_save.mat')

lista = [1,4,7,8,9,12,13,14,16,20,22,24,26,28,32,36];
allASD = 1:37;
lista_inv = allASD; lista_inv(lista)=[];
listb_inv = [2,13,23,29,43,45];
allTD = 1:45;
listb = allTD; listb(listb_inv)=[];
m = 0; n = m;
for i = 1:82
%     close all
    swipe=swipe_save{i};
%     figure;
    skip=0;
    if i<10
        ntext = '0';
    else
        ntext = '';
    end
    ASD_id = ['subject',ntext,num2str(i),'_ASD_sharing_touch.csv'];
    TD_id = ['subject',ntext,num2str(i),'_TD_sharing_touch.csv'];
    
    if isfile([file_loc,ASD_id])
        n=n+1;
%         if ismember(n,lista)
            colr = [0.8500, 0.3250, 0.0980];
            titlename = ['ID ',ntext,num2str(i),' - ASD'];
            filename = [file_loc,ASD_id];
            savename = ['ASD_subject',ntext,num2str(i)];
%         else
%             skip = 1;
%         end
    elseif isfile([file_loc,TD_id])
        m=m+1;
%         if ismember(m,listb)
            colr = [0, 0.4470, 0.7410];
            titlename = ['ID ',ntext,num2str(i),' - TD'];
            filename = [file_loc,TD_id];
            savename = ['TD_subject',ntext,num2str(i)];
%         else
%             skip = 1;
%         end
    else
        skip = 1;
    end
    
    if skip == 0
        tab=sortrows(readtable(filename),7);

%         figure
%         %% plot
%         plot(tab.x,tab.y,colr)
%         axis([0,1025,0,800])

%         %% animation
%         h = animatedline;
%         axis([0,1050,0,770])
%         title(titlename)
%         prev_tP = 0;
%         alternate = 0;
%         for m = 1 : length(swipe)
%             for n = 1 : length(swipe{m})
%                 k = swipe{m}(n);
%                 if n>1
%                     addpoints(h,tab.x(k),tab.y(k));
%                     drawnow
% %                     pause(.1)
%                 else
%                     clearpoints(h)
%                 end
%             prev_tP = tab.touchPhase(k);
%             end
%         end
        
        figure;
        title(titlename)
        axis([0,1050,0,770])
        for m = 1 : length(swipe)
            k = swipe{m};
            hold on
            plot(tab.x(k),tab.y(k),'color',colr,'linewidth',1.5)
%             pause(.1)
        end
   
        axis off
%         export_fig transparent_swipes -transparent
%         saveas(gcf,[savename,'.png'])
    end
end