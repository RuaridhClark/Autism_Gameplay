clear all
folder1 = 'H:\My Documents\MATLAB\Autism\EEG_eigalign_validate';
folder2 = 'H:\My Documents\MATLAB\Autism\EEG_eigalign_validate\functions';
folder3 = 'H:\My Documents\MATLAB\Autism\adj';
addpath(folder1,folder2,folder3)
file_loc = 'H:\My Documents\MATLAB\Autism\adj\';

A_num = 0;
T_num = 0;
ASD_adj=[];
TD_adj=[];
%% stack the adjs
for i = 1:82
    if i<10
        ntext = '0';
    else
        ntext = '';
    end
    ASD_id = ['ASD_subject',ntext,num2str(i),'.mat'];
    TD_id = ['TD_subject',ntext,num2str(i),'.mat'];
    
    if isfile([file_loc,ASD_id])
        A_num = A_num + 1;
        load(ASD_id)
%         adj=adj./sum(adj,2);
%         adj(isnan(adj))=0;
        ASD_adj(:,:,A_num)=adj;
        titlename = ['ID ',ntext,num2str(i),' - ASD'];
        savename = ['ASD_subject',ntext,num2str(i)];
    elseif isfile([file_loc,TD_id])
        T_num = T_num + 1;
        load(TD_id)
%         adj=adj./sum(adj,2);
%         adj(isnan(adj))=0;
        TD_adj(:,:,T_num)=adj;
        titlename = ['ID ',ntext,num2str(i),' - TD'];
        savename = ['TD_subject',ntext,num2str(i)];
    end
    
    
    row_l = 9;
    col_l=6;
    hm = zeros(col_l,row_l);
    for j = 1 : length(adj)
        x = rem(j,row_l);
        x(x==0)=row_l;
        y = ceil(j/row_l);
        hm(col_l+1-y,x) = sum(adj(j,:)); 
    end
    figure;
    heatmap(1:row_l,sort(1:col_l,'desc'),hm)
    title(titlename)
    saveas(gcf,['heatmap_',savename,'.png'])
end

