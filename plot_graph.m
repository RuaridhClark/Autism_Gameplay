% load('H:\My Documents\MATLAB\Autism_MAIN\adjs_110721\adj_obj_end\subject_GLA-EH-002.mat', 'adj')
adj=adj-diag(diag(adj));
Adj=[]

% allow=[2,4,5,6,7];
% for it = 1:16
%     if ~ismember(it,allow)
%         adj(it,4:7)=zeros(1,4);
%     end
% end

yl = 768;
obj(1,1:4)=[19, 227, yl-752, yl-525];
obj(2,1:4)=[234, 774, yl-678, yl-525];
obj(3,1:4)=[781, 1004, yl-752, yl-525];
obj(4,1:4)=[101, 101+131, yl-400-96, yl-400];   
obj(5,1:4)=[325,325+131, yl-400-96, yl-400];    
obj(6,1:4)=[551, 551+131, yl-400-96, yl-400];   
obj(7,1:4)=[776, 776+131, yl-400-96, yl-400];   
obj(8,1:4)=[19, 252, yl-190, yl-74];
obj(9,1:4)=[257, 420, yl-190, yl-27];
obj(10,1:4)=[426, 589, yl-190, yl-27];
obj(11,1:4)=[595, 758, yl-190, yl-27];
obj(12,1:4)=[770, 1005, yl-260, yl-52];
obj(13,1:4)=[19, 309, yl-523, yl-196];
obj(14,1:4)=[311, 531, yl-523, yl-196];
obj(15,1:4)=[533,753, yl-523, yl-196];
obj(16,1:4)=[755, 1004, yl-523, yl-264];

x=[];y=[];
for i = 1:16
    x(i) = obj(i,1)+(obj(i,2)-obj(i,1))/2;
    y(i) = obj(i,3)+(obj(i,4)-obj(i,3))/2;
end


G=digraph(adj);
figure
G.Edges.LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
h=plot(G,'XData',x,'YData',y,'LineWidth',G.Edges.LWidths,'EdgeColor','b','NodeColor','k');
labelnode(h,[8 9 10 11 12 13 14 15 16],{'12' '13' '14' '15' '16' '8' '9' '10' '11'})
axis off
title('TD GLA-EH-049')