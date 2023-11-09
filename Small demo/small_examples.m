adj = zeros(16,16);
% close all
figure()
legend_list = [{'Even plates'}, {'Zone 2 to 3'}, {'Zone 5 to 3'}, {'Inter-plate (Direct)'}, {'Inter-plate (Indirect)'}, {'Zone 5 to 2'}, {'Uneven plates'}];
all_legends=[];
list = [4:5];
for k = 1:length(list)
    i = list(k);
    adj = zeros(16,16);
    if i == 1
        %% example 1
        adj(2,4) = 9;
        adj(2,5) = 9;
        adj(2,6) = 9;
        adj(2,7) = 9;
        
        adj(3,1) = 2;

        mrkr = '+';
        titlename = 'Even swipes (Zone 3 to 1)';
    elseif i == 2
        %% example 1
        adj(2,4) = 8;
        adj(2,5) = 10;
        adj(2,6) = 10;
        adj(2,7) = 8;
        
        adj(2,3) = 2;

        mrkr = 'o';
        titlename = 'Misdirected food (Zone 2 to 3)';
    elseif i == 3
        %% example 1
        adj(2,4) = 8;
        adj(2,5) = 10;
        adj(2,6) = 10;
        adj(2,7) = 8;
        
        adj(5,3) = 2;

        mrkr = '*';
        titlename = 'Misdirected food (Zone 5 to 3)';
    elseif i == 4
        %% example 1
        adj(2,4) = 8;
        adj(2,5) = 10;
        adj(2,6) = 10;
        adj(2,7) = 8;

        adj(5,8) = 1;
        adj(6,11) = 1;

        adj(1,3) = 1;

        mrkr = 'v';
        titlename = 'Inter-plate swipes (Direct)';
    elseif i == 5
        %% example 1
        adj(2,4) = 8;
        adj(2,5) = 10;
        adj(2,6) = 10;
        adj(2,7) = 8;

        adj(5,4) = 1;
        adj(6,7) = 1;

        adj(1,3) = 1;

        mrkr = 'pentagram';
        titlename = 'Inter-plate swipes (Indirect)';
    elseif i == 6
        %% example 1
        adj(2,4) = 8;
        adj(2,5) = 10;
        adj(2,6) = 10;
        adj(2,7) = 8;

        adj(5,2) = 2;

        mrkr = 'square';
        titlename = 'Return food (Zone 5 to 2)';
    elseif i == 7
        %% example 1
        adj(2,4) = 8;
        adj(2,5) = 10;
        adj(2,6) = 10;
        adj(2,7) = 8;
        
%         adj(3,1) = 2;

        mrkr = 'x';
        titlename = 'Uneven swipes (Zone 3 to 1)';
    end

    weight = .01; %1/length(adj)^2;%0.01;%sum(adj(:))/length(adj)^2;
    [adj] = create_adj(adj,weight);

    pert = 0;%.179;
    for j = 4:7
        adj(j,j) = adj(j,j)-pert;
    end

    [eig_one] = first_eigvec(adj);

    tempadj = adj-ones(16,16)*adj(1,1);%0.01;
    tempadj = tempadj - diag(diag(tempadj));

    hold on
%     scatter(1:16,abs(eig_one),'Marker',mrkr)
    plot(1:16,abs(eig_one),['-',mrkr])

    all_legends = [all_legends,legend_list(i)];

    figure()
    [x,y] = network_layout();
    
    DG = digraph(tempadj);

    LWidths = 10*DG.Edges.Weight/max(DG.Edges.Weight);
    p = plot(DG,'EdgeLabel',round(DG.Edges.Weight,2),'LineWidth',LWidths);
    % To remove all labels
    p.NodeLabel = {};

    p.XData = x;
    p.YData = y;
    hold on
    scatter(x,y,0.01+(abs(eig_one)-min(abs(eig_one)))*1000,'LineWidth',2)

    title(titlename)
    legend('Swipe transition graph',['Eigenvector entry' newline '(deviation from minimum)'])

    axis equal

    saveas(gcf,['Figures/Graph_',num2str(i),'.png'])
    close gcf
end

xlabel('Zone/Vertex','fontsize', 14)
ylabel('Eigenvector entry','fontsize', 14)
legend(all_legends)
title(['$p_i = ',num2str(pert),', w = ',num2str(weight),'$'],'Interpreter','latex')

function [adj] = create_adj(adj,weight)
    adj = adj - diag(diag(adj));
    adj = (adj./sum(adj(:)));% Normalising
    

    temp = ones(length(adj),length(adj))*weight; 
    temp(isnan(temp))=0;
    adj=adj+temp;
    
end

function [eig_one] = first_eigvec(adj) 
    
    [V,D] = eig(adj');
    
    [~,I]=sort(diag(real(D)),'desc');

    eig_one = V(:,I(1));
end

% function [] = scatter
% %     plot(1:16,abs(V(:,I(1))))
%     hold on
%     scatter(1:16,abs(V(:,I(1))),'Marker',mrkr)
%     hold on
% end

function [x,y] = network_layout()

    yl = 768;
    obj(1,1:4)=[19, 227, yl-752, yl-525];
    obj(2,1:4)=[234, 774, yl-678, yl-525];
    obj(3,1:4)=[781, 1004, yl-752, yl-525];
    obj(4,1:4)=[106, 101+126, yl-400-86, yl-400]; 
    obj(5,1:4)=[330,325+126, yl-400-86, yl-400];  
    obj(6,1:4)=[556, 551+126, yl-400-86, yl-400]; 
    obj(7,1:4)=[781, 776+126, yl-400-86, yl-400];  
    obj(12,1:4)=[19, 252, yl-190, yl-74];
    obj(13,1:4)=[257, 420, yl-190, yl-27];
    obj(14,1:4)=[426, 589, yl-190, yl-27];
    obj(15,1:4)=[595, 758, yl-190, yl-27];
    obj(16,1:4)=[770, 1005, yl-190, yl-27];
    obj(8,1:4)=[19, 309, yl-523, yl-196];
    obj(9,1:4)=[311, 531, yl-523, yl-196];
    obj(10,1:4)=[533,753, yl-523, yl-196];
    obj(11,1:4)=[755, 1004, yl-523, yl-196];

    xdata = obj(:,1:2);
    ydata = obj(:,3:4);

    x = mean(xdata,2);
    y = mean(ydata,2);
end