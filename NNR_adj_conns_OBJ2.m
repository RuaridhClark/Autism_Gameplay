function [adj] = NNR_adj_conns_OBJ2(adj,bweight)
    %% create all to all connections
    temp = ones(length(adj),length(adj))*bweight; 
    temp=temp-diag(diag(temp));
    
%     temp(temp==2)=1;
%     for i = 1 : length(temp)
%         temp(i,:)=temp(i,:)./sum(temp(i,:));
%     end    
    temp(isnan(temp))=0;
    adj=adj+temp;
%     adj=adj-diag(diag(adj));

end