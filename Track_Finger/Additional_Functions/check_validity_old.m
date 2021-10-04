function [valid,invalid,starts,ends] = check_validity_old(tab,starts,ends)
    %% Check validity of touch recording
    v=1;    % swipe ID
    invalid=[];
    for j = 2 : length(tab.X)           % for each touch recording
        if (tab.X(j-1)-tab.X(j))== 0 && (tab.Y(j-1)-tab.Y(j))==0 ...
            && abs(tab.TouchPhase(j)-tab.TouchPhase(j-1))~=3 % Check for repeats and only allow them if a touch is ending
            invalid=[invalid,j];        % add touch to invalid list
            if ismember(j,starts)
                starts(starts==j)=[];   % remove from list of starts
            end
            if ismember(j,ends)
                ends(ends==j)=[];       % remove from list of ends
            end
        end
    end
    
    valid=1:length(tab.X);
    valid(ismember(valid,invalid))=[];  % remove all invalid touches from valid list
end