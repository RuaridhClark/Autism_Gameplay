function [sets] = set_allocate(subject_details,nam_save,saved)
    sets = {[] [] [] []};
%     OND = [];
    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        if max(saved(:,I))>0 & subject_details{i,2} < 7 %& strcmp(subject_details{i,4},'Male')%% removing subjects with no sharing game data
            if strcmp(subject_details{i,5},'TD')
                sets{1} = [sets{1},I];
            elseif strcmp(subject_details{i,5},'ASD')
                sets{2} = [sets{2},I];
            elseif strcmp(subject_details{i,5},'OND')
                if contains(nam_save,'Study3')
                    1
                end
                sets{3} = [sets{3},I];
%                 OND = [OND;subject_details{i,1}];
            elseif strcmp(subject_details{i,5},'ONDE')
                sets{4} = [sets{4},I];
            end
        end
    end

end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end