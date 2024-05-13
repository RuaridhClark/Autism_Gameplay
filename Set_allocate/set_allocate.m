function [sets] = set_allocate(subject_details,nam_save,saved)
    sets = {[] [] []};

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        if max(saved(:,i))>0 && subject_details{i,2} < 7 % removing subjects with no sharing game data
            if strcmp(subject_details{i,5},'TD')
                sets{1} = [sets{1},I];
            elseif strcmp(subject_details{i,5},'ASD')
                sets{2} = [sets{2},I];
            elseif strcmp(subject_details{i,5},'OND')
                sets{3} = [sets{3},I];
            end
        end
    end
end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end