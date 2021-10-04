function [sets] = set_allocate_TYPE_OND_AGE(subject_details,OND_details,nam_save,saved,age)
    if age == 1
        age_range = [2 + 5/12,3 + 9/12];
    elseif age == 2
        age_range = [3 + 8/12,4 + 11/12];
    elseif age == 3
        age_range = [4 + 10/12,6 + 1/12];
    end    
    sets = {[] [] [] [] []};

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        [J]=name_id(OND_details(:,2),subject_details{i,1});
        if max(saved(:,I))>0 & ~isempty(J)  %% removing subjects with no sharing game data
            yr = subject_details{i,2} + subject_details{i,3}/12;
            if  yr > age_range(1) && yr < age_range(2)
                if strcmp(subject_details{i,5},'OND')
                    if contains(OND_details{J,6},'ADHD')
                        sets{1} = [sets{1},I];
                    elseif contains(OND_details{J,6},'Down')
                        sets{2} = [sets{2},I];
                    elseif contains(OND_details{J,6},'Language')
                        sets{3} = [sets{3},I];
                    elseif contains(OND_details{J,6},'Global')
                        sets{4} = [sets{4},I];
                    else
                        sets{4} = [sets{4},I];
                    end
                end
            end
        end
    end
end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end