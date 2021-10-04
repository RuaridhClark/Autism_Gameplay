function [sets] = set_allocate_GENDER_AGE(subject_details,nam_save)
    sets = {[] [] [] [] [] [] [] []};
    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        if max(saved(:,I))>0  %% removing subjects with no sharing game data
            yr = subject_details{i,2} + subject_details{i,3}/12;
            if  yr > (4 + 10/12) && yr < (6 + 1/12)
                if strcmp(subject_details{i,4},'Female')
                    if strcmp(subject_details{i,5},'TD')
                        sets{1} = [sets{1},I];
                    elseif strcmp(subject_details{i,5},'ASD')
                        sets{2} = [sets{2},I];
                    elseif strcmp(subject_details{i,5},'OND')
                        sets{3} = [sets{3},I];
                    elseif strcmp(subject_details{i,5},'ONDE')
                        sets{4} = [sets{4},I];
                    end
                elseif strcmp(subject_details{i,4},'Male')
                    if strcmp(subject_details{i,5},'TD')
                        sets{5} = [sets{5},I];
                    elseif strcmp(subject_details{i,5},'ASD')
                        sets{6} = [sets{6},I];
                    elseif strcmp(subject_details{i,5},'OND')
                        sets{7} = [sets{7},I];
                    elseif strcmp(subject_details{i,5},'ONDE')
                        sets{8} = [sets{8},I];
                    end
                end
            end
        end
    end

end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end