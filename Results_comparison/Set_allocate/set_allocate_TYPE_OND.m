function [sets,other] = set_allocate_TYPE_OND(subject_details,name_save,saved)
    load('OND_details.mat','OND_details')

    sets = {[] [] [] []};
    other = [];

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},name_save);
        [J]=name_id(OND_details(:,1),subject_details{i,1});
        if max(saved(:,i))>0 && ~isempty(J)  % removing subjects with no sharing game data
            if strcmp(subject_details{i,5},'OND')
                if contains(OND_details{J,2},'ADHD')
                    sets{1} = [sets{1},I];
                elseif contains(OND_details{J,2},'Down')
                    sets{4} = [sets{4},I];
                elseif contains(OND_details{J,2},'Language')
                    sets{2} = [sets{2},I];
                else
                    sets{3} = [sets{3},I];
                    other = [other,{OND_details{J,2}}];
                end
            end
        end
    end
end

function [I] = name_id(name,name_save)
    I = find(strcmp(name_save,name));
end