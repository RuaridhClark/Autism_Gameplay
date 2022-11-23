function [sets,other] = set_allocate_TYPE_OND(subject_details,OND_details,nam_save,saved)
    sets = {[] [] [] []};
    other = [];
    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        [J]=name_id(OND_details(:,2),subject_details{i,1});
        if max(saved(:,I))>0 & ~isempty(J)  %% removing subjects with no sharing game data
            if strcmp(subject_details{i,5},'OND')
                if contains(OND_details{J,6},'ADHD')
                    sets{1} = [sets{1},I];
                elseif contains(OND_details{J,6},'Down')
                    sets{2} = [sets{2},I];
                elseif contains(OND_details{J,6},'Language')
                    sets{3} = [sets{3},I];
%                 elseif contains(OND_details{J,6},'Language') && contains(OND_details{J,6},'delay')
%                     sets{5} = [sets{5},I];
%                 elseif contains(OND_details{J,6},'Language') && contains(OND_details{J,6},'disorder')
%                     sets{6} = [sets{6},I];
%                 elseif contains(OND_details{J,6},'Global')
%                     sets{4} = [sets{4},I];
                else
%                     [OND_details{J,6},OND_details{J,7}]
                    sets{4} = [sets{4},I];
                    other = [other,{OND_details{J,6}}];
                end
                
            end
        end
    end
end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end