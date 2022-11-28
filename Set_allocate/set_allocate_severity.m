function [sets] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,gender)
    sets = {[] [] [] []};
    OND = [];
    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        [Ind] = name_id(subject_details{i,1},tab_sev.id_study_id(:));
        if ~isempty(Ind)
            [sev_num] = severity_score(tab_sev.clinical_diagnosis__asd_severity_level{Ind});
        else
            sev_num = 0;
        end
%         if ~isempty(Ind)
        if strcmp(gender,subject_details{i,4}) | strcmp(gender,'')
            if max(saved(:,i))>0 & subject_details{i,2} < 7 %& strcmp(subject_details{i,4},'Male')%% removing subjects with no sharing game data
                if strcmp(subject_details{i,5},'TD')
                    sets{1} = [sets{1},I];
                elseif strcmp(subject_details{i,5},'ASD')
                    if sev_num == 1
                        sets{2} = [sets{2},I];
                    elseif sev_num == 2
                        sets{3} = [sets{3},I];
                    elseif sev_num == 3
                        sets{4} = [sets{4},I];
                    end
                end
            end
        end
%         end
    end

end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end

function [sev_num] = severity_score(severity)
    if strcmp(severity,'1. Level 1 "Requiring support"')
        sev_num = 1;
%         sets{1} = [sets{1},m];
    elseif strcmp(severity,' 2. Level 2 "Requiring substantial support"')
        sev_num = 2;
%         sets{2} = [sets{2},m];
    elseif strcmp(severity,'3. Level 3 "Requiring very substantial support"')
        sev_num = 3;
%         sets{3} = [sets{3},m];
    else
        sev_num = '';
    end
end