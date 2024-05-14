function [months] = list_AGE(subject_details,name_save,saved)

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},name_save);
        if max(saved(:,I))>0  % removing subjects with no sharing game data
            yr = subject_details{i,2} + subject_details{i,3}/12; % age - years
            mnths = subject_details{i,2}*12 + subject_details{i,3}; % age - months
            months(I)=mnths;
        end
    end

end

function [I] = name_id(name,name_save)
    I = find(strcmp(name_save,name));
end