function [months] = list_AGE(subject_details,nam_save,saved)

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
        if max(saved(:,I))>0  %% removing subjects with no sharing game data
            yr = subject_details{i,2} + subject_details{i,3}/12;
            mnths = subject_details{i,2}*12 + subject_details{i,3};
            months(I)=mnths;
%         else
%             1
        end
    end

end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end