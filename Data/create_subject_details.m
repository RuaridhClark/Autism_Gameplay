% create subject_details

clear all
file_init = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\Krysiek_data\subject_data\';
D = dir(file_init);
subject_details = cell(length(D),5);
i = 0;
for f = 1:length(D)    
    file_id = D(f).name; % Get the current subdirectory name
    file_loc =[file_init,file_id];
    
    if isfile(file_loc)
        load(file_loc)
        i = i +1;
        subject_details{i,1} = erase(file_id,'.mat');
        subject_details{i,2} = floor(age/12);
        subject_details{i,3} = rem(age,12);
        if gender == 'f'
            subject_details{i,4} = 'Female';
        elseif gender == 'm'
            subject_details{i,4} = 'Male';
        end
        if strcmp(cat,'asd')
            subject_details{i,5} = 'ASD';
        elseif strcmp(cat,'td')
            subject_details{i,5} = 'TD';
        elseif strcmp(cat,'ond')
            subject_details{i,5} = 'OND';
        end
    end
end