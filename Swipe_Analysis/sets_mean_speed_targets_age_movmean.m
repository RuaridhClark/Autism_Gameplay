clear all
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\paths.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Track_Finger\Data\swipe_speeds.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\subject_details.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\swipes_all704.mat','nam_save')
folder1='H:\My Documents\MATLAB\Autism_MAIN\Set_allocate';
addpath(folder1)

trgt_speed = zeros(17,17);
trgt_swipes = zeros(17,17);
%%
[months] = find_ages(subject_details_776,nam_save);

%%
[sets] = set_allocate(subject_details_776,nam_save,saved);
exclude = find(ranked==0.5);
exclude = [exclude;find(months>75)];
for i = 1 : length(sets)
    rmv=find(ismember(sets{i},exclude')==1);
    sets{i}(rmv)=[];
end

%%
mn_speed_age = cell(max(months)-min(months)+1,1); 
speed_age = cell(max(months)-min(months)+1,1);
n_swipes_age = cell(max(months)-min(months)+1,1);
for i = 1 : length(speed_age)
    speed_age{i}=zeros(17,17);
    n_swipes_age{i}=zeros(17,17);
end



for m = 1 : 704 % collate speeds and number of swipes for all in age sets
    trgt_speed=speed_age{months(m)-min(months)+1};
    trgt_swipes=n_swipes_age{months(m)-min(months)+1};
    for i = 1 : length(speed{m})
        path = all_paths{m};
        if speed{m}(i,1)>0 && min(path(i,:))>0 && ~isinf(speed{m}(i,1)) 
            trgt_speed(path(i,1),path(i,2)) = trgt_speed(path(i,1),path(i,2)) + speed{m}(i,1);
            trgt_swipes(path(i,1),path(i,2))= trgt_swipes(path(i,1),path(i,2)) + 1;
        end
    end
    speed_age{months(m)-min(months)+1}=trgt_speed;
    n_swipes_age{months(m)-min(months)+1}=trgt_swipes;
end

rng = 2;
for i = 1 : length(speed_age)
    trgt_speed = zeros(17,17);
    trgt_swipes = zeros(17,17);
    
    for j = i-rng : i+rng
        if j >= 1 && j <= length(speed_age)
            trgt_speed = trgt_speed + speed_age{j};
            trgt_swipes = trgt_swipes + n_swipes_age{j};
        end
    end
    
    mn_speed = trgt_speed./trgt_swipes;
    mn_speed(isnan(mn_speed))=0;
    
    if max(mn_speed(:))>0
        mn_speed_age{i} = mn_speed;
    else
        mn_speed_age{i} = [];
    end
end

function [months] = find_ages(subject_details,nam_save)

    for i = 1 : size(subject_details,1)
        [I] = name_id(subject_details{i,1},nam_save);
%         if max(saved(:,I))>0  %% removing subjects with no sharing game data
            yr = subject_details{i,2} + subject_details{i,3}/12;
            mnths = subject_details{i,2}*12 + subject_details{i,3};
            months(I)=mnths;
%         end
    end

end

function [I] = name_id(name,nam_save)
    I = find(strcmp(nam_save,name));
end