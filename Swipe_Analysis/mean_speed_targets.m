clear all
load('H:\My Documents\MATLAB\Autism_MAIN\Create_adj\paths.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Track_Finger\Data\swipe_speeds.mat')

trgt_speed = zeros(17,17);
trgt_swipes = zeros(17,17);

for m = 1 : 704
    for i = 1 : length(speed{m})
        path = all_paths{m};
        if speed{m}(i,1)>0 && min(path(i,:))>0 && ~isinf(speed{m}(i,1))
            trgt_speed(path(i,1),path(i,2)) = trgt_speed(path(i,1),path(i,2)) + speed{m}(i,1);
            trgt_swipes(path(i,1),path(i,2))= trgt_swipes(path(i,1),path(i,2)) + 1;
        end
    end
end

mn_speed = trgt_speed./trgt_swipes;
mn_speed(isnan(mn_speed))=0;