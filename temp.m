figure;
for i = 1:length(sets{1,4})
swipes=zone2_swipes{sets{1,4}(i)};
hold on
plot(1:length(swipes),swipes)
end

load('C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\Create_adj\zone2_swipes_combine2.mat')
for i = 1 : 1097
if ~isempty(zone2_swipes{i})
n_swipes(i) = length(find(zone2_swipes{i}>200));
else
n_swipes(i)=0;
end
end