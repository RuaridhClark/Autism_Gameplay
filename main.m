% main pipeline

% Ranking_Main folder contains the main analysis
% Ranking_Correlations folder contains the correlation analysis

% adj_obj_end contains the adjacency matrices for swipe start to end links
% OBJ_end.mat contains the paper data swipe end locations

% updated swipe zones:
% Ranking_Main_11072 folder
% Ranking_Correlations_110721 folder
% adjs_110721 folder

% run fix_tapCount2
addpath([pwd,'\Track_Finger\Additional_Functions'])
fix_tapCount_errors

% identify swipes - run finger_swipes3
% create adj - create_adjs_end.m
% Rank subjects - subject_rank_CDI_...

