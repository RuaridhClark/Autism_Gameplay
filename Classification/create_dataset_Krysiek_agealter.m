dataset_Krysiek = create_data();

function dataset_Krysiek = create_data()
    num = 12;               % accurate = 16, snap-to = 12
    option = 1;             % 1 = n_swipes, 2 = sharing score, 3 = swipe accuracy ratio
    destination = 'plates';   % n_swipes for 'plates', 'food' or 'inter' (inter-plates) destinations
    gender = '';
    severity = '';
    
    [folder_loc,alt_folder_loc,file_loc,floc] = setup();
    tab_sev = readtable([floc,'\eCRF.csv']);
    
    addpath 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_krysiek'
    load('subject_details_combine_ond.mat')
    subject_details = subject_details_combine;
    
    extra = load([folder_loc,'\Ranking_Correlations\Data\snapto_2only_Krysiek.mat'],'nam_save','ranked');
    sharing_score_deliv = extra.ranked;
    nam_save = extra.nam_save;
    
    extra = load([folder_loc,'\Ranking_Correlations\Data\accurate_2only_Krysiek.mat'],'nam_save','ranked');
    sharing_score_accurate = extra.ranked;
    
    extra = load([folder_loc,'\Ranking_Correlations\Data\OBJ_snapto_redirect_Krysiek2.mat'],'nam_save','ranked'); % inter-plate sharing score
    sharing_score_inter = extra.ranked;
    
    %% no. of swipes
    [n_sw_food_deliv,~] = swipe_analysis(num,file_loc,nam_save,'plates');
    [n_sw_accurate,~] = swipe_analysis(16,file_loc,nam_save,'plates');
    [n_sw_food_zone,~] = swipe_analysis(num,file_loc,nam_save,'food');
    [n_sw_inter,~] = swipe_analysis(16,file_loc,nam_save,'inter');
    [n_sw_total,list] = swipe_analysis(num,file_loc,nam_save,'total');
    
    saved=ones(1,length(subject_details)); saved(list)=zeros(1,length(list)); % create list of 1s delete those without adjs
    [sets,months] = create_sets_months(subject_details,nam_save,saved,gender,tab_sev,severity);
    
    cat = []; age = []; gender = []; n_sw_d = []; n_sw_a = []; n_sw_f =[]; n_sw_i = []; n_sw_t = []; ss_d = []; ss_a = []; ss_i = [];
    for i = 1 : size(nam_save,2)
        [I] = find(strcmp(nam_save{1,i},subject_details(:,1)));
    
        if max(saved(:,I))>0 && ~strcmp(subject_details{I,5},'ONDE') %&& ~strcmp(subject_details{I,5},'OND') %&& strcmp(subject_details{I,4},'Male')
            age = [age;subject_details{I,2}*12 + subject_details{I,3}];
            cat = [cat;{subject_details{I,5}}];
            gender = [gender;{subject_details{I,4}}];
            
            m=1.221; b=4.357;
            n_sw_d = [n_sw_d;n_sw_food_deliv(i)/(m*age(end)+b)];
            m=1.236; b=-9.371;
            n_sw_a = [n_sw_a;n_sw_accurate(i)/(m*age(end)+b)];
            m=-0.1905; b=52.814;
            n_sw_f = [n_sw_f;n_sw_food_zone(i)/(m*age(end)+b)];
            m=0.0410; b=-0.5786;
            n_sw_i = [n_sw_i;n_sw_inter(i)/(m*age(end)+b)];
            m=0.5879; b=104.37;
            n_sw_t = [n_sw_t;n_sw_total(i)/(m*age(end)+b)];
            m=0.0023; b=-0.05786;
            ss_d = [ss_d;sharing_score_deliv(i)/(m*age(end)+b)];
            m=0.0025; b=-0.118;
            ss_a = [ss_a;sharing_score_accurate(i)/(m*age(end)+b)];
            m=0.00186; b=-0.0097;
            ss_i = [ss_i;sharing_score_inter(i)/(m*age(end)+b)];
        end
    end
    
    dataset_Krysiek = table(cat,age,gender,n_sw_d,n_sw_a,n_sw_f,n_sw_i,n_sw_t,ss_d,ss_a,ss_i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [folder_loc,alt_folder_loc,file_loc,floc] = setup()
    folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay';
    alt_folder_loc = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data';
    file_loc = [folder_loc,'\adjs\adj_obj_end_accurate\']; % should match zone type
    floc=[alt_folder_loc,'\IQ_severity'];

    folder1 = [folder_loc,'\Set_allocate'];
    folder2 = [folder_loc,'\Plots'];
    folder3 = [folder_loc,'\Create_adj'];
    folder4 = [folder_loc,'\adjs\adj_obj_end_accurate'];
    folder5 = [folder_loc,'\Data'];
    folder6 = folder_loc;
    addpath(folder1,folder2,folder3,folder4,folder5,folder6)
end

function [nam_save,saved,ranked,list] = load_dataset(option,num,folder_loc,destination)
    saved = []; list = [];
    if num == 16
        if strcmp(destination,'inter')
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_accurate.mat'],'nam_save','ranked') % inter-plate sharing score
        else
            load([folder_loc,'\Ranking_Correlations\Data\accurate_2only.mat'],'nam_save','ranked')
        end
    elseif num == 12
        if strcmp(destination,'inter')
            load([folder_loc,'\Ranking_Correlations\Data\OBJ_snapto_redirect2.mat'],'nam_save','ranked') % inter-plate sharing score
        else
            load([folder_loc,'\Ranking_Correlations\Data\snapto_2only.mat'],'nam_save','ranked')
        end
    end
%     end
end

function [sets,months] = create_sets_months(subject_details,nam_save,saved,gender,tab_sev,severity)
    [months] = list_AGE(subject_details,nam_save,saved);
    if size(months,1)<size(months,2)
        months=months';
    end
    if strcmp(severity,'on')
        [sets] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,gender);
    elseif strcmp(gender,'male')
        [tmp_sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
        sets = tmp_sets(5:8);
    elseif strcmp(gender,'female')
        [tmp_sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
        sets = tmp_sets(1:4);
    elseif strcmp(gender,'compare')
        [sets] = set_allocate_GENDER_TYPE(subject_details,nam_save,saved);
    else
        [sets] = set_allocate(subject_details,nam_save,saved);
%         [sets_sev] = set_allocate_severity(subject_details,nam_save,saved,tab_sev,gender);
%         sets{2} = setdiff(sets{2},sets_sev{4});
    end
end

function [n_swipes,list] = swipe_analysis(num,file_loc,nam_save,destination) 
    n_swipes = zeros(1,length(nam_save));
    f_num = 0;
    
    list=[];
    for jj = 1:length(nam_save)
        adj=zeros(num,num);
        file_id = ['subject_',nam_save{jj},'.mat'];
    
        if isfile([file_loc,file_id]) || isfile(['C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\GitHub\Autism_Gameplay\adjs\adj_krysiek\',file_id])
            f_num = f_num + 1;
            load(file_id,'adj')
            if num == 12
                [adj] = adj_snap2zones(adj,num);
            end
%             adj=adj(1:num,1:num);
        end
    
        if max(adj(:))==0
            list=[list,jj];
        end
    
        if strcmp(destination,'plates')
            n_swipes(jj) = sum(adj(2,[4,5,6,7]));
        elseif strcmp(destination,'food')
            n_swipes(jj) = sum(adj(2,2));
%             n_swipes(jj) = sum(adj(2,:));
%             temp=diag(adj);
% %             sv_2 = temp(2);
%             temp(2)=0;
%             n_swipes(jj) = sum(temp);
        elseif strcmp(destination,'inter')
            n_swipes(jj) = sum(adj(4,[5,6,7]))+sum(adj(5,[4,6,7]))+sum(adj(6,[4,5,7]))+sum(adj(7,[4,5,6]));
        elseif strcmp(destination,'total')
            n_swipes(jj) = sum(adj(:));
        end
    end
end

function [adj] = adj_snap2zones(adj,num)
    % rewire adjacency from 16 to 12 zones
    init=adj;
%     allow=[2,4,5,6,7];
    for it = 2
        adj(it,4:7)=adj(it,4:7)+adj(it,13:16);    % reconnect to 4-7
        adj(it,13:16)=zeros(1,4);                     % remove non-food connections
    end
%     adj = adj(1:12,1:12);
    %% remove zn 4-7 incoming except from 2
    allow=[2];%,4,5,6,7];
    for it = 1:num
        if ~ismember(it,allow)
            adj(it,13:16) = adj(it,4:7)+adj(it,13:16);
            adj(it,4:7)=zeros(1,4);                     % remove non-food connections
        end
    end
%     adj=adj-diag(diag(adj));
    
    bweight=1;
%     [adj] = NNR_adj_conns_OBJ2(adj,bweight);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
