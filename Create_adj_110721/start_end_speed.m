% close all
clear all
file_loc = 'I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\';
load('swipes_all704.mat')
load('H:\My Documents\MATLAB\Autism_MAIN\Track_Finger\Data\swipe_speeds.mat')

file_adj = 'H:\My Documents\MATLAB\Autism_MAIN\adjs\';

yl = 768;
obj(1,1:4)=[69, 69+150, yl-527-125, yl-527];
minus(1,1:4)=[144,144+75,yl-527-70,yl-527];
obj(2,1:4)=[270, 270+474, yl-526-170, yl-526];
obj(3,1:4)=[790, 790+150, yl-540-112, yl-540];
minus(3,1:4)=[790,790+75,yl-540-71,yl-540];
obj(4,1:4)=[825, 825+158, yl-684-84, yl-684];
obj(5,1:4)=[101, 101+131, yl-400-96, yl-400];
obj(6,1:4)=[325,325+131, yl-400-96, yl-400];
obj(7,1:4)=[551, 551+131, yl-400-96, yl-400];
obj(8,1:4)=[776, 776+131, yl-400-96, yl-400];
obj(9,1:4)=[113, 113+112, yl-283-84, yl-283];
obj(10,1:4)=[340, 340+112, yl-283-84, yl-283];
obj(11,1:4)=[566, 566+112, yl-283-84, yl-283];
obj(12,1:4)=[790, 790+112, yl-283-84, yl-283];
obj(13,1:4)=[21, 21+234, yl-28-208, yl-28];
obj(14,1:4)=[260, 260+158, yl-28-140, yl-28];
obj(15,1:4)=[428, 428+158, yl-28-140, yl-28];
obj(16,1:4)=[597, 597+158, yl-28-140, yl-28];
obj(17,1:4)=[772, 772+234, yl-54-208, yl-54];

all_paths = cell(1,704);
for i = 1:704
    close all
    swipe=swipe_save{i};
    skip=0;
    
    file_id = [nam_save{i},'/',nam_save{i},'.Sharing.TouchData.typed.csv'];
    
    if isfile([file_loc,file_id])
        sw_speed = speed{i};
        mn_speed = mean(sw_speed(sw_speed>0));
        
        colr = 'b';
        titlename = ['ID ',nam_save{i}];
        filename = [file_loc,file_id];
        savename = ['subject_',nam_save{i}];
        
        tab=sortrows(readtable(filename),4); % sort table according to time
        if iscell(tab.X(1)) %% Convert strings to doubles for X and Y
            tab=sortrows(readtable(filename),14); % sort table according to time
            tab.X=strrep(tab.X,',','.');
            tab.Y=strrep(tab.Y,',','.');
            tab.X=str2double(tab.X);
            tab.Y=str2double(tab.Y);
        end
              
        prev_tP =0;
        path = zeros(length(swipe),2);
        for m = 1 : length(swipe)
            prev_n = [];
            n=1;
            while isempty(prev_n) && n<=length(swipe{m})

                k = swipe{m}(n);
                nn=[];
                for ii = 1 : 17
                    if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                        if ii == 1 || ii == 3   % zone 1 and 3 have areas to be minused from shape
                            if minus(ii,1)<tab.X(k) && tab.X(k)<minus(ii,2) && minus(ii,3)<tab.Y(k) && tab.Y(k)<minus(ii,4)
                            else
                                nn=ii;
                            end
                        else
                            nn = ii; 
                        end
                    end
                end

                if n>1 && ~isempty(prev_n) && ~isempty(nn)

                elseif isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
                n=n+1;
            end

            %% start from the back
            nn=[];
            n=-1;
            while isempty(nn) && ~isempty(prev_n)
                n=n+1;
                k = swipe{m}(end-n);

                for ii = 1 : 17
                    if obj(ii,1)<tab.X(k) && tab.X(k)<obj(ii,2) && obj(ii,3)<tab.Y(k) && tab.Y(k)<obj(ii,4)
                        if ii == 1 || ii == 3   % zone 1 and 3 have areas to be minused from shape
                            if minus(ii,1)<tab.X(k) && tab.X(k)<minus(ii,2) && minus(ii,3)<tab.Y(k) && tab.Y(k)<minus(ii,4)
                            else
                                nn=ii;
                            end
                        else
                            nn = ii; 
                        end
                    end
                end

                if ~isempty(prev_n) && ~isempty(nn)
                    path(m,:)=[prev_n,nn];
                elseif isempty(prev_n) && ~isempty(nn)
                    prev_n = nn;
                end
            end
        end
    end
    
    all_paths{i} = path;
end