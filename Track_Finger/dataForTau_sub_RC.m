close all
clear all
clc

format long

%%

xy2mm = 200.66/5*3/768;

xmin  = 0.5; xmax  = 1024;
ymin  = 0.5; ymax  = 768;

start_x1  = 330-193/2+15;
start_x2  = 676+192/2-15;
start_y1  = 121;
start_y2  = 299-67/2-5;

end_x1 = 173-111;
end_x2 = 841+111;
end_y1 = 299-67/2+8;
end_y2 = 520;
        
STARTbox  = [start_x1, start_y2;...
             start_x1, start_y1;...
             start_x2, start_y1;...
             start_x2, start_y2;...
             start_x1, start_y2];

ENDbox = [end_x1, end_y1;...
          end_x1, end_y2;...
          end_x2, end_y2;...
          end_x2, end_y1;...
          end_x1, end_y1];

% plate3 = [618, 299];
% watermelon4 = [600, 184+62/2];

%%
% [file,path] = uigetfile('*.csv');
% filedir = [path, file];
% filedir ='/Users/xnb18110/Desktop/Projects/tau/data_study0_sharingTouch/POL-AA-001-sharingTouch.csv';
filedir ='I:\Engineering\EEE\RESEARCH\SPACE\MALCOLMSPACE\2013_RuaridhClark\Research\Project\Autism\PlayCare\PlayCare\GLA-AM-001\GLA-AM-001.Sharing.TouchData.typed.csv';
% data = importdata(filedir,',');

% touch1 = [];
touch2 = [];
X = [];
Y = [];
tp = [];
T = [];
tp0 = [];
tp3 = [];
mt1 = [];
mt2 = [];
mt_temp = [];
multit = [];
mt3 = [];
mt4 = [];
mt5 = [];
mt6 = [];
mt7 = [];
mt8 = [];
food2pltchi =[];
     
%% IDENTIFY SWIPES

% touch1 = data.data;
% touch2 = sortrows(touch1,1); % sort by time
touch2=sortrows(readtable(filedir),'Time'); % sort table according to time
if iscell(touch2.X(1)) %% Convert strings to doubles for X and Y
    touch2.X=str2double(strrep(touch2.X,',','.'));
    touch2.Y=str2double(strrep(touch2.Y,',','.'));
end

 X = touch2.X;
 Y = touch2.Y;
tp = touch2.TouchPhase; % touch phase
 T = touch2.Time;

tp0 = find(tp==0); % find the start of a movement
for i = 1:length(tp0)-1
    if isempty(find(tp(tp0(i):tp0(i+1),:)==3)) == 1 % does this swipe have NO recorded end
        tp3(i,1) = 0;
    else
        tp3(i,1) = find(tp(tp0(i):tp0(i+1),:)==3,1)+(tp0(i)-1); % tp3 is location of end point in touch data
    end
end

%%% note: checking last swipe in the same way as before
if isempty(find(tp(tp0(end):end) ==3)) == 1
    tp3(length(tp0),1) = 0;
else
    tp3(length(tp0),1) = find(tp(tp0(end):end)==3,1)+(tp0(length(tp0))-1);
end

mt1 = [tp0 tp3];
mt2 = mt1;
mt2(find(tp3==0),:)=[]; % remove those without ending point

celltouch2=table2cell(touch2);
for i = 1:length(mt2(:,1))
    mt_temp = touch2(mt2(i,1):mt2(i,2),:);
    if isempty(find(diff(mt_temp.Time)==0)) == 1 % check diff between times to see if multi touch occured
        multit(i,1)=1;
    else
        multit(i,1)=0; % find multiple touch
    end
end

mt3 = mt2;
mt3(find(multit==0),:) = []; % remove multiple touch 

mt4(:,1:2) = mt3; % start-frame and end-frame of a movement 
mt4(:,3) = mt3(:,2)-mt3(:,1)+1; % total frame number of a movement % 2frames:tap % >2frames:swipe
mt5 = mt4;
mt6 = mt5(find(mt4(:,3)>4),:); % find movement >4frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mt7 = mt6;
mt7(:,4) = X(mt6(:,1)); % startX
mt7(:,5) = Y(mt6(:,1)); % startY
mt7(:,6) = X(mt6(:,2)); %   endX
mt7(:,7) = Y(mt6(:,2)); %   endY

% make sure the swipes start and end within the right zones (food and
% plates/children)
food2pltchi = find(mt7(:,4)>start_x1 & mt7(:,4)<start_x2 &...% start point x boundary
                   mt7(:,5)>start_y1 & mt7(:,5)<start_y2 &...% start point y boundary
                   mt7(:,6)>  end_x1 & mt7(:,6)< end_x2 &... %  end point x boundary
                   mt7(:,7)>  end_y1 & mt7(:,7)< end_y2);    %  end point y boundary

mt8 = mt7(food2pltchi,:); % exclude swipes that don't start and end in these zones

%% 
if length(mt8(:,1))>0 % F2PC

    r0 = []; % r0 = [0;0;0;0];
    r  = []; % r  = [0;0;0;0];
    
    xa = [];
    ya = [];
    
    for i = 1:length(mt8(:,1))  % one swipe at a time
        
        t0 = [];
        x0 = [];
        y0 = [];
        pe = [];
        pj = [];
        rj = [];
        
        t0 = T(mt8(i,1):mt8(i,2),1);
        x0 = X(mt8(i,1):mt8(i,2),1)*xy2mm;  % convert xy coords to mm
        y0 = Y(mt8(i,1):mt8(i,2),1)*xy2mm;
        
        for j = 1:length(x0)
            pe = [x0(end) y0(end)]; % last position
            pj = [x0(j)   y0(j)  ]; % current position
            rj(j,1) = norm(pe-pj);  % norm of last - current vectors
        end
        r0 = [r0;rj];
        
        t = [];
        x = [];
        y = [];
        pe = [];
        pk = [];
        rk = [];
        
        fs = 60;
        [x,t] = resample(x0,t0,fs,1,1); % resample at a constant rate
        [y,t] = resample(y0,t0,fs,1,1);
        
        xa = [xa; x];
        ya = [ya; y];
        
        for k = 1:length(x)
            pe = [x(end) y(end)];   % last position
            pk = [x(k)   y(k)  ];   % current position
            rk(k,1) = norm(pe-pk);
        end
        r = [r;rk];
        
        %%% SUB-PROGRAM
        % dataForTau_sub2_plotOneSwipe
        
        % dlmwrite([subject '_' group '.txt'],r)
        
        % csvwrite([filename(1:23) '-x.txt'],xa)
        % csvwrite([filename(1:23) '-y.txt'],ya)
        % csvwrite([filename(1:23) '-r.txt'],r)
    end
     
else
end
