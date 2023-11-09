% close all
% clear all

file_device = 'device_202106201434.csv';
file_session = 'session_202106201636.csv';
file_touch = 'touchdata_202106201443.csv';

% touch contains session ID
% session references session ID and device ID
% device links to disorder,age,gender

% device = readtable(['Data\',file_device]);
% session = readtable(['Data\',file_session]);
% touch = readtable(['Data\',file_touch]);

i=1;
while i <= height(touch)
    session_id = touch.session_id{i};
    Ind_sess = find(strcmp(session_id,session.session_id));
    touch_list = find(strcmp(session_id,touch.session_id));
    if ~isempty(Ind_sess)
        if strcmp(session.app_name{Ind_sess},'SharingFoodGame')
            device_id = session.device_id{Ind_sess};
            Ind_dev = find(strcmp(device_id,device.device_id));
            if contains(device.source(Ind_dev),'PILOT') %device.dev_id(Ind_dev)<1067 % data beyond not categorised %%%%%%%%%%%%%%%%%%%% why PILOT?
                [cat] = neurodev_category(device.disorder_type{Ind_dev},device.comment{Ind_dev},device.disorder(Ind_dev));
                if ~strcmp(cat,'other')
%                 if strcmp(cat,'ond') 
                    gender = lower(device.gender{Ind_dev});
                    age = device.age(Ind_dev);
                    
                    subject_table = touch(touch_list,:);
                    
                    subject_name = ['subject_data/',device.source{Ind_dev},'-',num2str(device.dev_id(Ind_dev)),'-',num2str(session.sess_id(Ind_sess))];
                    save(subject_name,'subject_table','age',"gender","cat",'device_id')
                end
            end
        end
    end
    i = touch_list(end)+1;% find next subject
end



function [cat] = neurodev_category(in1,in2,disorder)
    input1 = lower(in1);
    input2 = lower(in2);
    if disorder == 1 || contains(input1,'asd') || contains(input1,'autyzm') || contains(input1,'asperger') || contains(input1,'asprgera') || contains(input2,'asd') || contains(input2,'autyzm') || contains(input2,'asperger') || contains(input2,'asprgera')
        cat = 'asd';
    elseif disorder == 0
        cat = 'td';
    elseif contains(input1,'intellectual disability') || contains(input1,'other') || contains(input2,'down syndrom') || contains(input2,'zespol downa') || contains(input2,'Ca?o?ciowe zaburzenia')
        cat = 'ond';
    elseif disorder == 2
        cat = 'ond2';
    elseif disorder == 3
        cat = 'ond3';
    elseif disorder == 4
        cat = 'ond4';
    else
        cat = 'other';
    end
end