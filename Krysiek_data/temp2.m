% close all
% clear all

file_device = 'device_202106201434.csv';
file_session = 'session_202106201636.csv';
file_touch = 'touchdata_202106201443.csv';

% touch contains session ID
% session references session ID and device ID
% device links to disorder,age,gender

device = readtable(['Data\',file_device]);
session = readtable(['Data\',file_session]);

list = [];

for Ind_sess = 1 : 2204
    if strcmp(session.app_name{Ind_sess},'SharingFoodGame')
        device_id = session.device_id{Ind_sess};
        Ind_dev = find(strcmp(device_id,device.device_id));
        [cat] = neurodev_category(device.disorder_type{Ind_dev},device.comment{Ind_dev},device.disorder(Ind_dev));
        if strcmp('ond',cat)
            if device.age(Ind_dev)>=30 && device.age(Ind_dev)<=72
                list = [list;[device.dev_id(Ind_dev),session.sess_id(Ind_dev),device.source(Ind_dev),device.disorder_type{Ind_dev},device.comment{Ind_dev},device.disorder(Ind_dev),device.age(Ind_dev)]];
            end
        end
    end
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