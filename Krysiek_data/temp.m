for Ind_dev = 1 : 862
    if contains(device.source(Ind_dev),'PILOT') 
        [cat] = neurodev_category(device.disorder_type{Ind_dev},device.comment{Ind_dev},device.disorder(Ind_dev));
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