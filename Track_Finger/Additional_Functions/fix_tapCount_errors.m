% check for mismatches in the touch count (0 for start etc)
function [] = fix_tapCount_errors()

    file_init = 'C:\Users\pxb08145\OneDrive - University of Strathclyde\Documents\Research\Autism\Data\PlayCare\';
    
    name_save={};tC_save={};
    m = 0; n = m;
    %% cycle through folders
    D = dir(file_init); % A is a struct ... first elements are '.' and '..' used for navigation.
    g=0;
    for f = 1:length(D) % avoid using the first ones?
        currD = D(f).name; % Get the current subdirectory name
        file_loc =[file_init,currD,'\'];
        for i = 1:1
            skip=0;
            if i<10
                ntext = '0';
            else
                ntext = '';
            end
            sub_id = [currD,'.Sharing.TouchData.typed.csv'];
    
            if isfile([file_loc,sub_id])
                n=n+1;
                filename = [file_loc,sub_id];
                savename = ['currD',ntext,num2str(i)];
            else
                skip = 1;
            end
    
            if skip == 0
                tab=sortrows(readtable(filename),4);
                if iscell(tab.X(1)) %% Convert strings to doubles for X and Y
                    tab=sortrows(readtable(filename),14);
                    tab.X=strrep(tab.X,',','.');
                    tab.Y=strrep(tab.Y,',','.');
                    tab.X=str2double(tab.X);
                    tab.Y=str2double(tab.Y);
                end
                    
                
                tC = zeros(length(tab.X),1);
    
                starts = find(tab.TouchPhase==0);
                ends = find(tab.TouchPhase==3);
    
                if starts(1)~=1
                    starts = [1;starts];
                end
                
                [valid,invalid,starts,ends] = check_validity(tab,starts,ends);
                
                starts_tC = 1:1:length(starts);
                tC(starts)=starts_tC;
                pos=cell(1,max(starts_tC));
    
                viable = starts_tC;
                for j = 1 : length(valid)
                    vj = valid(j);
                    if tC(vj)==0
                        prevj=j;
                        current_Phase = tab.TouchPhase(vj);
                        nm_srch = 750;
                        dist_all=ones(nm_srch,1)*10000;
                        tC_check=[];
                        for jj = 1:nm_srch          % check up to nm_srch prior touches for best match
                            j=j-1;  % dodgy coding but goes back through the previous nm_srch touches
                            if j>0 && ~ismember(tC(valid(j)),tC_check) %%&& ~ismember(j,invalid)
                                vj = valid(j);
                                if size(pos{tC(vj)},1)>1
                                    xx=tab.X(vj);
                                    if length(pos{tC(vj)}(:,1))>4 %% spline susceptible to error use linear as well
                                        try
                                            yy = spline(pos{tC(vj)}(end-4:end,1),pos{tC(vj)}(end-4:end,2),xx);
                                        catch
                                            xx=2*pos{tC(vj)}(end,1)-pos{tC(vj)}(end-1,1);
                                            yy=2*pos{tC(vj)}(end,2)-pos{tC(vj)}(end-1,2);
                                        end
                                    else
                                        try
                                            yy = spline(pos{tC(vj)}(:,1),pos{tC(vj)}(:,2),xx);
                                        catch
                                            xx=2*pos{tC(vj)}(end,1)-pos{tC(vj)}(end-1,1);
                                            yy=2*pos{tC(vj)}(end,2)-pos{tC(vj)}(end-1,2);
                                        end
                                    end
                                    dist_all(jj)=sqrt((tab.X(valid(prevj))-xx)^2+(tab.Y(valid(prevj))-yy)^2);
                                else
                                    dist_all(jj)=sqrt((tab.X(valid(prevj))-tab.X(vj))^2+(tab.Y(valid(prevj))-tab.Y(vj))^2);
                                end
                                tC_check =[tC_check,tC(vj)];
                            end
                        end
                        [~,I]=min(dist_all); % closest match
                        while ~ismember(tC(valid(prevj-I)),viable) || tab.TouchPhase(valid(prevj-I))==3 % check if swipe has already ended
                            dist_all(I)=10000;
                            [~,I]=min(dist_all);
                            if min(dist_all) == 10000 %%|| prevj-I<=0
                                break
                            end
                        end
                        if min(dist_all) == 10000 % no connection found in previous nm_srch touches, create new start 
                            tC(valid(prevj))=max(tC)+1;
                            starts=[starts;valid(prevj)];
                            pos{tC(valid(prevj))}=[tab.X(valid(prevj)),tab.Y(valid(prevj))];
                            viable=[viable,max(viable)+1];
                            [num2str(valid(prevj)),' is now a start']
                        else
                            tC(valid(prevj))=tC(valid(prevj-I));
                            pos{tC(valid(prevj))}=[pos{tC(valid(prevj))};tab.X(valid(prevj)),tab.Y(valid(prevj))];
                        end
    
    
                        if tab.TouchPhase(valid(prevj))==3
                            viable(find(viable==tC(valid(prevj))))=[]; % remove swipe viability if it ends
                        end
                    else
                        pos{tC(vj)}=[pos{tC(vj)};tab.X(vj),tab.Y(vj)];
                    end
                end
                g=g+1;
                tC_save{g}=tC;
                nam_save{g}=currD;
            end
            
        end
    end
end