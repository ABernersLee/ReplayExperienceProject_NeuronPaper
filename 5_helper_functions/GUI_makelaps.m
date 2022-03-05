function [laps_coverspace,laps_twoarms,laps_singlepass,headingarm,error_correct] = GUI_makelaps(behavior, behave_change_log, pos,armpos, linposcat,params)

%%%%%
%%%%% GUI that makes laps a few different ways (see below, mostly for
%%%%% multi-arm enviornments, MidPoint is used for liner ones) and then
%%%%% makes sure it all looks good
%%%%%

%first, make laps_coverspace, for linear laps this is just back and forth,
    %for multi-arm laps this is going to all 3 arms
%then, make laps_twoarms, for linear laps this is just back and forth,
    %for multi-arm laps this is going to 2 arms
for laptype = 1:2
    laplog = NaN(length(behave_change_log),1);
    clear laps
    if params.Track_Type == 1
        lap_starts = find(behave_change_log(:,1));
        otherside = find(behave_change_log(:,end));            
        otherside(otherside<lap_starts(1,1),:) = [];
        for o = 1:length(otherside)
            lpend = lap_starts(find(lap_starts>otherside(o,1),1,'first'))-1;
            if o == length(otherside) && isempty(lpend)
                continue
            else
                lap_ends(o,1) = lpend;
            end
        end
        lap_ends = unique(lap_ends);
        laps = sort([lap_starts(1:end-1) lap_ends]);
        clear lap_ends lap_starts otherside o lap_ends
    else
        [~,thisarm] = find(~isnan(behavior(find(~isnan(nansum(behavior,2)),1,'first'),:)));
        i = 1;
        laps(i,1) = find(~isnan(behavior(:,thisarm)),1,'first');
        clear thisarm
        stop = false;
        while max(laps)<length(behavior) && ~stop
            if laptype == 1
                ending = laps(i,:)+find(all((cumsum(behavior(laps(i,:):end,:)==4)>0),2),1,'first'); % has to have gone to each end at least once
            elseif laptype == 2
                    ending = laps(i,:)+find(sum((cumsum(behavior(laps(i,:):end,:)==4)>0),2)>=2,1,'first'); % has to have gone to at least two
            end
            i = i+1;
            if isempty(find(behavior(ending:end,:)==1,1,'first')) || ending>=length(behavior)
                stop = true;
            else                
            testind = find(diff(armpos(ending:end))~=0,1,'first'); %then has to move arm
            laps(i,1) = testind+ending;
            end
        end    
        clear i stop ending testind
        laps(:,2) = [laps(2:end)-1;length(behavior)];

        %if he didn't finish the last lap, cut it
        if sum(all((cumsum(behavior(laps(end,1):laps(end,2),:)==4)>0),2))==0
            laps(end,:) = [];
        end
    end

    for ilap = 1:size(laps,1)
        laplog(laps(ilap,1):laps(ilap,2),1) = ilap;        
    end    
    clear ilap laps


    vc = varycolor(max(laplog)); 
    vc = vc(randperm(size(vc,1)),:); 
    figure; hold on;
    for ilap = 1:max(laplog)
         plot(pos(laplog==ilap,1),linposcat(laplog==ilap),'.','MarkerEdgeColor',vc(ilap,:))
    end
    axis tight


    %checks if this is okay
    prompt = 'Is this okay? y/n';
    title(prompt)
    yn = input([prompt ':  '],'s');
    if strcmp(yn,'n')
        return
    else
        close gcf
        clear vc ilap prompt yn
        if laptype == 1
            laps_coverspace = laplog; clear laplog
        elseif laptype == 2
            laps_twoarms = laplog; clear laplog
        end
    end
end

%finally, make laps_singlepass, for linear laps this is just back and forth,
    %for multiple arms its in on one arm, out to the next arm

if params.Track_Type == 2
    ilap = 1;
    beh = nansum(behavior,2);
    laps = NaN(size(behavior,1),1);
    naht = true(size(laps));
    st = find(beh==4,1,'first');
    starm = armpos(st);
    nd = find(beh==4 & armpos~=starm,1,'first');
    laps(st:nd-1) = ilap; ilap = ilap+1;
    naht(st:nd-1) = false;
    while ~isempty(nd) && ~isempty(st)
       st = find(beh==4 & naht,1,'first');
       starm = armpos(st);
       nd = find(beh==4 & armpos~=starm & naht,1,'first');
       laps(st:nd-1) = ilap;
       ilap = ilap+1;
       naht(st:nd-1) = false;
    end              
    clear ilap beh naht ns st 

    vc = varycolor(nanmax(laps)); 
    vc = vc(randperm(size(vc,1)),:); 
    figure; hold on;
    for ilap = 1:max(laps)
         plot(pos(laps==ilap,1),linposcat(laps==ilap),'.','MarkerEdgeColor',vc(ilap,:))
    end
    axis tight


    %checks if this is okay
    prompt = 'Is this okay? y/n';
    title(prompt)
    yn = input([prompt ':  '],'s');
    if strcmp(yn,'n')
        return
    else
        close gcf        
        laps_singlepass = laps;
        clear vc ilap prompt yn laps
    end


    headingarm = NaN(size(laps_singlepass,1),1);   
    for ilap = 1:max(laps_singlepass)
        headingarm(laps_singlepass==ilap) = armpos(find(laps_singlepass==ilap,1,'last'));
    end
    clear ilap 
else
    headingarm = NaN(size(linposcat,1),1);
    laps_singlepass = headingarm;
end

error_correct = NaN(size(laps_singlepass,1),1);

close all