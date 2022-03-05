function [behavior,behave_change_log, behave_ind] = GUI_make_behavior(linposcat,armpos,linposcatnan,params)

%%%% makes behavior epoch marticies (middle, stem, neck, lick area)

okay = 0;
while okay == 0
    %make Track_Places
    if params.Track_Type==1
            Track_Places = NaN(4,1);
            [h,c] = hist(linposcat-min(linposcat),100);
            figure; hold on
            plot(h,'k-')
            title('Click on where the bottom or left lick area ends')   
            [Track_Places(1),~] = getpts;  
            title('Click on where the bottom or left neck area ends')     
            [Track_Places(2),~] = getpts;  
            title('Click on where the top or right neck area starts')     
            [Track_Places(3),~] = getpts;  
            title('Click on where the top or right lick area starts')     
            [Track_Places(4),~] = getpts;
            hh = hist(Track_Places,1:length(c));
            Track_Places = c(hh==1)';
    elseif params.Track_Type == 2
         Track_Places = NaN(3,size(params.arms,1));
         for iarm = 1:size(params.arms,1)        
            [h,c] = hist(linposcat(armpos==iarm)-min(linposcat(armpos==iarm)),100);
            figure; hold on
            plot(h,'k-')
            Track_Places1 = NaN(3,1);
            title(['Arm ' num2str(iarm) ', Click on where the middle area ends'])   
            [Track_Places1(1),~] = getpts;  
            title(['Arm ' num2str(iarm) ', Click on where the neck area starts'])     
            [Track_Places1(2),~] = getpts;  
            title(['Arm ' num2str(iarm) ', Click on where the lick area starts'])     
            [Track_Places1(3),~] = getpts;  
            hh = hist(Track_Places1,1:length(c));
            Track_Places(:,iarm) = c(hh==1)';
        end

    end

    cutoff = bsxfun(@plus,Track_Places,((min(linposcatnan(:,2:end)))));
    close gcf

    % get 
    behavior = ones(size(linposcatnan,1),size(Track_Places,2));
    for j = 1:size(Track_Places,2) %arm
        for i = 1:size(Track_Places,1) %bevahior
            if params.Track_Type == 2
                index = nansum(linposcatnan(:,2:end),2)>cutoff(i,j) & armpos==j;
            else
                index = nansum(linposcatnan(:,2:end),2)>cutoff(i,j);
            end
            behavior(index,j) = i+1;        
        end
        behavior(armpos~=j,j) = NaN;
    end

    c = varycolor(size(Track_Places,1)+1);
    figure; hold on
    for j = 1:size(Track_Places,1)+1   %behavior
        for i = 1:size(Track_Places,2) %arm
            pos2 = linposcatnan(:,1+i);
            plot(linposcatnan(behavior(:,i)==j,1),pos2(behavior(:,i)==j),'.','MarkerEdgeColor',c(j,:))
        end    
    end
    axis tight
    prompt = 'Is this okay? y/n';
    title(prompt)
    yn = input([prompt ':  '],'s');
    if strcmp(yn,'n')
        okay = 0;
    elseif strcmp(yn,'y')
        okay = 1;
    end
end


%this makes epochs
behave_time = cell(max(max(behavior))*2-2,1);
if size(params.Run_Times,1)>1
    ind = linposcatnan(:,1)>=params.Run_Times(params.Run,1) & linposcatnan(:,1)<=params.Run_Times(params.Run,2);
else
    ind = linposcatnan(:,1)>=params.Run_Times(1,1) & linposcatnan(:,1)<=params.Run_Times(1,2);
end
p = linposcatnan(ind,:);
b = nansum(behavior(ind,:),2);
% positive change is going out, negative change is going back in
chng=[diff(b);0];
laps = false(length(p),max(max(behavior))*2-2);
ilaps = NaN(length(p),1);    
for i = 1:max(max(behavior))-1
    laps(:,i) = (chng>0 & b==i) | laps(:,i);
    ilaps(chng>0 & b==i,1) = i;
end
for i = max(max(behavior)):max(max(behavior))*2-2
    laps(:,i) = (chng<0 & b==(max(max(behavior))*2-i)) | laps(:,i);
    ilaps(chng<0 & b==(max(max(behavior))*2-i),1) = i;        
end
clear i chng


ilaps2 = NaN(length(ilaps),1);
%     ilapsch = ilaps2;
ind = find(~isnan(ilaps),1,'first'); 
s = ilaps(ind);
ilaps2(1:ind) = s;
while ind<=length(ilaps)
    s2 = s+1;
    if s2>max(max(behavior))*2-2
        s2 = 1;
    end
    ind1 = find(ilaps(ind+1:end)==s2,1,'first');
    notfulllap = false;
    %this makes it so that from when they enter to the last time they
    %exit the lick zone, that counts as lick time. This works well.
    if params.Track_Type ==2 && s == 5
        ind_back = find(b(1:ind)==4,1,'last');
        ind_back2 = find(ilaps2(1:ind)==3,1,'last');
        if ~isempty(ind_back) && ~isempty(ind_back2)                
            ilaps2(ind_back2:ind_back) = 3;
        end
    end            
    if params.Track_Type ==2 && s == 6
        ind_back = find(b(1:ind)==3,1,'last');
        ind_back2 = find(ilaps2(1:ind)==4,1,'last');
        if ~isempty(ind_back) && ~isempty(ind_back2)                
            ilaps2(ind_back2:ind_back) = 4;
        end
    end

    if params.Track_Type == 1 
       if s==7 
           ind_back = find(b(1:ind)==4,1,'last');
           ind_back2 = find(ilaps2(1:ind)==5,1,'last');
           if ~isempty(ind_back) && ~isempty(ind_back2)                
                ilaps2(ind_back2:ind_back) = NaN; %was 5
           end
           ind_otherend = find(b(ind:end)==5 | b(ind:end)==4,1,'first');
           ind_end = find(b(ind:end)==1,1,'first');
           if ind_end>ind_otherend
               notfulllap= true;
           end
       elseif s == 3
           ind_back = find(b(1:ind)==2,1,'last');
           ind_back2 = find(ilaps2(1:ind)==1,1,'last');
           if ~isempty(ind_back) && ~isempty(ind_back2)                
                ilaps2(ind_back2:ind_back) = NaN; %was 1
           end

           ind_otherend = find(b(ind:end)==1 | b(ind:end)==2,1,'first');
           ind_end = find(b(ind:end)==5,1,'first');
            if ind_end>ind_otherend
               notfulllap= true;
           end
       elseif s== 5
           ind_back = find(b(1:ind)==4,1,'last');
           ind_back2 = find(ilaps2(1:ind)==3,1,'last');
           if ~isempty(ind_back) && ~isempty(ind_back2)                
                ilaps2(ind_back2:ind_back) = NaN; %was 3 
           end
       elseif s== 1
           ind_back = find(b(1:ind)==2,1,'last');
           ind_back2 = find(ilaps2(1:ind)==7,1,'last');
           if ~isempty(ind_back) && ~isempty(ind_back2)                
                ilaps2(ind_back2:ind_back) = NaN; %was 7 
           end               
       end
    end
    if ~notfulllap
        ilaps2(ind+1:ind1+ind) = s;
    elseif notfulllap
        ilaps2(ind+1:ind1+ind) = NaN;
    end
%         ilapsch(ind1+ind) = s;
    s = s2; ind = ind1+ind;
end

if params.Track_Type == 2
    %here I post-hoc make sure that segments that were assigned worngly are righted
    ilaps2_new = ilaps2;
    ilaps2_new(b==1,1) = 6;        

    problem = [0; b==2 & ilaps2_new~=1 & ilaps2_new~=5];
    if sum(problem)>0
        diffprob = diff(problem);
        st_prob = find(diffprob==1);
        nd_prob = find(diffprob==-1);
        if length(nd_prob)<length(st_prob)
            nd_prob = [nd_prob;size(diffprob,1)];
        end
        stnd = [st_prob nd_prob];
        stnd(diff(stnd')<50,:) = [];
        pdiff = nansum(p(stnd(:,2),2:end),2)-nansum(p(stnd(:,1),2:end),2);
        for iproblem = 1:length(pdiff)
            if pdiff(iproblem,1)>30
                ilaps2_new(stnd(iproblem,1):stnd(iproblem,2),1) = 1;
            elseif pdiff(iproblem,1)<-30
                ilaps2_new(stnd(iproblem,1):stnd(iproblem,2),1) = 5;
            else
                dat = nansum(p(stnd(:,1):stnd(:,2),2:end),2);
                pos_ind = find(max(dat)==dat,1,'first');
                ilaps2_new(stnd(:,1):stnd(:,1)+pos_ind,1) = 1;
                ilaps2_new(stnd(:,1)+pos_ind+1:stnd(:,2),1) = 5;
            end
        end
        clear problem diffprob st_prob nd_prob stnd pdiff iproblem
    else clear problem
    end

    problem = [0; b==3 & ilaps2_new~=2 & ilaps2_new~=4];
    if sum(problem)>0
        diffprob = diff(problem);
        st_prob = find(diffprob==1);
        nd_prob = find(diffprob==-1);
        if length(nd_prob)<length(st_prob)
            nd_prob = [nd_prob;size(diffprob,1)];
        end
        stnd = [st_prob nd_prob];
        stnd(diff(stnd')<50,:) = [];
        pdiff = nansum(p(stnd(:,2),2:end),2)-nansum(p(stnd(:,1),2:end),2);
        for iproblem = 1:length(pdiff)
            if pdiff(iproblem,1)>30
                ilaps2_new(stnd(iproblem,1):stnd(iproblem,2),1) = 2;
            elseif pdiff(iproblem,1)<-30
                ilaps2_new(stnd(iproblem,1):stnd(iproblem,2),1) = 4;
            else
                dat = nansum(p(stnd(:,1):stnd(:,2),2:end),2);
                pos_ind = find(max(dat)==dat,1,'first');
                ilaps2_new(stnd(:,1):stnd(:,1)+pos_ind,1) = 2;
                ilaps2_new(stnd(:,1)+pos_ind+1:stnd(:,2),1) = 4;
            end
        end
        clear problem diffprob st_prob nd_prob stnd pdiff iproblem
    else clear problem
    end
    ilaps2_new(ilaps2==3) = 3; %in case any top bits got erroded by the #4
    ilaps2 = ilaps2_new; 
end

%to change/make ilapsch --> matters for laps
ilapsch = NaN(length(ilaps2),1);
ilapsch([0;diff(ilaps2)]~=0) = ilaps2([0;diff(ilaps2)]~=0,1);


laps2 = false(length(p),max(max(behavior))*2-2);
laps2ch = laps2;
for i = 1:size(laps2,2)
    laps2(ilaps2==i,i) = i;
    laps2ch(ilapsch==i,i) = i;
    behave_time{i,1} = [behave_time{i,1};[p(ilaps2==i) i*ones(sum(ilaps2==i),1)]];
end  

behave_change_log = laps2ch;
behave_change_ind = ilapsch;
behave_log = laps2;
behave_ind = ilaps2;

plotind = 1:size(p,1);
c = varycolor(size(behave_log,2));
figure; hold on
for j = 1:size(behave_log,2)   %behavior
    plot(plotind(behave_ind==j),p(behave_ind==j,2),'.','MarkerEdgeColor',c(j,:))
end

axis tight
prompt = 'Is this okay? y/n';
title(prompt)
yn = input([prompt ':  '],'s');
if strcmp(yn,'n')
   return
elseif strcmp(yn,'y')
    close gcf
end

close all