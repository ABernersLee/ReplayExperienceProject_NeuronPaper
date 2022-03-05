function params = extract_nvt(dirs,params)
cd(params.daydir)

videostring ='VT1.nvt';
Rat_Name = params.Rat_Name;
Date = params.Date;
Track_Type = params.Track_Type;
% NovelDay = params.NovelDay;
RewardChangeDay = params.RewardChangeDay;

%extract raw nvt file
if ismac
    eval(sprintf('[ptimes,Xpos,Ypos,HD]=Nlx2MatVT_v3(''%s'',[1 1 1 1 0 0],0,1);',videostring));    
elseif ispc
    eval(sprintf('[ptimes,Xpos,Ypos,HD]=Nlx2MatVT(''%s'',[1 1 1 1 0 0],0,1);',videostring));
end
clear videostring

%convert times to seconds
ptimes = ptimes/1e6;

%concatonate into position data matrix
pos1=[ptimes',Xpos',Ypos',HD'];
clear ptimes Xpos Ypos HD

%if times are not in order, sort them
[~,index]=sortrows(pos1(:,1));
if sum(diff(index)<0)>0        
    pos1 = pos1(index,:);
    disp('There are position times out of order')    
end
clear index



%extract the run times
st = (repmat(pos1(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(pos1(:,1))]);
nd = (repmat(pos1(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(pos1(:,1))]);
rundat = (st&nd)';
clear st nd

%save each run seperately
for irun = 1:size(rundat,2)
    rawpos = pos1(rundat(:,irun),:);     

    if RewardChangeDay~=0 && irun>1 
        if (RewardChangeDay==-1 && irun==2) || (RewardChangeDay==1 && irun==3)
                RewardChangeRun = -1;
        elseif (RewardChangeDay==1 && irun==3) || (RewardChangeDay==1 && irun==2)
                RewardChangeRun = 1;
        end
    else RewardChangeRun = 0;
    end

    NovelRun = params.Novel;
    Run = params.Run;
    params.RewardChangeRun = RewardChangeRun;
    params.ident = [params.Rat_Name '_' num2str(params.Date) '_Run' num2str(Run)];
    
    if exist([dirs.spikedatadir '\' params.ident '.mat'],'file')>0
        save([dirs.spikedatadir '\' params.ident '.mat'],'rawpos','Rat_Name','Date','Run','Track_Type','NovelRun','RewardChangeRun','params','-append')
    else
        save([dirs.spikedatadir '\' params.ident '.mat'],'rawpos','Rat_Name','Date','Run','Track_Type','NovelRun','RewardChangeRun','params')
    end
    clear rawpos
end