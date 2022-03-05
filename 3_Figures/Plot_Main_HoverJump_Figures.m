function Plot_Main_HoverJump_Figures(basedir,coveragecutoff)
%%%%% These figures plot quantification of aspects of replay
%%%%% changing over laps

%%


%%%%%
%%%%%
%%%%% set up directories and parameters
%%%%%
%%%%%

addlab = '_pearson_'; 
toexclude = false;
wc_cutoff = .6;
jd_cutoff = .4;
figlab = ['LapByLap_' addlab];
titlab = ' - Lap Fields';
load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')
cd(dirs.spikedatadir)

figfolder = [dirs.figdir 'LapByLap\'];
finalfigfolder = [dirs.figdir 'Final\'];

if ~isfolder(figfolder)
    mkdir(figfolder)
end
if ~isfolder(finalfigfolder)
    mkdir(finalfigfolder)
end

List=dir('*.mat');
numses = size(List,1);

if ~toexclude    
    excl_ind = true(numses,1);
else
    figlab = ['Excl_' figlab];
    titlab = [titlab ' - Sessions Excluded'];
end

% get basic numbers from each session to look at and use later
a1=0; a2 = 0;
numcells = NaN(size(List,1),8);
sessionnonSig = []; sessionSig = [];
for ListNo= 1:numses 
    if excl_ind(ListNo)==0
        continue
    end
    load(List(ListNo).name,'hp_cells','hpinterneurons')
    if exist('hp_cells','var')
        numcells(ListNo,1) = sum(~ismember(hp_cells,hpinterneurons));
        clear hp_cells hpinterneurons
    end
    load(List(ListNo).name,'CandCorr','CandDis','OutFR','CandStepSize','params','MidTime','-mat')
    numcells(ListNo,2)=sum(abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff);
    numcells(ListNo,3) = params.Novel;
    numcells(ListNo,4) = params.Track_Type;
    if contains(List(ListNo).name,'Chesapeake')
        numcells(ListNo,5) = false;
    else
        numcells(ListNo,5) = true;
    end
    numcells(ListNo,6) = params.Run;
    numcells(ListNo,8) = params.RewardChangeRun;
    numcells(ListNo,7) = length(MidTime)-1;
    
    CandPassCrit=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; %new
    sessionSig = [sessionSig; ones(sum(CandPassCrit),1)*ListNo];
    sessionnonSig = [sessionnonSig; ones(sum(~CandPassCrit),1)*ListNo];
    a1=a1+sum(CandPassCrit);
    a2=a2+sum(~CandPassCrit);    
    
end
excl_ind = numcells(:,1)>=20;



%%

%%%%%
%%%%%
%%%%% get data out for replays
%%%%%
%%%%%


% Lap No, Duration, RS Slope, Corr, Max Move Distance, Step number, Step Distance (range), Step Distance (sum),
% Proportion of Track Covered, Step Duration,
%  'Kendall R';'Spearman R';'Linear Fit';'Cubic Fit';'Cubic/Linear';
% gamma frequency, gamma power, gamma phase, gamma MRV, 
% ListNo, Run, Track Type, Novelty, RewardChange

AllSig=NaN(a1,25); 
AllnonSig=NaN(a2,25);
RipHzSig = NaN(a1,2);
RipHznonSig = NaN(a2,2);
hoverA=[];% lap, duration, ListNo, track type, Novelty, RewardChange
moveA=[];% lap, duration, dis, ListNo, track type, Novelty, RewardChange
hoverSmA=[];% lap, duration, ListNo, track type, Novelty, RewardChange
moveSmA=[];% lap, duration, dis, ListNo, track type, Novelty, RewardChange
hoverfwdA=[];% lap, duration, ListNo, track type, Novelty, RewardChange
movefwdA=[];% lap, duration, dis, ListNo, track type, Novelty, RewardChange
GammaS = []; %Degree, dis, lap, ListNo, track type, Novelty, RewardChange
count1=0; count2=0;

for ListNo=1:numses 
    if excl_ind(ListNo)==0
        continue
    end
    load(List(ListNo).name,'Radjusted','KendallR','SpearmanR','OutFRlaps','InFRlaps',...
        'hover*','move*','StepGamma','replayparams','GammaProp','GammaPropHigh','SeqGamma',...
        'CandStepSize','Spike','PosShuffleCorr','PosShuffleDis','CandCorr','CandDis','OutFR','CandSeq','CandRS','MidTime','params',...
        'hover_sm*','move_sm*','hover_fwd*','move_fwd*','hovertime2','movetime2','rp_freq','-mat')

    stepdur = NaN(size(CandSeq,1),1);
    for ii = 1:size(CandSeq,1)
          stepdur(ii,1) = mean(diff([movetime2(movetime2(:,3)==ii,1)]))*5;
    end
    
    maxmove = CandDis/size(OutFR,2);
    pc = ((sum(abs(PosShuffleCorr)>=abs(CandCorr),2)+1)./(size(PosShuffleCorr,2)+1))<.05;
    pd = ((sum(PosShuffleDis<=CandDis,2)+1)./(size(PosShuffleDis,2)+1))<.05;
    CandPassCrit = abs(CandCorr)>wc_cutoff & maxmove<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; % & pc & pd; %new

    [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);
    if sum(I==0)>0
        I(I==0) = NaN;
    end
    I(I>length(MidTime)) = length(MidTime);
    
    t = CandPassCrit;                                                %1/26/2022 change* need to change unites of slope!! 5*old#                  %step distance, 3 is range/time, 4 is mean(abs(moves)), 5 is the proportion of the track it covers
    AllSig(count1+1:count1+sum(t),:)=[I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 abs(CandCorr(t)) maxmove(t) CandStepSize(t,2,2) CandStepSize(t,3,2)*replayparams.Bin_Size ...
        CandStepSize(t,4,2)*replayparams.Bin_Size CandStepSize(t,5,2) stepdur(t) ... 
        abs(KendallR(t)) abs(SpearmanR(t)) Radjusted(t,1) Radjusted(t,3) Radjusted(t,3)-Radjusted(t,1) ... %     'Kendall R';'Spearman R';'Linear Fit';'Cubic Fit';'Cubic/Linear';
        GammaProp(t,1) zscore(GammaProp(t,2)) SeqGamma(t,:) params.RewardChangeDay*ones(sum(t),1) ListNo*ones(sum(t),1) ...
        params.Run*ones(sum(t),1) params.Track_Type*ones(sum(t),1) params.Novel*ones(sum(t),1) params.RewardChangeRun*ones(sum(t),1)];
    RipHzSig(count1+1:count1+sum(t),:) = [I(t) rp_freq(t)];
    count1=count1+sum(t); 
    
    t = ~CandPassCrit;                                                  %1/26/2022 change* need to change unites of slope!! 5*old# 
    AllnonSig(count2+1:count2+sum(t),:)=[I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 abs(CandCorr(t)) maxmove(t) CandStepSize(t,2,2) CandStepSize(t,3,2)*replayparams.Bin_Size ...
        CandStepSize(t,4,2)*replayparams.Bin_Size CandStepSize(t,5,2) stepdur(t) ...
        abs(KendallR(t)) abs(SpearmanR(t)) Radjusted(t,1) Radjusted(t,3) Radjusted(t,3)-Radjusted(t,1)... %     'Kendall R';'Spearman R';'Linear Fit';'Cubic Fit';'Cubic/Linear';
        GammaProp(t,1) zscore(GammaProp(t,2)) SeqGamma(t,:) params.RewardChangeDay*ones(sum(t),1) ListNo*ones(sum(t),1) ...
        params.Run*ones(sum(t),1) params.Track_Type*ones(sum(t),1) params.Novel*ones(sum(t),1) params.RewardChangeRun*ones(sum(t),1)];
    RipHznonSig(count2+1:count2+sum(t),:) = [I(t) rp_freq(t)];
    count2=count2+sum(t); 
          
    %gamma phase, step distance (max, mean), event
    GammaS = [GammaS;[StepGamma(:,1) StepGamma(:,3) I(StepGamma(:,4)) CandPassCrit(StepGamma(:,4)) ListNo*ones(size(StepGamma,1),1) params.Run*ones(size(StepGamma,1),1) params.Track_Type*ones(size(StepGamma,1),1) params.Novel*ones(size(StepGamma,1),1) params.RewardChangeRun*ones(size(StepGamma,1),1)]];    
    
    hover2(hover2(:,5)~=1,:) = []; %use only hovers in the middle of the replay (exclude start and end)
    
    hoverA=[hoverA;[I(hover2(:,4)) CandPassCrit(hover2(:,4)) hover2(:,1)*5 ListNo*ones(size(hover2,1),1) params.Run*ones(size(hover2,1),1) params.Track_Type*ones(size(hover2,1),1) params.Novel*ones(size(hover2,1),1) params.RewardChangeRun*ones(size(hover2,1),1)]];    
    moveA=[moveA; [I(move2(:,4)) CandPassCrit(move2(:,4)) move2(:,1)*5 move2(:,3)*replayparams.Bin_Size ListNo*ones(size(move2,1),1) params.Run*ones(size(move2,1),1) params.Track_Type*ones(size(move2,1),1) params.Novel*ones(size(move2,1),1) params.RewardChangeRun*ones(size(move2,1),1)]];   

    hoverSmA=[hoverSmA;[I(hover_sm2(:,4)) CandPassCrit(hover_sm2(:,4)) hover_sm2(:,1)*5 ListNo*ones(size(hover_sm2,1),1) params.Run*ones(size(hover_sm2,1),1) params.Track_Type*ones(size(hover_sm2,1),1) params.Novel*ones(size(hover_sm2,1),1) params.RewardChangeRun*ones(size(hover_sm2,1),1)]];    
    moveSmA=[moveSmA; [I(move_sm2(:,4)) CandPassCrit(move_sm2(:,4)) move_sm2(:,1)*5 move_sm2(:,3)*replayparams.Bin_Size ListNo*ones(size(move_sm2,1),1) params.Run*ones(size(move_sm2,1),1) params.Track_Type*ones(size(move_sm2,1),1) params.Novel*ones(size(move_sm2,1),1) params.RewardChangeRun*ones(size(move_sm2,1),1)]];   
    
    hoverfwdA=[hoverfwdA;[I(hover_fwd2(:,4)) CandPassCrit(hover_fwd2(:,4)) hover_fwd2(:,1)*5 ListNo*ones(size(hover_fwd2,1),1) params.Run*ones(size(hover_fwd2,1),1) params.Track_Type*ones(size(hover_fwd2,1),1) params.Novel*ones(size(hover_fwd2,1),1) params.RewardChangeRun*ones(size(hover_fwd2,1),1)]];    
    movefwdA=[movefwdA; [I(move_fwd2(:,4)) CandPassCrit(move_fwd2(:,4)) move_fwd2(:,1)*5 move_fwd2(:,3)*replayparams.Bin_Size ListNo*ones(size(move_fwd2,1),1) params.Run*ones(size(move_fwd2,1),1) params.Track_Type*ones(size(move_fwd2,1),1) params.Novel*ones(size(move_fwd2,1),1) params.RewardChangeRun*ones(size(move_fwd2,1),1)]];   
    clear hover1 move1 move_sm1 hover_sm1 move_fwd1 hover_fwd1
    clear hover2 move2 move_sm2 hover_sm2 move_fwd2 hover_fwd2
end
clear hover move
hover = hoverA; move = moveA;
hoverSm = hoverSmA; moveSm = moveSmA;
hoverfwd = hoverfwdA; movefwd = movefwdA;

 %zscore gamma power
for ListNo=1:max([max(AllSig(:,20)) max(AllnonSig(:,20))])
    if excl_ind(ListNo)==0
        continue
    end
    CandPassCrit=find(AllSig(:,20)==ListNo);
    AllSig(CandPassCrit,17)=nanzscore(AllSig(CandPassCrit,17),[],1); 
    
    CandPassCrit=find(AllnonSig(:,20)==ListNo);
    AllnonSig(CandPassCrit,17)=nanzscore(AllnonSig(CandPassCrit,17),[],1); 
end




%%

%%%%%
%%%%%
%%%%% get data out for the cells
%%%%%
%%%%%

clear FRA FRB FR FRI FRAI FRBI FRAH FRBH FRAM FRBM
sz = NaN(2,32,2); FRind = []; FRindH = []; FRindM = []; FRindI = []; 
FRlaps = []; FRlapsINT = []; FRindINT = []; FRlapsProp = [];
for iSig = 1:2
    clear FR1 FR2 FR3 FR4 FR5 FR6 FR7 FR8 FR9 FR10 ...
        FR11 FR12 FR13 FR14 FR15 FR16 FR17 FR18 FR19 ...
        FR20 FR21 FR22 FR23 FR24 FR25 FR26 FR27 FR28 FR29 FR30 FR31 FR32
    for itype = 0:32
        
        for ListNo=1:numses 
            if excl_ind(ListNo)==0
                continue
            end
            load(List(ListNo).name,'params','-mat');
            
            if itype == 0
                if iSig == 1
                load(List(ListNo).name,'OutFRlaps','InFRlaps','OutFRlapsINT','InFRlapsINT')
                FRlaps1 = squeeze(max(cat(2,OutFRlaps,InFRlaps),[],2)); %max
                FRlapsINT1 = squeeze(mean(cat(2,OutFRlapsINT,InFRlapsINT),2)); %mean
                if size(FRlapsINT1,2)~=size(OutFRlapsINT,3)
                    FRlapsINT1 = FRlapsINT1';
                end
                if ~isempty(FRlaps)                    
                   if size(FRlaps,2)>size(FRlaps1,2)
                       FRlaps1 = cat(2,FRlaps1,NaN(size(FRlaps1,1),size(FRlaps,2)-size(FRlaps1,2)));
                   elseif size(FRlaps,2)<size(FRlaps1,2)                       
                       FRlaps = cat(2,FRlaps,NaN(size(FRlaps,1),size(FRlaps1,2)-size(FRlaps,2)));
                   end                       
                end
                FRlaps  = cat(1,FRlaps,FRlaps1);
                
                 if ~isempty(FRlapsINT)
                   if size(FRlapsINT,2)>size(FRlapsINT1,2)
                       FRlapsINT1 = cat(2,FRlapsINT1,NaN(size(FRlapsINT1,1),size(FRlapsINT,2)-size(FRlapsINT1,2)));
                   elseif size(FRlapsINT,2)<size(FRlapsINT1,2)
                       FRlapsINT = cat(2,FRlapsINT,NaN(size(FRlapsINT,1),size(FRlapsINT1,2)-size(FRlapsINT,2)));
                   end                       
                 end
                if ~isempty(FRlapsINT1)
                    FRlapsINT  = cat(1,FRlapsINT,FRlapsINT1);
                    FRindINT = [FRindINT; ListNo*ones(size(FRlapsINT1,1),1) ...
                        params.Run*ones(size(FRlapsINT1,1),1) params.Track_Type*ones(size(FRlapsINT1,1),1) ...
                        params.Novel*ones(size(FRlapsINT1,1),1) params.RewardChangeRun*ones(size(FRlapsINT1,1),1)];
                end
                end
            else
                if itype == 1
                    load(List(ListNo).name,'SpikePerEvent')
                    data=SpikePerEvent;
                    clear SpikePerEvent
                elseif itype == 2
                    load(List(ListNo).name,'ActRatio')
                     data=ActRatio;
                     clear ActRatio
                elseif itype == 3
                    load(List(ListNo).name,'SpikePerAct')
                     data=SpikePerAct;
                     clear SpikePerAct
                elseif itype ==4
                    load(List(ListNo).name,'FRPerCell')
                       data=FRPerCell;
                       clear FRPerCell      
                elseif itype ==5
                    load(List(ListNo).name,'ActRatioRate')
                     data=ActRatioRate;
                     clear ActRatioRate
                elseif itype == 6
                    load(List(ListNo).name,'SpikePerActRate')
                     data=SpikePerActRate;
                     clear SpikePerActRate                
                elseif itype == 7
                    load(List(ListNo).name,'SpikePerEventHover')
                    data=SpikePerEventHover;
                    clear SpikePerEventHover
                elseif itype == 8
                    load(List(ListNo).name,'ActRatioHover')
                     data=ActRatioHover;
                     clear ActRatioHover                                  
                elseif itype == 9
                    load(List(ListNo).name,'SpikePerActHover')
                     data=SpikePerActHover;
                     clear SpikePerActHover
                elseif itype == 10
                    load(List(ListNo).name,'FRPerCellHover')
                       data=FRPerCellHover;
                       clear FRPerCellHover
                elseif itype == 11
                    load(List(ListNo).name,'ActRatioRateHover')
                     data=ActRatioRateHover;
                     clear ActRatioRateHover   
                elseif itype == 12
                    load(List(ListNo).name,'SpikePerActRateHover')
                     data=SpikePerActRateHover;
                     clear SpikePerActRateHover
                elseif itype == 13
                    load(List(ListNo).name,'SpikePerEventMove')
                    data=SpikePerEventMove;
                    clear SpikePerEventMove
                elseif itype == 14
                    load(List(ListNo).name,'ActRatioMove')
                     data=ActRatioMove;
                     clear ActRatioMove
                elseif itype == 15
                    load(List(ListNo).name,'SpikePerActMove')
                     data=SpikePerActMove;
                     clear SpikePerActMove
                elseif itype == 16
                    load(List(ListNo).name,'FRPerCellMove')
                       data=FRPerCellMove;
                       clear FRPerCellMove
                elseif itype == 17
                    load(List(ListNo).name,'ActRatioRateMove')
                     data=ActRatioRateMove;
                     clear ActRatioRateMove
                elseif itype == 18                
                    load(List(ListNo).name,'SpikePerActRateMove')
                     data=SpikePerActRateMove;
                     clear SpikePerActRateMove
                elseif itype ==19
                    load(List(ListNo).name,'FRPerCellFirst80')
                    data=FRPerCellFirst80(:,:,:,1);
                    clear FRPerCellFirst80
                elseif itype ==20
                    load(List(ListNo).name,'FRPerCellFirst80')
                    data=FRPerCellFirst80(:,:,:,2);
                    clear FRPerCellFirst80
                elseif itype ==21
                    load(List(ListNo).name,'ActRatioRateFirst80')
                    data=ActRatioRateFirst80(:,:,:,1);
                    clear ActRatioRateFirst80
                elseif itype ==22
                    load(List(ListNo).name,'ActRatioRateFirst80')
                     data=ActRatioRateFirst80(:,:,:,2);
                     clear ActRatioRateFirst80
                elseif itype == 23
                    load(List(ListNo).name,'SpikePerActRateFirst80')
                     data=SpikePerActRateFirst80(:,:,:,1);
                     clear SpikePerActRateFirst80
                elseif itype == 24
                    load(List(ListNo).name,'SpikePerActRateFirst80')
                     data=SpikePerActRateFirst80(:,:,:,2);
                     clear SpikePerActRateFirst80
                elseif itype == 25
                    load(List(ListNo).name,'SpikePerEventIN')
                    data=SpikePerEventIN;
                    clear SpikePerEventIN
                elseif itype == 26
                    load(List(ListNo).name,'ActRatioIN')
                     data=ActRatioIN;
                     clear ActRatioIN
                elseif itype == 27
                    load(List(ListNo).name,'SpikePerActIN')
                     data=SpikePerActIN;
                     clear SpikePerActIN
                elseif itype ==28
                    load(List(ListNo).name,'FRPerCellIN')
                       data=FRPerCellIN;
                       clear FRPerCellIN      
                elseif itype ==29
                    load(List(ListNo).name,'ActRatioRateIN')
                     data=ActRatioRateIN;
                     clear ActRatioRateIN
                elseif itype == 30
                    load(List(ListNo).name,'SpikePerActRateIN')
                     data=SpikePerActRateIN;
                     clear SpikePerActRateIN]
                 elseif itype ==31
                    load(List(ListNo).name,'FRall') 
                    data=FRall;
                     clear FRall
                 elseif itype == 32
                    load(List(ListNo).name,'FRallIN') 
                     data=FRallIN;
                     clear FRallIN
                end

                if itype == 1 && iSig==1
                   FRind = [FRind; ListNo*ones(size(data,1),1) params.Run*ones(size(data,1),1) params.Track_Type*ones(size(data,1),1) params.Novel*ones(size(data,1),1) params.RewardChangeRun*ones(size(data,1),1)];
                end

                if itype == 7 && iSig==1
                   FRindH = [FRindH; ListNo*ones(size(data,1),1) params.Run*ones(size(data,1),1) params.Track_Type*ones(size(data,1),1) params.Novel*ones(size(data,1),1) params.RewardChangeRun*ones(size(data,1),1)];
                end

                if itype == 13 && iSig==1
                   FRindM = [FRindM; ListNo*ones(size(data,1),1) params.Run*ones(size(data,1),1) params.Track_Type*ones(size(data,1),1) params.Novel*ones(size(data,1),1) params.RewardChangeRun*ones(size(data,1),1)];
                end

                if itype == 26 && iSig==1
                   FRindI = [FRindI; ListNo*ones(size(data,1),1) params.Run*ones(size(data,1),1) params.Track_Type*ones(size(data,1),1) params.Novel*ones(size(data,1),1) params.RewardChangeRun*ones(size(data,1),1)];
                end

                if ~exist((['FR' num2str(itype)]),'var') 
                    eval(['FR' num2str(itype) '=data(:,:,(iSig-1)*3+1);'])            
                    eval(['sz(:,itype,iSig)=size(FR' num2str(itype) ');'])
                else            
                    if isnan(sz(:,itype,iSig))
                        eval(['sz(:,itype,iSig)=size(FR' num2str(itype) ');'])
                    end
                    if sz(2,itype,iSig)>size(data,2)
                        a=sz(2,itype,iSig)-size(data,2);
                        eval(['FR' num2str(itype) '=[FR' num2str(itype) ';[data(:,:,(iSig-1)*3+1) NaN*ones(size(data,1),a)]];'])
                        eval(['sz(:,itype,iSig)=size(FR' num2str(itype) ');'])
                    elseif sz(2,itype,iSig)<size(data,2)
                        a=size(data,2)-sz(2,itype,iSig);
                        eval(['FR' num2str(itype) '=[FR' num2str(itype) ' NaN*ones(sz(1,itype,iSig),a)];'])
                        eval(['FR' num2str(itype) '=[FR' num2str(itype) ';data(:,:,(iSig-1)*3+1)];'])
                        eval(['sz(:,itype,iSig)=size(FR' num2str(itype) ');'])
                    else
                        eval(['FR' num2str(itype) '=[FR' num2str(itype) ';data(:,:,(iSig-1)*3+1)];'])
                        eval(['sz(:,itype,iSig)=size(FR' num2str(itype) ');'])
                    end
                end
            end
        end
    end
    if iSig == 1
        FRA = cat(3,FR1,FR2,FR3,FR4,FR5,FR6,FR19,FR20,FR21,FR22,FR23,FR24,FR31);
        FRAI = cat(3,FR25,FR26,FR27,FR28,FR29,FR30,FR32);
        FRAH = cat(3,FR7,FR8,FR9,FR10,FR11,FR12);
        FRAM = cat(3,FR13,FR14,FR15,FR16,FR17,FR18);
        clear FR1 FR2 FR3 FR4 FR5 FR6 FR7 FR8 FR9 FR10 FR11 FR12 ...
            FR13 FR14 FR15 FR16 FR17 FR18 FR19 FR20 FR21 FR22 FR23 ...
            FR24 FR25 FR26 FR27 FR28 FR29 FR30 FR31 FR32
    elseif iSig == 2
        FRB = cat(3,FR1,FR2,FR3,FR4,FR5,FR6,FR19,FR20,FR21,FR22,FR23,FR24,FR31);
        FRBI = cat(3,FR25,FR26,FR27,FR28,FR29,FR30,FR32);
        FRBH = cat(3,FR7,FR8,FR9,FR10,FR11,FR12);
        FRBM = cat(3,FR13,FR14,FR15,FR16,FR17,FR18);
        clear FR1 FR2 FR3 FR4 FR5 FR6 FR7 FR8 FR9 FR10 FR11 FR12 FR13 ...
            FR14 FR15 FR16 FR17 FR18 FR19 FR20 FR21 FR22 FR23 FR24 FR25...
            FR26 FR27 FR28 FR29 FR30 FR31 FR32
    end
end
FR = cat(4,FRA,FRB);
FRI = cat(4,FRAI,FRBI);
FRH = cat(4,FRAH,FRBH);
FRM = cat(4,FRAM,FRBM);

clearvars -except finalfigfolder figfolder sessionSig sessionnonSig ...
    numses RipHzSig RipHznonSig FRindINT FRlapsINT gammaS hoverSm ...
    moveSm hoverfwd movefwd addlab coveragecutoff toexclude ids excl_ind...
    FRindI FRI List titlab AllSig AllnonSig FR FRH FRM FRindH FRindM ...
    hover move FRind numcells wc_cutoff jd_cutoff ids ...
    downsamplab figlab FRlaps

%%

%%%%%
%%%%%
%%%%% things for plotting like yaxes and labels of the aspects to be
%%%%% plotted
%%%%%
%%%%%

dat.AllSig = AllSig;
dat.AllnonSig = AllnonSig;
dat.nlab= {'Familiar';'Novel';'All Novel and Familiar'};
dat.rlab = {'No Change in Reward';'Change In Reward';...
    'Increase in Reward';'Decrease in Reward';''};
dat.tlab = {'Linear Track';'Multi-Arm Maze';'Open Field'};
col_label = {'Duration, second';'Slope, meters per sec';'Correlation';...
    'Max Jump Dist';'Number of Steps';...
    'Step Distance (est by range dt), cm';...
    'Step Distance mean(abs(move)), cm';'Proportion of Track Covered';...
    'Step Duration, ms';'Kendall R';'Spearman R';'Linear Fit';...
    'Cubic Fit';'Cubic - Linear';
    'Gamma Frequency, Hz';'z-scored gamma power';...
    'Gamma Phase of Spikes';'Gamma MRV';...
    'Hover Duration, ms';'Move Duration, ms';'Move Distance, cm';...
    'Spike number per event';'Proportion Cells Recruited';...
    'Spike number when Active';'Firing rate Per cell, Hz'; ...
    'Prop. cells recruited per sec';'FR Per Cell When Active';...
    'Spike number per Hover event';...
    '% Cells Recruited to Hover';'Spike number when Active in Hover';...
    'Firing rate Per cell in Hover, Hz';...
    '% Cells Recruited to Hover Per Second';...
    'FR Per Cell When Active in Hover'; ...
    'Spike number per Jump event';'% Cells Recruited to Jump';...
    'Spike number when Active in Jump';...
    'Firing rate Per cell in Jump, Hz'; ...
    '% Cells Recruited to Jump Per Second';...
    'FR Per Cell When Active in Jump'...
    ;'Firing rate first 80ms, Hz';'Firing rate after 80ms, Hz';...
    '% Cells Recruited Per Sec first 80ms';...
    '% Cells Recruited Per Sec after 80ms'...
    ;'FR Per Cell When Active first 80ms';...
    'FR Per Cell When Active after 80ms';
    'IN Spike number per event';'Proportion IN Recruited';...
    'Spike number IN when Active';'Firing rate Per IN, Hz'; ...
    'Prop IN Recruited Per Second';'FR Per IN When Active'; ...
    'Hover Duration (sm)';'Move Duration (sm)';'Move Distance (sm)';...
    'Hover Duration(fwd)';'Move Duration (fwd)';'Move Distance (fwd)'};
siglab = {'nonSig&Sig';'SigOnly'};
yl1 = [0.127331663439379,0.187165712948113;1.88636691649093,8;0.243654196820355,1;0.166193038165101,0.698593465965816;4.43597102854015,7.14047020523272;2.76172985746776,5.41484137929250;22.8826170297274,30.1114128421843;35.4369369305897,38.8109671397487;-0.504054377976601,0.373343774774513;-200,200;0.0692591028796730,0.151608324916799;8,13.3770216045783;11.5525851444245,16.6184843010216;20,67.1882246241458;0.197336585767740,1.20000000000000;0.0777402627726523,0.450444681637543;2,3.49482038571412;1.49629967658803,8;0.112942815743844,0.860456220739988;3.10727235917155,15.2776505254286;0.0280627722930466,0.121772121083554;0.0191109115321186,0.0862170415193409;1.19069476095638,1.84197304189955;1.58196173322910,8.14101021519188;0,0.781442798581782;2.76398971692514,14.5497650013382;0.0231137410592409,0.171194055125773;0.0187078949229797,0.120000000000000;1.17031615166199,1.51893829697684;1.43349498150120,8.23475299733518;0,1;2.64137784554730,14.4850394577371];
yl2 = [0.127331663439379,0.187165712948113;1.88636691649093,8;0.243654196820355,1;0.166193038165101,0.698593465965816;4.43597102854015,7.14047020523272;2.76172985746776,5.41484137929250;22.8826170297274,30.1114128421843;35.4369369305897,38.8109671397487;-0.504054377976601,0.373343774774513;-200,200;0.0692591028796730,0.151608324916799;8,13.3770216045783;11.5525851444245,16.6184843010216;20,67.1882246241458;0.197336585767740,1.20000000000000;0.0777402627726523,0.450444681637543;2,3.49482038571412;1.49629967658803,8;0.112942815743844,0.860456220739988;3.10727235917155,15.2776505254286;0.0280627722930466,0.121772121083554;0.0191109115321186,0.0862170415193409;1.19069476095638,1.84197304189955;1.58196173322910,8.14101021519188;0,0.781442798581782;2.76398971692514,14.5497650013382;0.0231137410592409,0.171194055125773;0.0187078949229797,0.120000000000000;1.17031615166199,1.51893829697684;1.43349498150120,8.23475299733518;0,1;2.64137784554730,14.4850394577371];
yl3 = [0.124648004036648,0.186134089413620;1.88060444083018,7.19085689607295;0.268158290478303,1;0.185290509894633,0.741180055873105;4,8;2.52043163541446,5;20,40.6026079284530;35,40.7030212999717;-0.680015300736440,0.436165934139866;-200,200;0.0659690831151955,0.153325743040444;8.93599197497508,15;9.74455310966830,16.0521485029523;20,83.3495669135557;0.172061531441926,1.11388615966216;0.0637950713608366,0.431481071655911;2,3.35696099275019;1.33657747819794,8;0.104132093607276,1;3.71432882883999,20;0.0293534640578060,0.150000000000000;0.0187122278342719,0.0903758949896064;1.20000000000000,2;1.39069608057638,8;0,0.810663762488145;3.69769721100795,20;0,0.154074790933140;0.0153409646188725,0.106323118782410;1.16784753354907,1.44079246509011;1.28980773123192,8.09641228027566;0,1;3.41003829579974,22.4914234181444];
yl4 = [0.120000000000000,0.190000000000000;1.50000000000000,9;0.700000000000000,0.880000000000000;0.140000000000000,0.320000000000000;4,7.50000000000000;1.50000000000000,5;20,42;34,43;-1,0.800000000000000;-300,150;0.0500000000000000,0.180000000000000;8,20;9.50000000000000,14.5000000000000;18,40;0.200000000000000,0.910000000000000;0.0600000000000000,0.360000000000000;2.18000000000000,3.60000000000000;1.50000000000000,6.80000000000000;0.100000000000000,0.750000000000000;4.50000000000000,22;0.0200000000000000,0.160000000000000;0.0200000000000000,0.0800000000000000;1.29000000000000,1.75000000000000;1,7.20000000000000;0.100000000000000,0.900000000000000;4,16;0.0200000000000000,0.140000000000000;0.0100000000000000,0.0900000000000000;1.20000000000000,1.45000000000000;1.50000000000000,6.90000000000000;0.100000000000000,1;3.50000000000000,18];
yl5 = [0.124648004036648,0.200000000000000;0,7.39894549808760;0.251828859484305,1;0,0.762476094114175;3.73717111874132,8;1.97484734930291,5;20,50.4171181866496;33.0591164297214,42.2549123354542;-1,0.486519453125691;-200,200;0.0597651615760427,0.200000000000000;8.93599197497508,22.8248890991915;8.35162826961141,25.9641941385921;0,100;0.155658795907176,1.11388615966216;0.0587078899132009,0.431481071655911;1.85304718546532,4;0,8;0.0789761081509952,1;3.71432882883999,25.0768442805042;0,0.158591167465330;0.0161541659197696,0.0903758949896064;1.20000000000000,2.13731014309849;1.17184252320011,8;0,0.810663762488145;3.69769721100795,21.1407250243807;0,0.154074790933140;0,0.106323118782410;1,1.59703784337777;0,8.09641228027566;0,1;0,32.6238926117088];
yl = [min([yl1(:,1) yl2(:,1) yl3(:,1) yl4(:,1) yl5(:,1)],[],2) max([yl1(:,2) yl2(:,2) yl3(:,2) yl4(:,2) yl5(:,2)],[],2)];
yl = [yl; zeros(6,2)];
yl = zeros(size(col_label,1),2);




%%

%%%%%
%%%%%
%%%%% Plot all aspects across laps - 1 big figure but also small subplots.
%%%%% This also produces additional figures that aren't used in the paper.
%%%%%
%%%%%

% plots each type of session, with both significant and nonsignificant
% replays on them

close all
plotsmaller =true;
dat.sz = [6,10];
ilap = 14;
dat.Lap = ilap;
dat.Lapskip = 1;
ylims = NaN(2,size(col_label,1),12);
%numcol = size(col_label,1);
in = 1;
TrackType = 1;

for sigonly = 0:1 
    for isNovel = 0:2 %0:2      
        for isDValue = [1 4] %1:4 

            if (sigonly==0 && isNovel~=1) || (sigonly==0 && isDValue~=4) || (isDValue==4 && isNovel==2) || (isDValue<4 && isNovel~=2)
                continue
            end

            if coveragecutoff>0 && ~(sigonly==1 && isDValue==4 && isNovel==1)
                continue
            end

            if sigonly==1 && isDValue==4 && isNovel==1
                toplotsmall = [1:9 15:16 18:27];
            elseif (sigonly==0 && isDValue==4 && isNovel==1) || (sigonly==1 && isDValue==4 && isNovel==0) || (sigonly==1 && isDValue==1 && isNovel==2)
                toplotsmall = [1:2];
            else
                toplotsmall = [];
            end

            f = figure; hold on
            hps = NaN(size(col_label,1)+1,2,2);
            for col = 2:size(col_label,1)+1
                if col<18
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFAlldat_SN(f,col,dat,col_label{col-1},yl(col-1,:),[],isNovel,isDValue,TrackType,0,sigonly,finalfigfolder);
                elseif col==18
                     [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFAlldat_phase(f,col,dat,col_label{col-1},yl(col-1,:),[],isNovel,isDValue,TrackType,0,sigonly);
                elseif col==19
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFAlldat_SN(f,col,dat,col_label{col-1},yl(col-1,:),[],isNovel,isDValue,TrackType,0,sigonly,finalfigfolder);
                elseif col==20
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFhoverjump_SN(f,col-18,dat,col_label{col-1},yl(col-1,:),col,hover,isNovel,isDValue,TrackType,0,sigonly);
                elseif col>20 && col<23
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFhoverjump_SN(f,col-19,dat,col_label{col-1},yl(col-1,:),col,move,isNovel,isDValue,TrackType,0,sigonly);
                elseif col>=23 && col<29
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TF_FR_SN(f,col-22,dat,col_label{col-1},yl(col-1,:),col,FR,isNovel,isDValue,FRind,TrackType,0,false,sigonly);
                elseif col>=29 && col<35
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TF_FR_SN(f,col-28,dat,col_label{col-1},yl(col-1,:),col,FRH,isNovel,isDValue,FRindH,TrackType,0,false,sigonly);
                elseif col>=35 && col<41
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TF_FR_SN(f,col-34,dat,col_label{col-1},yl(col-1,:),col,FRM,isNovel,isDValue,FRindM,TrackType,0,false,sigonly);
                elseif col>=41 && col<47
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TF_FR_SN(f,col-40,dat,col_label{col-1},yl(col-1,:),col,FR,isNovel,isDValue,FRind,TrackType,0,false,sigonly);
                elseif col>=47 && col<53
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TF_FR_SN(f,col-46,dat,col_label{col-1},yl(col-1,:),col,FRI,isNovel,isDValue,FRindI,TrackType,0,false,sigonly);
                elseif col==53
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFhoverjump_SN(f,col-51,dat,col_label{col-1},yl(col-1,:),col,hoverSm,isNovel,isDValue,TrackType,0,sigonly);
                elseif col>53 && col<56
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFhoverjump_SN(f,col-52,dat,col_label{col-1},yl(col-1,:),col,moveSm,isNovel,isDValue,TrackType,0,sigonly);
                elseif col==56
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFhoverjump_SN(f,col-54,dat,col_label{col-1},yl(col-1,:),col,hoverfwd,isNovel,isDValue,TrackType,0,sigonly);
                elseif col>56 && col<59
                    [hps(col,1,:),hps(col,2,:)] = plotsubplot_TFhoverjump_SN(f,col-55,dat,col_label{col-1},yl(col-1,:),col,movefwd,isNovel,isDValue,TrackType,0,sigonly);
                end
                ylims(:,col-1,in) = get(gca,'ylim');
            end
            set(f,'Position',[ 30        -495        2893        1486])

            in = in+1;
            if isDValue<4
                suptitle([dat.nlab{isNovel+1} ' ' dat.tlab{TrackType} ' with ' dat.rlab{isDValue+1} titlab])                    
                thissavefig = [figlab siglab{sigonly+1} ' ' num2str(dat.Lap) '-' num2str(dat.Lapskip) '-laps Linear Decoding WC' num2str(wc_cutoff) ' JD' num2str(jd_cutoff) 'coveragecutoff ' num2str(coveragecutoff) ' ' [dat.nlab{isNovel+1} ' ' dat.tlab{TrackType} ' with ' dat.rlab{isDValue+1}]];
                helper_savefig([figfolder thissavefig])
            else
                suptitle([dat.nlab{isNovel+1} ' ' dat.tlab{TrackType} titlab])
                thissavefig = [figlab siglab{sigonly+1} ' ' num2str(dat.Lap) '-' num2str(dat.Lapskip) '-laps Linear Decoding WC' num2str(wc_cutoff) ' JD' num2str(jd_cutoff) 'coveragecutoff ' num2str(coveragecutoff) ' ' [dat.nlab{isNovel+1} ' ' dat.tlab{TrackType}]];
                helper_savefig([figfolder thissavefig])
            end
            if plotsmaller
                pp = size(f.Children):-1:1;
             for ip = 1:size(f.Children)        
                if ~ismember(ip,toplotsmall)
                    continue
                end
                if strcmp(f.Children(pp(ip)).HandleVisibility,'on')
                    g(pp(ip)) = figure; hold on;                            
                    ps = g(pp(ip)).Children(end).Position;
                    axis off
                    II = copyobj(f.Children(pp(ip)),g(pp(ip)));
                    set(g(pp(ip)).Children(end-1),'Position',ps)
                    set(g(pp(ip)).Children(end-1),'FontSize',26)
                    set(gcf,'Position',[ 2063         229         752         551])
                    set(g(pp(ip)).Children(end-1),'OuterPosition',[0.0530    0.0118    0.8016    0.9854])
                    set(g(pp(ip)).Children(end-1),'PlotBoxAspectRatio',[ 1.0000    0.7538    0.7538])
                    g(pp(ip)).Children(end-1).Children(1).LineWidth = 3;
                    if size(g(pp(ip)).Children(end-1).Children,1)==2
                        g(pp(ip)).Children(end-1).Children(2).LineWidth = 3;
                    end
                    g(pp(ip)).Children(end-1).TitleFontWeight = 'normal';
                    helper_saveandclosefig([finalfigfolder thissavefig '_smplot_' g(pp(ip)).Children(end-1).YLabel.String])
                end                            
             end
            end

            close all

        end
    end
end


