function plot_ripple_events_over_laps_replays(basedir)

%%%%%
%%%%%
%%%%% This looks at changes of replays over laps by averaging each session
%%%%% and also does some ANOVAs to determine changes over laps 
%%%%%
%%%%%

load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')
cd(dirs.spikedatadir)
dat.Lap = 14;
dat.Lapskip = 1;
toplot = true;
if ~isfolder([dirs.figdir '\ReplayOverLaps\'])
    mkdir([dirs.figdir '\ReplayOverLaps\'])
end

%%%%% plotting replay aspects (col) across types of replay (itype)
for col = 2:6
    if col == 2; label = 'Replay Duration';
    elseif col == 3; label = 'Replay Slope';
    elseif col == 4; label = 'Replay Slope (range)';
    elseif col == 5; label = 'Number of Steps';
    elseif col == 6; label = 'Proportion of the track';
    end
    d2 = dir('*.mat');
    wc_cutoff = .6;
    jd_cutoff = .4;
    coveragecutoff = 0; %.4;
    addlab = '';
    AllDatAll = [];
    
    
    titlab = {'Novel Tracks'; 'Familiar Tracks';...
        'Novel Tracks - Reverse';'Novel Tracks - Forward';...
        'RewChange Tracks - Reverse';'RewChange - Forward'; 'RewChange'};    
    for itype = 1:7
        if itype==1
            novel=1; fwdrev = NaN; rpass = true; rewardchange = NaN;
        elseif itype==2
            novel=0; fwdrev = NaN; rpass = true; rewardchange = NaN; 
            %changed to reward change==0 and rpass is false, 4/17/2020
        elseif itype==3
            novel=1; fwdrev = 1; rpass = true; rewardchange = NaN;

        elseif itype==4
            novel=1; fwdrev = 2; rpass = true; rewardchange = NaN;
        elseif itype==5
            novel=0; rewardchange = [-1 1]; fwdrev = 1; rpass = false;
        elseif itype==6
            novel=0; rewardchange = [-1 1]; fwdrev = 2; rpass = false;
        elseif itype ==7
            novel = 0; rewardchange = [-1 1]; rpass = false; fwdrev = NaN;
        end
        


        AllDat = [];
        for id = 1:size(d2,1)

            load(d2(id).name,'CandSeq','OutFR','CandDis','CandCorr','MidTime','params','pos','CandRS','CandStepSize','replayparams')

            if (novel==params.Novel) && (rpass || ismember(params.RewardChangeRun,rewardchange)) % || novel==2
                maxmove = CandDis/size(OutFR,2);
                CandPassCrit = abs(CandCorr)>wc_cutoff & maxmove<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; %new
                if fwdrev==1
    %                 %CandRS(:,4), positive numbers here is out, which is =
    %                 %dirdat==1 (not 0), which is going up the track, so
    %                 %positive slopes (CandRS(:,2) or CandCorr)
                    reversereplay = (CandCorr<0 & CandRS(:,4)>0) | (CandCorr>0 & CandRS(:,4)<0);
                    CandPassCrit = CandPassCrit & reversereplay;
                elseif fwdrev==2
                    fwdreplay = (CandCorr>0 & CandRS(:,4)>0) | (CandCorr>0 & CandRS(:,4)>0);
                    CandPassCrit = CandPassCrit & fwdreplay;
                end

                CandSeq = CandSeq(CandPassCrit,:);

                [~,IR]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(pos(:,1))]);
                if sum(IR==0)>0
                    IR(IR==0) = NaN;
                end
                IR(IR>length(MidTime)) = length(MidTime);
                
                
                 %1/26/2022 change* need to change unites of slope!! 5*old#     *replayparams.Bin_Size/100*1000/5
                newdat = [IR diff(CandSeq(:,1:2),[],2) ...
                    abs(CandRS(CandPassCrit,2))*replayparams.Bin_Size/100*1000/5 ...
                    CandStepSize(CandPassCrit,3,2)*replayparams.Bin_Size...
                    CandStepSize(CandPassCrit,2,2) ...
                    CandStepSize(CandPassCrit,5,2) id*ones(size(IR))]; 
                AllDat = cat(1,AllDat,newdat); 
                
            end
            
            
        end
        AllDatAll = cat(1,AllDatAll,[AllDat(AllDat(:,1)<=dat.Lap,col) ...
            ones(size(AllDat(AllDat(:,1)<=dat.Lap,:),1),1)*itype AllDat(AllDat(:,1)<=dat.Lap,1)]); % each replay


        if toplot 
            figure; hold on

            A=NaN(dat.Lap,2);
            for j=1:dat.Lapskip:dat.Lap    
                t=find(AllDat(:,1)==j | AllDat(:,1)==j+dat.Lapskip-1);
                A(j,1)=nanmean(AllDat(t,col));
                A(j,2)=nanstd(AllDat(t,col))/sqrt(sum(~isnan(AllDat(t,col))));        
            end    
            A(sum(~isnan(A),2)==0,:) = [];


            errorbar(A(:,1),A(:,2),'k','LineWidth',2);


            xlabel('Lap number')
            set(gca,'XTick',[0:5:dat.Lap])
            t=find(AllDat(:,1)>0 & AllDat(:,1)<=dat.Lap);
            
            if itype>1 && col<6               
                if col==2 || col==5
                    [r,p]=corr(AllDat(t,1),AllDat(t,col),'rows','complete','tail','right');
                else                    
                    [r,p]=corr(AllDat(t,1),AllDat(t,col),'rows','complete','tail','left');
                end
                rtlab = 'OneTail';                
            elseif itype==1 || col==6                
                [r,p]=corr(AllDat(t,1),AllDat(t,col),'rows','complete');
                rtlab = 'TwoTail';
            end

            ylabel(label)    
            xlim([0 (dat.Lap/dat.Lapskip)+1])
            if p<.05 
                title({['\color{black}' titlab{itype} '\color{red}, ' rtlab ' r=' num2str(round(r,2,'significant')) ...
                         ', p=' num2str(round(p,2,'significant'))]},'FontSize',12); 
            else
                title({['\color{black}' titlab{itype} '\color{gray}, ' rtlab ' r=' num2str(round(r,2,'significant')) ...
                         ', p=' num2str(round(p,2,'significant'))]},'FontSize',12);
            end
            helper_saveandclosefig([dirs.figdir '\ReplayOverLaps\Replays ' label ' ' titlab{itype} ' Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab])
            
             A=NaN(dat.Lap,2);
             AllDatSes = [];
            for j=1:dat.Lapskip:dat.Lap    
                t=find(AllDat(:,1)==j | AllDat(:,1)==j+dat.Lapskip-1);                
                for ises = 1:max(AllDat(:,end))
                    if sum(ises==AllDat(:,end))>0
                        AllDatSes = cat(1,AllDatSes,[j nanmean(AllDat(intersect(t,find(ises==AllDat(:,end))),2:end),1)]);
                    end
                end                
                A(j,1)=nanmean(AllDatSes(AllDatSes(:,1)==j,col));
                A(j,2)=nanstd(AllDatSes(AllDatSes(:,1)==j,col))/sqrt(sum(~isnan(AllDatSes(AllDatSes(:,1)==j,col))));        
            end    
            
            
            figure; hold on
            A(sum(~isnan(A),2)==0,:) = [];
            
            errorbar(A(:,1),A(:,2),'k','LineWidth',2);
            xlabel('Lap number')
            set(gca,'XTick',[0:5:dat.Lap])
            t=find(AllDatSes(:,1)>0 & AllDatSes(:,1)<=dat.Lap);
            [r,p]=corr(AllDatSes(t,1),AllDatSes(t,col),'rows','complete');  
            ylabel(label)    
            xlim([0 (dat.Lap/dat.Lapskip)+1])            
            if p<.05 
                title({['\color{black}' titlab{itype} ' Sessions, N =  ' num2str(size(AllDatSes,1)) ',\color{red} r=' num2str(round(r,2,'significant')) ...
                         ', p=' num2str(round(p,2,'significant'))]},'FontSize',12); 
            else
                title({['\color{black}' titlab{itype} ' Sessions, N =  ' num2str(size(AllDatSes,1)) ',\color{gray} r=' num2str(round(r,2,'significant')) ...
                         ', p=' num2str(round(p,2,'significant'))]},'FontSize',12);
            end
            if itype==1
                helper_savefig([dirs.figdir '\Final\Replays ' label ' ' titlab{itype} ' Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab '_Sessions'])
            end
            
            helper_saveandclosefig([dirs.figdir '\ReplayOverLaps\Replays ' label ' ' titlab{itype} ' Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab '_Sessions'])
            
            
        end
    end
    
    if itype==7
        y = NaN(4,3);
        AllDatAll2 = AllDatAll(AllDatAll(:,2)==1 | AllDatAll(:,2)==2,:);
        [y(1,:),tbl] = anovan(AllDatAll2(:,1),{AllDatAll2(:,2),AllDatAll2(:,3)},'model','interaction','continuous',2,'display','off');
        AllDatAll2 = AllDatAll(AllDatAll(:,2)==3 | AllDatAll(:,2)==4,:);
        y(2,:) = anovan(AllDatAll2(:,1),{AllDatAll2(:,2),AllDatAll2(:,3)},'model','interaction','continuous',2,'display','off');
        AllDatAll2 = AllDatAll(AllDatAll(:,2)==5 | AllDatAll(:,2)==6,:);
        y(3,:) = anovan(AllDatAll2(:,1),{AllDatAll2(:,2),AllDatAll2(:,3)},'model','interaction','continuous',2,'display','off');
        AllDatAll2 = AllDatAll(AllDatAll(:,2)==1 | AllDatAll(:,2)==7,:);
        [y(4,:),tbl2] = anovan(AllDatAll2(:,1),{AllDatAll2(:,2),AllDatAll2(:,3)},'model','interaction','continuous',2,'display','off');

        figure; hold on;
        plab = {'Group';'Lap';'Interaction'};
        iy = 1;
        text(0,.9,[titlab{(iy-1)*2+1} ' vs ' titlab{(iy-1)*2+2}])
        for ip = 1:3
            text(0,ip/4,[plab{ip} ': F = ' num2str(tbl{ip+1,6}) ' df = ' num2str(tbl{ip+1,3}) ',' num2str(sum(AllDatAll(:,2)==1 | AllDatAll(:,2)==2)-tbl{ip+1,3}) ' , p = ' num2str(y(iy,ip))])
        end
        set(gcf,'Position',[1919         133         861         403])
        axis off
        if col<4
            helper_savefig([dirs.figdir '\Final\Replays ' label ' AllGroupsStats_ANOVA_Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab '_NovelVFam'])
        end
        helper_saveandclosefig([dirs.figdir '\ReplayOverLaps\Replays ' label ' AllGroupsStats_ANOVA_Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab '_NovelVFam'])



        figure; hold on;
        for iy = 1:3
            text(0,iy/4,[titlab{(iy-1)*2+1} ' vs ' titlab{(iy-1)*2+2} ': group: ' num2str(y(iy,1)) ', laps: ' num2str(y(iy,2)) ', interaction: ' num2str(y(iy,3))])
        end
        iy = 4;
        text(0,iy/4,[titlab{1} ' vs ' titlab{7} ': group: ' num2str(y(iy,1)) ', laps: ' num2str(y(iy,2)) ', interaction: ' num2str(y(iy,3))])
        set(gcf,'Position',[1919         133         861         403])
        axis off
        helper_saveandclosefig([dirs.figdir '\ReplayOverLaps\Replays ' label ' AllGroupsStats Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab])


        figure; hold on;
        plab = {'Group';'Lap';'Interaction'};
        iy = 4;
        text(0,.9,[titlab{1} ' vs ' titlab{7}])
        for ip = 1:3
            text(0,ip/4,[plab{ip} ': F = ' num2str(tbl2{ip+1,6}) ' df = ' num2str(tbl2{ip+1,3}) ',' num2str(sum(AllDatAll(:,2)==1 | AllDatAll(:,2)==7)-tbl2{ip+1,3}) ' , p = ' num2str(y(iy,ip))])
        end
        set(gcf,'Position',[1919         133         861         403])
        axis off
        helper_saveandclosefig([dirs.figdir '\ReplayOverLaps\Replays ' label ' AllGroupsStats_ANOVA_Laps' num2str(dat.Lap) ' Lapskips' num2str(dat.Lapskip) addlab '_NovelVRew'])
    end

%Fernandez-Ruiz, which has opposite effect:
% To detect ripples, the wide-band signal was band-pass filtered 
% (difference-of-Gaussians; zero-lag, linear phase FIR), and instantaneous 
% power was computed by clipping at 4 SD, rectified and low-pass filtered (20). 
% The low-pass filter cut-off was at a frequency corresponding to p cycles of the mean bandpass 
% (for 80-250 Hz band-pass, the low-pass was 55 Hz). The mean and SD of baseline LFP were computed 
% based on the power during non-REM sleep. Subsequently, the power of the 
% non-clipped signal was computed, and all events exceeding 4 SD from the mean were detected. 
% Events were then expanded  
% until the (non-clipped) power fell below 1 SD; short events (< 15 ms) were discarded. 
    

% I use 150-250Hz and 2SD, back to 1SD anf exclude any <50ms. with the
% baseline being whole session, no sleep.
end
