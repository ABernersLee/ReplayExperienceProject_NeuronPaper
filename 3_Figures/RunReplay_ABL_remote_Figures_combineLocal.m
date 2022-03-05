function RunReplay_ABL_remote_Figures_combineLocal(basedir)

%%%%% Remote vs local replays across 2 sessions (two different novel
%%%%% tracks)


%%
%%%%%
%%%%% dirs and params
%%%%%

load([basedir 'dirs_linear_remote'],'dirs')

cd(dirs.spikedatadir)
d2 = dir('*.mat');
testing = extractfield(d2,'name');
testing = erase(testing,'_1.mat'); testing = erase(testing,'_2.mat'); testing = erase(testing,'_3.mat');
d1 = unique(testing);
lapnum = 14; lapskip = 2;
jd_cutoff = .4; wc_cutoff = .6; coveragecutoff = 0; 
replays1 = []; replays2 = []; replays3 = [];
replabC = {'Local replays';'Remote replays'};
collab = {'Duration';'Slope';'Range over time';'mean(stepssize)'};


numcells = NaN(size(d1,2),2);
lapsall = [];
for idd = 1:size(d1,2)    
     lps = NaN(1,2);
    for itrack = 1:2
     load([d1{idd} '_' num2str(itrack)],'hp_cells','hpinterneurons','MidTime')
        if exist('hp_cells','var')
            numcells(idd,itrack) = sum(~ismember(hp_cells,hpinterneurons));
            clear hp_cells hpinterneurons
        end
        lps(itrack) = length(MidTime)-1;
    end
    lapsall = cat(1,lapsall,lps);
end

exl = find(sum(lapsall<14,2)>0);

%%

%%%%%
%%%%% Get out replays
%%%%%

replaysS = NaN(size(d1,2),lapnum,6,3);
for idd = 1:size(d1,2)
    if ismember(idd,exl)
        continue
    end
    for itrack = 1:2        
        load([d1{idd} '_' num2str(itrack)],'CandCorr','OutFR','CandDis','CandRS','CandSeq','MidTime','Spike','CandStepSize','replayparams')
        t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff;
        
        [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);
        I(I>length(MidTime)) = length(MidTime);
        t(I==0) = false;
        
                 %1/26/2022 change* need to change unites of slope!! 5*old#     *replayparams.Bin_Size/100*1000/5
        dat = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams.Bin_Size CandStepSize(t,4,2)*replayparams.Bin_Size idd*ones(size(I(t))) ];
        datS = NaN(lapnum,size(dat,2));
        for ilap = 1:lapskip:lapnum
            if ilap==lapnum && lapskip>1;continue;end
            if lapskip==1
                datS(ilap,:) = mean(dat(dat(:,1)==ilap,:),1);
            elseif lapskip==2                
                datS(ilap,:) = mean(dat(dat(:,1)==ilap | dat(:,1)==ilap+1,:),1);
            end
        end
        if itrack == 1
            replays1 = cat(1,replays1,dat);
            replaysS(idd,:,:,1) = datS;
        elseif itrack == 2
            replays2 = cat(1,replays2,dat);
            replaysS(idd,:,:,2) = datS;
        end
    end
    
    load([d1{idd} '_3'],'CandCorr','CandDis','CandRS','CandSeq')
     t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff;
     t(I==0) = false;           
     
                 %1/26/2022 change* need to change unites of slope!! 5*old#     *replayparams.Bin_Size/100*1000/5
    dat = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams.Bin_Size CandStepSize(t,4,2)*replayparams.Bin_Size idd*ones(size(I(t)))];
    replays3 = cat(1,replays3,dat);
    datS = NaN(lapnum,size(dat,2));
    for ilap = 1:lapskip:lapnum            
        if ilap==lapnum && lapskip>1;continue;end
        if lapskip==1
            datS(ilap,:) = mean(dat(dat(:,1)==ilap,:),1);
        elseif lapskip==2                
            datS(ilap,:) = mean(dat(dat(:,1)==ilap | dat(:,1)==ilap+1,:),1);
        end
    end
    replaysS(idd,:,:,3) = datS;
end

replaysC1 = cat(1,[replays1 ones(size(replays1,1),1)],[replays2 2*ones(size(replays2,1),1)]);
replaysC3 = [replays3 ones(size(replays3,1),1)];

if (size(replays1,1)>= size(replays2,1)) && (size(replays1,1)>=size(replays3,1))
    replays2 = cat(1,replays2,NaN(size(replays1,1)-size(replays2,1),size(replays2,2)));
    replays3 = cat(1,replays3,NaN(size(replays1,1)-size(replays3,1),size(replays3,2)));
elseif (size(replays2,1)>= size(replays1,1)) && (size(replays2,1)>=size(replays3,1))
    replays3 = cat(1,replays3,NaN(size(replays2,1)-size(replays3,1),size(replays3,2)));
    replays1 = cat(1,replays1,NaN(size(replays2,1)-size(replays1,1),size(replays1,2)));
elseif (size(replays3,1)>= size(replays2,1)) && (size(replays3,1)>=size(replays1,1))
    replays2 = cat(1,replays2,NaN(size(replays3,1)-size(replays2,1),size(replays2,2)));
    replays1 = cat(1,replays1,NaN(size(replays3,1)-size(replays1,1),size(replays1,2)));
end

if (size(replaysC1,1)>= size(replaysC3,1))
    replaysC3 = cat(1,replaysC3,NaN(size(replaysC1,1)-size(replaysC3,1),size(replaysC3,2)));
elseif (size(replaysC3,1)>= size(replaysC1,1)) 
    replaysC1 = cat(1,replaysC1,NaN(size(replaysC3,1)-size(replaysC1,1),size(replaysC1,2)));
end

replaysC = cat(3,replaysC1,replaysC3);


%% 

%%%%%
%%%%% plot across different replay aspects (col) and local vs remote
%%%%% (irep) replay type
%%%%%

colmult = [.1 .85 .85 .85];
for col = 2:3 %4
        
    yls = NaN(3,2);
    
    for irep = 1:2
        a(irep) = figure; hold on
        A = NaN(lapnum,2);
        for ilap = 1:lapskip:lapnum    
            if lapskip==1
                t = replaysC(:,1,irep)==ilap;
            elseif lapskip==2
                t = replaysC(:,1,irep)==ilap | replaysC(:,1,irep)==ilap+1;
            end
            A(ilap,1)=nanmean(replaysC(t,col,irep));
            A(ilap,2)=nanstd(replaysC(t,col,irep))/sqrt(sum(~isnan(replaysC(t,col,irep))));        
        end    
        A(sum(~isnan(A),2)==0,:) = [];        
        errorbar(A(:,1),A(:,2),'k','LineWidth',2);    
        yls(irep,:) = get(gca,'ylim');
        title(replabC{irep})
        xlabel('Lap')
        ylabel(collab{col-1})        
    end
    
    
    ny = [min(yls(:,1)) max(yls(:,2))]; 
    for irep = 1:2
        set(a(irep).Children,'ylim',ny) 
        axes(a(irep).Children)
        
        if col==2 %if duration, hypothesis is that correlation is greater than zero        
            [r,p] = corr(replaysC(replaysC(:,1,irep)<=lapnum,1,irep),replaysC(replaysC(:,1,irep)<=lapnum,col,irep),'rows','complete','tail','right');         
        elseif col>2 %if slope, hypothesis this that correlation is less than zero            
            [r,p] = corr(replaysC(replaysC(:,1,irep)<=lapnum,1,irep),replaysC(replaysC(:,1,irep)<=lapnum,col,irep),'rows','complete','tail','left');         
        end
        text(0,ny(1)+colmult(col-1)*range(ny),['1tail, N = ' num2str(sum(~isnan(replaysC(replaysC(:,1,irep)<=lapnum,1,irep)))) ', r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))],'FontSize',12)
        set(a(irep).Children,'xlim',[0 lapnum/lapskip+1])
        if lapskip>1            
            set(gca,'xtick',[1:2:lapnum/2],'xticklabel',2*[1:2:lapnum])
        end
        set(gcf,'Position',[2173         762         310         228])
        if col <4
            helper_savefig([basedir 'Figures\Final\2ways_' collab{col-1} '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_laps' num2str(lapnum) '_lapskip' num2str(lapskip) '_' replabC{irep}])
        end
        helper_saveandclosefig([basedir 'Figures\remote\2ways_' collab{col-1} '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_laps' num2str(lapnum) '_lapskip' num2str(lapskip) '_' replabC{irep}])
    end
    
        
    dat = [replaysC(replaysC(:,1,1)<=lapnum,col,1);replaysC(replaysC(:,1,2)<=lapnum,col,2)];
    lab1 = [replaysC(replaysC(:,1,1)<=lapnum,1,1);replaysC(replaysC(:,1,2)<=lapnum,1,2)];
    lab2 = [ones(size(replaysC(replaysC(:,1,1)<=lapnum,1,1)));2*ones(size(replaysC(replaysC(:,1,2)<=lapnum,1,2)))];
    [p,tbl] =anovan(dat,{lab1,lab2},'model','interaction','continuous',1,'display','off');

    figure; hold on;
    plab = {'Group';'Lap';'Interaction'};
    text(0.1,.9,[collab{col-1} '; ' replabC{1} ' and ' replabC{2} ':'])
    for ip = 1:3
        text(0,ip/4,[plab{ip} ': F = ' num2str(tbl{ip+1,6}) ' df = ' num2str(tbl{ip+1,3}) ',' num2str(sum(~isnan(dat))-tbl{ip+1,3}) ' , p = ' num2str(p(ip))])
    end
    set(gcf,'Position',[1919         133         861         403])
    axis off

    if col <4
        helper_savefig([basedir 'Figures\Final\2ways_' collab{col-1} '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_laps' num2str(lapnum) '_lapskip' num2str(lapskip) '_Interaction'])
    end
    helper_saveandclosefig([basedir 'Figures\remote\2ways_' collab{col-1} '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_laps' num2str(lapnum) '_lapskip' num2str(lapskip) '_Interaction'])


end
