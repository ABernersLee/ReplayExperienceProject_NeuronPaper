%%%% get difference between remote and local replays

% Using novel track data I orginally got (later adding a few more sessions
% I found and clustered for RemoteSleep)

clearvars
basedir = 'E:\test\';
load([basedir 'dirs_linear_remote'],'dirs')

%% Adding orginal novel sessions
if 1 % to get data out and save
cd(dirs.spikedatadirOrig)
d2 = dir('*.mat');
cd(dirs.remotedatadir)
d1 = dir('*e*');
eids = [4 5; 5 6; 13 14; 24 25; 32 33; 35 36; 35 37; 36 37; 42 43; 51 52; 53 54; 53 55; 54 55; 57 58; 57 59; 79 80];

for idd = 1:size(d1,1)
    d2dirs = eids(idd,:);
    
    try         
    tic
        
    %copy over positions and laps from both
    
    for itrack = 1:2
        cd(dirs.spikedatadirOrig)        
        load(d2(d2dirs(itrack)).name,'linposcat','MidTime','pos','vel','dirdat','params','cm_conv') 
        params1.Novel = params.Novel;
        
        clear params        
        
        % make params
        params.Novel = params1.Novel;
        params.ident = [d1(idd).name '_' num2str(itrack)];
        params.converted = d1(idd).name;
        if idd== 3 || idd >4
            params.converted = [d1(idd).name '/Run' num2str(itrack)];
        end
        params.daydir = dirs.remotedatadir;
        params.Run_Times = [min(pos(:,1)) max(pos(:,1))];
        params.HPtt = 1:40;
        params.Othertt = [];        
        params.Rat_Name = d2(d2dirs(1)).name(1);
        params.Track_Type = 1;
        save([dirs.spikedatadir '\' params.ident '.mat'],'params')
        
        % get clusters out and into a spikedata with times and clusters
        extract_spikes_from_xclust(dirs,params)
        load([dirs.spikedatadir '\' params.ident '.mat'],'rawspikedata')

        %indexes the position bins of each spike, seperates types of neurons
        [spikedata,hp_cells,other_cells,hpinterneurons] = GUI_identify_interneurons(rawspikedata,pos,params);
        save([dirs.spikedatadir params.ident],'spikedata','hp_cells','other_cells','hpinterneurons',...
            'linposcat','MidTime','pos','vel','dirdat','cm_conv','-append')
            
    end
    
     tt = toc;
    disp(['Done with # ' num2str(idd) ' ' d2(idd).name ', done in ' num2str(round(tt/60,3,'significant')) ' minutes'])
    
    catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(idd)
        disp('****************************************Error Occured')
    end
end
end
%% Then adding the extra Remote Sleep sessions I found and clustered
if 1
    
clearvars -except dirs
cd(dirs.remotesleepload)
d2 = dir('*.mat');
ses = [1 3:6 9]; %sessions to copy over
for id = 1:length(ses)
    idd = ses(id);
    
    load(d2(idd).name,'spikedata','hp_cells','other_cells','hpinterneurons',...            
            'linposcat','MidTime','pos','vel','dirdat','params','cm_conv') 
        
    st = (repmat(pos(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(pos(:,1))]);
    nd = (repmat(pos(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(pos(:,1))]);
    rundat = (st&nd)';
    pos1 = pos; vel1 = vel; linposcat1 = linposcat; dirdat1 = dirdat;
    spikedata1 = spikedata; MidTime1 = MidTime;        
    params1 = params;
    
    figure; hold on
    %break up into two sessions like the other format
    for itrack = 1:2    
        
        clear pos vel linposcat dirdat spikedata MidTime params
        % make params
        params.Novel = 1;
        params.ident = [d2(idd).name(1:end-4) '_' num2str(itrack)];        
        params.daydir = dirs.remotesleepload;
        params.Run_Times = params1.Run_Times(itrack,:);
        params.HPtt = 1:40;
        params.Othertt = [];        
        params.Rat_Name = d2(idd).name(1);
        params.Track_Type = 1;
        
        % pos, vel, linposcat, dirdat
        pos = pos1(rundat(:,itrack),:);        
        vel = vel1(rundat(:,itrack),:);      
        linposcat = linposcat1(rundat(:,itrack),:);
        dirdat = dirdat1(rundat(:,itrack),:);
        % spikedata
        spikedata = spikedata1(ismember(spikedata1(:,3),find(rundat(:,itrack))),:);
        spikedata(:,3) = spikedata(:,3)-find(rundat(:,itrack),1,'first')+1;
        % MidTime
        MidTime = MidTime1(MidTime1>=min(pos(:,1)) & MidTime1<=max(pos(:,1)));
%         if exist([dirs.spikedatadir params.ident '.mat'],'file')
%             error('overwrite')
%         end
        subplot(2,1,itrack);
        hold on;
        plot(pos(:,1),linposcat(:,1),'k'); hold on; plot(pos(spikedata(:,3),1),linposcat(spikedata(:,3),1),'r.')
         yl = get(gca,'ylim');
         plot([MidTime MidTime]',[yl.*ones(size(MidTime,1),2)]','b-')
       
        title([num2str(itrack) ' - ' num2str(length(MidTime)-1) ' laps'])

        save([dirs.spikedatadir params.ident '.mat'],'spikedata','hp_cells','other_cells','hpinterneurons',...
            'linposcat','MidTime','pos','vel','dirdat','cm_conv','params')
    end
    suptitle([params.ident])
    helper_saveandclosefig([dirs.homedir 'Figures\remote\behavior_' params.ident(1:end-2)])
end
end

%% Then running decoding on all of them
if 1


clearvars -except dirs
cd(dirs.spikedatadir)
d2 = dir('*.mat');
testing = extractfield(d2,'name'); 
testing = erase(testing,'_1.mat'); 
testing = erase(testing,'_2.mat'); 
testing = erase(testing,'_3.mat');
d1 = unique(testing);
torun = [1 7 9 14]; %[1 5:7 9 14];
for idd = 1:size(d1,2)
    if ~ismember(idd,torun)
        disp(['Skip ' num2str(idd)])
        continue
    end
    
    try 
    tic
    
    if 1 %%%% plot to check
    cd(dirs.spikedatadir)
    
    load([d1{idd} '_1'],'hp_cells')
    [max(hp_cells) size(hp_cells)]
    load([d1{idd} '_2'],'hp_cells')
    [max(hp_cells) size(hp_cells)]
    load([d1{idd} '_1'],'spikedata')
    spikedata1 = spikedata;
    load([d1{idd} '_2'],'spikedata')
    figure; histogram(spikedata1(:,2),'Normalization','probability'); 
    hold on; histogram(spikedata(:,2),'Normalization','probability')
    
    end
    
    if 1
    for itrack = 1:2
        cd(dirs.spikedatadir)
        load([dirs.spikedatadir '\' d1{idd} '_' num2str(itrack) '.mat'],'params')   
        
        %decode 2 runs (two combinations of epoch and fields)
        decode_spikedensity_events(params.ident,false)
        shuffle_wc_replayscore(params.ident)

        get_hover_jumps(params.ident)
    end   
    end
    
    


    if 1
        %%%%% decode 3rd combinations of epoch and fields 
        %%%%% (e1f1, e2f2 inside loop above, now do e2f1)
        cd(dirs.spikedatadir)
        load([d1{idd} '_2'],'params')
        
        %%%%% saves out 'OutFRlaps','InFRlaps','CandSeq','OutFR',
        %%%%% 'InFR','Spike','Index','OutMatrix','InMatrix','Cell_Number'
        decode_spikedensity_events_remote(d1{idd})
        disp('Done with decode_spikedensity_events')
        
        %%%%% saves out 'CandCorr','PosShuffleCorr','CandDis',
        %%%%% 'PosShuffleDis','PosShift','CandRS'
        shuffle_wc_replayscore([d1{idd} '_3'])        
        disp('Done with shuffle_wc_replayscore')
        

        %%%%% saves out  'hover1','move1','hover2','move2','CandStepSize',
        %%%%% 'hovertime1', 'movetime1', 'hovertime2', 'movetime2'
        get_hover_jumps([d1{idd} '_3'])
        disp('Done with get_hover_jumps')
    end
    
    tt = toc;
    disp(['Done with # ' num2str(idd) ' ' d1{idd}...
        ', done in ' num2str(round(tt/60,3,'significant')) ' minutes'])
    
    catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(idd)
        disp('****************************************Error Occured')
    end
end

end
