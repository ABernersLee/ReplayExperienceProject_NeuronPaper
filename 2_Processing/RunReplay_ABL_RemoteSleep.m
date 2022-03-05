%%%%%
%%%%% 
%%%%% This processes the remote sleep data
%%%%% ran ExtractRawData_ABL_RemoteSleep before
%%%%%
%%%%%


%% process data

clearvars
basedir = 'E:\test\';
load([basedir 'dirs_remotesleep'],'dirs')
cd(dirs.spikedatadir)
d2 = dir('*.mat');

for id = 1:size(d2,1)    
    daydir = d2(id).name;
    
    if 0
    load(daydir,'params','rawpos','rawspikedata','pos')    
    params.Track_Type = 1;
    params.armslength = 180;
    cm_conv = 1;
    
    %%%% adds arms and armslength to params and makes pos that has had 
    %%%% incorrect points deleted and linearly interpolated points 
    %%%% replacing them
    if ismember(id,[5:8 10:13])
        [pos,params,vel] = GUI_make_arms_cleanup_posdata_day_cm(rawpos,params,cm_conv);
        save([dirs.spikedatadir daydir],'pos','params','vel','-append')
    end
    
    if id==3 && 1 % only do once
        disp('fixing')
        load([dirs.spikedatadir daydir],'pos','params','vel')
        st = (repmat(pos(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(pos(:,1))]);
        nd = (repmat(pos(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(pos(:,1))]);
        rundat = (st&nd)';
        pos2 = pos(rundat(:,2),:);        
        d1 = diff(pos2(:,2));
        d2 = diff(pos2(:,3));
        skipped = pos2(:,3)>435 | pos2(:,2)<256;
        pos2(skipped,2:3) = NaN;        
        pos2(:,2) = fillmissing(pos2(:,2),'movmedian',80); pos2(:,3) = fillmissing(pos2(:,3),'movmedian',80);       
        for i = 1:16
            d1 = diff(pos2(:,2));
            d2 = diff(pos2(:,3));
            skipped = [d1>15 | d2>15; false];
            % to check with figure
%             if i>=10
%                 figure; plot(pos2(:,1),pos2(:,3),'k.'); hold on; 
%                 plot(pos2(skipped,1),pos2(skipped,3),'r.');
%                 title(num2str(i))
%             end
            pos2(skipped,2:3) = NaN; 
            pos2(:,2) = fillmissing(pos2(:,2),'movmedian',80); pos2(:,3) = fillmissing(pos2(:,3),'movmedian',80);            
        end
        
        Filter =fspecial('gaussian',[20 1],10); %new 30 10
        figure; plot(pos2(:,1),pos2(:,3),'k.');
        flt = filter2(Filter,pos2(:,3));
        flt([1:10 size(flt,1)-10:size(flt,1)]) = pos2([1:10 size(flt,1)-10:size(flt,1)],3);
        pos2(:,3) = flt;
        flt = filter2(Filter,pos2(:,2));
        flt([1:10 size(flt,1)-10:size(flt,1)]) = pos2([1:10 size(flt,1)-10:size(flt,1)],2);
        pos2(:,2) = flt;
        pos(rundat(:,2),:) = pos2;
        params.arms{2} = params.arms{2}';
        save([dirs.spikedatadir daydir],'pos','params','-append')
    end
    
    if id==1 && 1 % only do once
        load([dirs.spikedatadir daydir],'pos','params','vel')
        params.arms{2} = params.arms{2}';
    end
    
    %linearizes position data and assigns arms
    [armpos,linposcat,linposnorm,linposcatnan,dirdat,cm_conv,rundat] = GUI_make_linear_position_v2(pos,params);
    save([dirs.spikedatadir daydir],'armpos','dirdat','rundat','linposcat','linposnorm','linposcatnan','cm_conv','-append')
    
    
    %makes behavior epoch marticies (middle, stem, neck, lick area)
    params.Run = 1;
    [behavior1, behave_change_log1, behave_ind1] = GUI_make_behavior(linposcat,armpos,linposcatnan,params);
    params.Run = 2;
    [behavior2, behave_change_log2, behave_ind2] = GUI_make_behavior(linposcat,armpos,linposcatnan,params);
    behavior = [behavior1;behavior2];
    behave_change_log = [behave_change_log1;behave_change_log2];
    behave_ind = [behave_ind1;behave_ind2];
    save([dirs.spikedatadir daydir],'behavior','behave_change_log','behave_ind','-append')
        
    %makes laps, a few different ways for track type 2, see script
    [laps_coverspace,laps_twoarms,laps_singlepass,headingarm,error_correct] = GUI_makelaps(behavior, behave_change_log, pos,armpos, linposcat,params);
    save([dirs.spikedatadir daydir],'laps_coverspace','laps_twoarms','laps_singlepass','headingarm','error_correct','-append')
        
    end
    
    if 1
        %indexes the position bins of each spike, seperates types of neurons
        load(daydir,'params','pos','rawspikedata')   
        params.HPtt = 1:40; params.Othertt = [];
        [spikedata,hp_cells,other_cells,hpinterneurons] = GUI_identify_interneurons(rawspikedata,pos,params);
        save([dirs.spikedatadir daydir],'spikedata','hp_cells','other_cells','hpinterneurons','-append')
    
    end
    clearvars -except dirs d2 id
    disp(num2str(id))
end

cd ../


%% decode replays two ways


clearvars
basedir = 'E:\test\';
load([basedir 'dirs_remotesleep'],'dirs')
cd(dirs.spikedatadir)
d2 = dir('*.mat');

for id = 1:size(d2,1)
    
    try         
    tic
    daydir = d2(id).name;
    if 0 %make MidTime semi-manually
    
    
    cd(dirs.spikedatadir) % using d2 with dirs.spikedatadir
    load(d2(id).name,'pos','linposcat','params','rundat')    
    if id ==1
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',500);    
    elseif id ==5 
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-4,'MINPEAKDISTANCE',300);        
    elseif id ==8
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-1.5,'MINPEAKDISTANCE',800);    
    elseif id ==11
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',600);    
    elseif id ==13
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-4,'MINPEAKDISTANCE',300);    
    elseif id > 8 
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-1.5,'MINPEAKDISTANCE',500);    
    else
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',300);    
    end
    MidTime = pos(locs,1);       
%     if id == 1
%         [~,locs1]=findpeaks(-abs(linposcat(rundat(:,1))-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',200);
%         pos2 = pos(rundat(:,2),:);
%         [~,locs2]=findpeaks(-abs(linposcat(rundat(:,2))-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',500);
%         MidTime = [pos(locs1,1); pos2(locs2,1)];       
%     end
    
    
    exl = [];
    
    if id==1       
      exl = [size(MidTime,1)-2 size(MidTime,1)-7];
    elseif id==2
      exl = [17 24 45 48 51 52 55];
    elseif id==3
      exl = [19 27 33];
    elseif id==4
      exl = [4 9 12 15 20:22];
    elseif id==5    
      exl = [54 74]; 
    elseif id==6
      exl = [35];
    elseif id==7
      exl = []; 
    elseif id==8
      exl = [];
    elseif id==9
      exl = [];
    elseif id==10
      exl = [39 54]; %9 14 41 42 57];
    elseif id==11
      exl = [51:53]; %29 31 52 54:57 64:66];
    elseif id==12
      exl = [50 55 56]; %47 51:52 54:56 59:60 62:65]; %66
    elseif id==13
      exl = [];
    end
    
    MidTime(exl) = [];
    
    figure; hold on
    aa = subplot(2,1,1); hold on;    
    plot(pos(:,1),linposcat,'k'); hold on; plot(MidTime(:),min(linposcat)+range(linposcat)/2,'r*')
    bb = subplot(2,1,2); 
    plot(pos(:,1),-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'k')
    linkaxes([aa bb],'x')

    figure; hold on; plot(pos(:,1),linposcat,'k'); hold on; plot(MidTime(:),min(linposcat)+range(linposcat)/2,'r*')
    set(gcf,'Position',[ 1065         285        1921         366])
    helper_saveandclosefig([dirs.figdir '\laps_' daydir])    
    save(d2(id).name,'MidTime','-append')
    close all
    end
    
    
    
    
    if 0  % this is running replays on it
        cd(dirs.spikedatadir)        
        load(d2(id).name,'params')
        daydir = d2(id).name;
        
        if 1 %change to having lapbylap seperately for each 'map' 
            % CHANGE TO NONoverlapping bins? 
            disp('decode_spikedensity_events Run 1')
            [CandSeq1,OutFR1,InFR1,Index1 ,OutMatrix1, InMatrix1] = decode_spikedensity_events_remote_sleep(daydir,1);
            
            disp('decode_spikedensity_events Run 2')
            [CandSeq2,OutFR2,InFR2,Index2 ,OutMatrix2, InMatrix2] = decode_spikedensity_events_remote_sleep(daydir,2);
            
            if sum(sum(CandSeq1~=CandSeq2))>0 || sum(Index1~=Index2)
                disp('Error with decoding with each run field')
            end
            OutFR = cat(3,OutFR1,OutFR2);
            InFR = cat(3,InFR1,InFR2);            
            OutMatrix = cat(3,OutMatrix1,OutMatrix2);
            InMatrix = cat(3,InMatrix1,InMatrix2);
            Index = Index1; CandSeq = CandSeq1;
            save(daydir,'InFR','OutFR','InMatrix','OutMatrix','Index','CandSeq','-append')
        end
        
        
        disp('Starting shuffle_wc_replayscore')
        [CandCorr1,CandDis1,CandRS1] = shuffle_wc_replayscore_remotesleep(daydir,1);
        [CandCorr2,CandDis2,CandRS2] = shuffle_wc_replayscore_remotesleep(daydir,2);
        CandCorr = cat(2,CandCorr1,CandCorr2);
        CandDis = cat(2,CandDis1,CandDis2);
        CandRS = cat(3,CandRS1,CandRS2);
        save(daydir,'CandCorr','CandDis','CandRS','-append')
        
        disp('Starting get_hover_jumps')
        [hover1,move1,hovertime1,movetime1,CandStepSize1] = get_hover_jumps_remotesleep(daydir,1); %using move2 for all of them
        [hover2,move2,hovertime2,movetime2,CandStepSize2] = get_hover_jumps_remotesleep(daydir,2);
        hover{1,1} = hover1;
        hover{1,2} = hover2;
        move{1,1} = move1;
        move{1,2} = move2;
        hovertime{1,1} = hovertime1;
        hovertime{1,2} = hovertime2;
        movetime{1,1} = movetime1;
        movetime{1,2} = movetime2;
        
        CandStepSize = cat(3,CandStepSize1,CandStepSize2);
        save(daydir,'hover','move','hovertime','movetime','CandStepSize','-append')
        
%         disp('Starting shuffle_wc_replayscore')
%         shuffle_wc_replayscore_nonoverlap(daydir)                       
%         disp('doneday')
%         
%         disp('Starting get_hover_jumps')
%         get_hover_jumps_nonoverlap(daydir)
        
    end   

    tt = toc;
    disp(['Done with # ' num2str(id) ' ' d2(id).name ', done in ' num2str(round(tt/60,3,'significant')) ' minutes'])
    
    catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(id)
        disp('****************************************Error Occured')
    end
end

disp('Done')

%% Figures


clearvars
basedir = 'E:\test\';
plot_paper_figures = true;
plot_other_figures = true;
RunReplay_ABL_RemoteSleep_Figures(basedir,plot_paper_figures,plot_other_figures)