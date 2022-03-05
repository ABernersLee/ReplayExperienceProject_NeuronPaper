function RunReplay_ABL_grosmark(basedir,plot_paper_figures,plot_other_figures)

%%%%% bring grosmark data in to look at INT

load([basedir 'dirs_linear_grosmark'],'dirs')
cd(dirs.spikedatadir)
d2 = dir('*.mat');

%% get data and process it, saving it out and checking with figures

%%%% circle tracks: 2, 5, 8. Others are linear.

for id = 1:size(d2,1)          
    if id == 2 || id==5 || id==8
        continue
    end
    
    try         
    tic
        
    if 0 % THIS WRITES OVER, this is changing grosmark data into my format
        cd(dirs.spikedatadirOrig)               
        load(d2(id).name)
        clear params        
        % make params
        params.Novel = true;
        params.Track_Type = 1;
        params.ident = d2(id).name(1:end-13);                        
        params.Run_Times = sessInfo.Epochs.MazeEpoch;
        %get data into my format                        
        rawpos = [sessInfo.Position.TimeStamps' sessInfo.Position.TwoDLocation.*100]; %m to cm        
        linposcat = sessInfo.Position.OneDLocation*100; %m to cm
        params.armslength = range(linposcat);   
        params.MazeType = sessInfo.Position.MazeType;linposcat = sessInfo.Position.OneDLocation*100; %m to cm
        params.armslength = range(linposcat);   
        %make vel
        [pos,params,vel,deltm] = GUI_make_arms_cleanup_posdata(rawpos,params);
        pos(:,2) = fillmissing(pos(:,2),'linear');
        pos(:,3) = fillmissing(pos(:,3),'linear');         
        vel = get_velocity(pos);                             
        linposcat(deltm) = [];        
        linposcat = fillmissing(linposcat,'nearest');
        
        skip = 50;
        movement = NaN(size(linposcat,1),1);
        for j = 1:skip-1       
            dat = linposcat(j:skip:end-j-1);
            movement(j:skip:end-j-1)=diff([dat;dat(end,:)]);                        
        end

        dirdat2 = NaN(size(linposcat,1),1);
        dirdat2(movement>0) = 2; %true; %new
        dirdat2(movement<0) = 1; %false; %new

        dirdat2(isnan(dirdat2)) = 1.5; %new
        [FilterA] =fspecial('gaussian',[20 1],10); %new 30 10
        Smoothed_Velocity=filtfilt(FilterA,1,dirdat2); %new    
        vel2 = round(Smoothed_Velocity)==2; %new
        dirdat = NaN(size(linposcat));
        dirdat(~isnan(linposcat)) = vel2(~isnan(linposcat));
        
        
        figure; hold on; 
        subplot(4,1,1)
        plot(pos(:,2),pos(:,3),'.k')
        title(num2str(id))
        subplot(4,1,2); hold on
        plot(pos(:,1),linposcat,'.k')
        plot(pos(vel<5,1),linposcat(vel<5),'.r')
        subplot(4,1,3); hold on
        plot(pos(:,1),vel,'k')
        plot([min(pos(:,1)) max(pos(:,1))],5,'r--')
        plot(pos(vel<5,1),vel(vel<5),'b.')
        subplot(4,1,4); hold on;
        plot(pos(dirdat==1,1),linposcat(dirdat==1),'r.')
        plot(pos(dirdat==0,1),linposcat(dirdat==0),'b.')
        plot(pos(isnan(dirdat),1),linposcat(isnan(dirdat)),'k.')
        set(gcf,'Position',[1921        -799        1080        1803])
        helper_saveandclosefig(['D:\test\Figures\grosmark\behavior' params.ident])
        
        
        
        %spikes
        hp_cells = sessInfo.Spikes.PyrIDs;
        hpinterneurons = sessInfo.Spikes.IntIDs;
        spikedata = [sessInfo.Spikes.SpikeTimes sessInfo.Spikes.SpikeIDs];
        other_cells = [];        
        params.Run_Times = [min(pos(:,1)) max(pos(:,1))];
        spikedata(spikedata(:,1)<params.Run_Times(1) | spikedata(:,1)>params.Run_Times(2),:) = [];
        Time = pos(:,1);
        [~,i] = histc(spikedata(:,1),[min([Time(1) min(spikedata(:,1))])-0.001; Time(1:end-1)+diff(Time)/2 ; max([Time(end) max(spikedata(:,1))])+0.001]);
        spikedata = [spikedata i];
                   
        save([dirs.spikedatadir '\' params.ident],'spikedata','params','hp_cells','other_cells','hpinterneurons',...
            'linposcat','pos','vel','dirdat')
    end
        
    if 0 %make MidTime semi-manually
        cd(dirs.spikedatadir) % using d2 with dirs.spikedatadir
        load(d2(id).name,'pos','linposcat','params')    
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',60);    
        if id==3
            locs = [locs; 27300];
            locs = sort(locs,'ascend');
        end
        MidTime = pos(locs,1);         
        if id==3       
            exl = 1:3;
        elseif id==4
            exl = [3,4,6,19,20,27,54,69];
        elseif id==5
            exl = 77;
        else
            exl = []; 
        end
        MidTime(exl) = [];
    %     figure; hold on
    %     aa = subplot(2,1,1); hold on;    
    %     plot(pos(:,1),linposcat,'k'); hold on; plot(MidTime(:),80,'r*')
    %     bb = subplot(2,1,2); 
    %     plot(pos(:,1),-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'k')
    %     linkaxes([aa bb],'x')
        figure; hold on; plot(pos(:,1),linposcat,'k'); hold on; plot(MidTime(:),80,'r*')
        set(gcf,'Position',[ 1065         285        1921         366])
        helper_saveandclosefig(['D:\test\Figures\grosmark\laps' params.ident])    
        save(d2(id).name,'MidTime','-append')
    end
        
    if 1  % this is running replays on it
        cd(dirs.spikedatadir)        
        load(d2(id).name,'params')
        cm_conv = 1;
        save(d2(id).name,'cm_conv','-append')
        decode_spikedensity_events(params.ident,false)
        shuffle_wc_replayscore(params.ident)
        get_hover_jumps(params.ident)
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
%% plot

RunReplay_ABL_grosmark_Figures(basedir,plot_paper_figures,plot_other_figures)