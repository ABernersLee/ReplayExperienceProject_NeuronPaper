%%%%%  This script brings Brad's one lap data in to look at whether 
%%%%%  replays get longer without experience

%%%%%  Only run this script once to get data out and do primary analysis
%%%%%  (e.g. position), then run figures for paper and/or other ones

%% Gets raw data out
if 0
    
clearvars
basedir = 'E:\test\';
load([basedir 'dirs_linear_OneLap'],'dirs')

cd(dirs.spikedatadirOrig)
d2 = dir('*_*');
for id = 1:size(d2,1)   
    clear params
    params.Date = d2(id).name(end-7:end);
    params.Rat_Name = d2(id).name(1:end-9);
    params.converted = 'Converted';
    params.Track_Type = 1;
    params.daydir = [dirs.spikedatadirOrig d2(id).name];
    if id ==1 
        params.TTname = '234';
        params.Sleep_Times = [9426 12840; 12970 16640; 17970 20203];
        params.Run_Times = [12850 12970; 16650 17970];
    elseif id==2
        params.Sleep_Times = [780 4595; 5115 9453; 10636 13400];
        params.Run_Times = [4595 5115; 9463 10636];
    elseif id==3
    params.Sleep_Times = [9815 12315; 12380 16615; 17170 18150];
        params.Run_Times = [12325 12380; 16625 17170];
    elseif id == 4
        params.Sleep_Times = [19270 21355; 21413 26557; 27600 29500];
        params.TTname = '1234';
        params.Run_Times = [21365 21413; 26567 27600];
    elseif id==5
        params.Sleep_Times = [2222 3952; 3995 8050; 8690 12400];
        params.TTname = '1234';
        params.Run_Times = [3962 3995; 8060 8690];        
    end
   cd(dirs.spikedatadir)
   save(d2(id).name,'params','-append')
end
end

%% Saves out position and replays

clearvars
basedir = 'E:\test\';
load([basedir 'dirs_linear_OneLap'],'dirs')

cd(dirs.spikedatadir)
d2 = dir('*.mat');

nn = NaN(size(d2,1),5);
%report neuron numbers and timings
for id = 1:size(d2,1)
    load(d2(id).name,'params','hp_cells','hpinterneurons')
    nn(id,1) = length(hp_cells(~ismember(hp_cells,hpinterneurons)));
    nn(id,2) = range(params.Sleep_Times(1,:))./60;
    nn(id,3) = range(params.Run_Times(1,:))./60;
    nn(id,4) = range(params.Sleep_Times(2,:))./60;
    nn(id,5) = range(params.Run_Times(2,:))./60;
end

for id = 1:size(d2,1)          
    
    try         
    tic
     
    if 1
        if 1 %saves over things, position

            cd(dirs.spikedatadir)
            load(d2(id).name,'params')
            params.ident = d2(id).name(1:end-4);
            params.Othertt = [];
            params.HPtt = 1:40;

            cd([params.daydir])   
            if id~=1 % id ==2 || id== 3 %%used to be only 2 and 3, now all but 1 %ABL 11/18/19
                load('Original_Position_Data')
                if sum(Position_Times<0)>0
                    Position_Times(Position_Times<0) = NaN;
                    Position_Times = fillmissing(Position_Times,'linear');
                end
                Position_Data = [Position_Times'./1e6  X_Positions'./2.6 Y_Positions'./2.6];
            else
                load('Position_Data')
            end

            cd(dirs.spikedatadir)

            %get data into my format                        
            rawpos = Position_Data(:,1:3);
            params.armslength = 160;               

            pos = []; vel = []; linposcat = []; dirdat = [];
            for irun = 1:2
                rawposr = rawpos(rawpos(:,1)>params.Run_Times(irun,1) & rawpos(:,1)<params.Run_Times(irun,2),:);
                [pos1,~,vel1] = GUI_make_arms_cleanup_posdata(rawposr,params);
                  
                if id == 3 || id == 5
                    linposcat1 =fillmissing(pos1(:,2),'linear');   
                else
                    linposcat1 = fillmissing(pos1(:,3),'linear');
                end
                lind = pos1(:,1)>=params.Run_Times(irun,1) & pos1(:,1)<=params.Run_Times(irun,2);
                vel1 = get_velocity([pos1(lind,1) linposcat1(lind) ones(sum(lind),1)]);
                 
                skip = 50;
                movement = NaN(size(linposcat1,1),1);
                for j = 1:skip-1       
                    dat = linposcat1(j:skip:end-j-1);
                    movement(j:skip:end-j-1)=diff([dat;dat(end,:)]);                        
                end

                dirdat2 = NaN(size(linposcat1,1),1);
                dirdat2(movement>0) = 2; %true; %new
                dirdat2(movement<0) = 1; %false; %new

                dirdat2(isnan(dirdat2)) = 1.5; %new
                [FilterA] =fspecial('gaussian',[20 1],10); %new 30 10
                Smoothed_Velocity=filtfilt(FilterA,1,dirdat2); %new    
                vel2 = round(Smoothed_Velocity)==2; %new
                dirdat1 = NaN(size(linposcat1));
                dirdat1(~isnan(linposcat1)) = vel2(~isnan(linposcat1));
                
                if id == 2 
                    [FilterA] =fspecial('gaussian',[20 1],5);
                    linposcat2 = linposcat1;
                    while sum(abs(diff(linposcat2))>20)>0
                        linposcat2(find(abs(diff(linposcat2))>20)+1) = NaN;
                        linposcat2 = fillmissing(linposcat2,'previous');
                    end
                    linposcat2 = filtfilt(FilterA,1,linposcat2);
                    figure; 
                    plot(linposcat2,'k.'); hold on; 
                    plot(find(abs(diff(linposcat2))>20),linposcat2(find(abs(diff(linposcat2))>20)),'g*'); 
                    plot(find(abs(diff(linposcat2))>20)+1,linposcat2(find(abs(diff(linposcat2))>20)+1),'r*')                        
                    linposcat1 = linposcat2;
                end
                if (id==3 && irun==1)
                    linposcat2 = linposcat1;
                     while sum(abs(diff(linposcat2))>10)>0
                         linposcat2(find(abs(diff(linposcat2))>10)+1) = NaN;
                         linposcat2(find(abs(diff(linposcat2))>10)) = NaN;
                          linposcat2 = fillmissing(linposcat2,'movmean',30);
                     end

                      figure; 
                    plot(linposcat2,'k.'); hold on; 
                    plot(find(abs(diff(linposcat2))>10),linposcat2(find(abs(diff(linposcat2))>10)),'g*'); 
                    plot(find(abs(diff(linposcat2))>10)+1,linposcat2(find(abs(diff(linposcat2))>10)+1),'r*')    
                    linposcat1 = linposcat2;                 
                end
                pos = [pos;pos1]; vel = [vel;vel1]; linposcat = [linposcat;linposcat1]; dirdat = [dirdat;dirdat1];
            end




            figure; hold on; 
            subplot(4,1,1)
            plot(pos(:,2),pos(:,3),'.k')
            title(num2str(id))
            a(1) = subplot(4,1,2); hold on
            plot(pos(:,1),linposcat,'.k')
            plot(pos(vel<5,1),linposcat(vel<5),'.r')
            a(2) = subplot(4,1,3); hold on
            plot(pos(:,1),vel,'k')
            plot([min(pos(:,1)) max(pos(:,1))],5,'r--')
            plot(pos(vel<5,1),vel(vel<5),'b.')
            a(3) = subplot(4,1,4); hold on;
            linkaxes(a,'x')
            plot(pos(dirdat==1,1),linposcat(dirdat==1),'r.')
            plot(pos(dirdat==0,1),linposcat(dirdat==0),'b.')
            plot(pos(isnan(dirdat),1),linposcat(isnan(dirdat)),'k.')
            set(gcf,'Position',[1921        -799        1080        1803])
            helper_saveandclosefig(['D:\test\Figures\OneLap\behavior' params.ident])
            save([dirs.spikedatadir '\' params.ident],'params','linposcat','pos','dirdat','vel','-append')
        end

        if 1 %spikes, saves over  
            %spikes
            load(d2(id).name,'params')
            disp('Starting extract_spikes_from_xclust')
            extract_spikes_from_xclust(dirs,params)
        end
        
        if 1 % seperates interneurons
            cd(dirs.spikedatadir)
            load(d2(id).name,'rawspikedata','pos','params')
            disp('Starting GUI_identify_interneurons')
            [spikedata,hp_cells,other_cells,hpinterneurons] = GUI_identify_interneurons(rawspikedata,pos,params);

            save([dirs.spikedatadir '\' params.ident],'spikedata','hp_cells','other_cells','hpinterneurons','-append')
        end

        if 1 %make MidTime semi-manually
        cd(dirs.spikedatadir) % using d2 with dirs.spikedatadir
        load(d2(id).name,'pos','linposcat','params')    
        [~,locs]=findpeaks(-abs(linposcat-max(linposcat)/2-min(linposcat)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',60);    
        if id==3
            locs = [locs; 12358]; %; 27300];
            locs = sort(locs,'ascend');
        end
        MidTime = pos(locs,1);       

        if id==1       
            exl = [];
        elseif id==2
            exl = [1:4 54:58]; %1:6]; % [1:6 11:12 58:63]; 
        elseif id==3
            exl = []; %25; 
        elseif id==4
            exl = [12];
        elseif id==5    
            exl = [50 54 56]; %55; 
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
        helper_saveandclosefig(['D:\test\Figures\OneLap\laps' params.ident])    
        save(d2(id).name,'MidTime','-append')
        close all
        end
        
        if 1 % plot velocity to check
           load(d2(id).name,'linposcat','params','pos')
           vel = [];
           for irun = 1:2
               lind = pos(:,1)>=params.Run_Times(irun,1) & pos(:,1)<=params.Run_Times(irun,2);
               pre_thinned_dat = [pos(lind,1) linposcat(lind) ones(sum(lind),1)];
%               vel1 = get_velocity(pre_thinned_dat);
              
                [FilterA] =fspecial('gaussian',[100 1],10);
                AA = diff(sqrt(pre_thinned_dat(:,2).^2+pre_thinned_dat(:,3).^2))./diff(pre_thinned_dat(:,1));
                vel11=filtfilt(FilterA,1,[AA(1:100); AA; AA(end-99:end)]);
                vel1 = abs(vel11(101:end-100));
                
%                 figure; plot(vel1,'.')
                vel = cat(1,vel,[vel1; vel1(end)]);
           end
           save(d2(id).name,'vel','-append')
            
        end

    end
    
    if 1  % this is running replays on it and saves them out
         cd(dirs.spikedatadir)        
        load(d2(id).name,'params')
        if 0            
            cm_conv = 1;
            save(d2(id).name,'cm_conv','-append')
        end
        
        if 1 % All non-overlapping for significance        
            
            disp('decode_spikedensity_events_OneLap')
            decode_spikedensity_events_OneLap(params.ident) %makes the spikedensity seperatley for each epoch using full run (run2) fields     

            disp('Starting shuffle_wc_replayscore')
            shuffle_wc_replayscore_nonoverlap(params.ident)     
        
            disp('decode_spikedensity_events_OneLap_Epoch1')
            decode_spikedensity_events_OneLap_Epoch1(params.ident) %makes the spikedensity seperatley for each epoch using first run (run1) fields %          makes 'FREpochOne','MatrixEpoch1'                                            

            disp('Starting shuffle_wc_replayscore_Epoch1')
            shuffle_wc_replayscore_Epoch1(params.ident)        
        end
       
        if 1 % overlapping bins
            
            disp('Starting get_hover_jumps_Epoch1')
            get_hover_jumps_Epoch1(params.ident)

            %makes the spikedensity seperatley for each epoch using full run (run2) fields  
            disp('decode_spikedensity_events_OneLap_overlap')
            decode_spikedensity_events_OneLap_overlap(params.ident) 
                        
            %makes the spikedensity seperatley for each epoch using first run (run1) fields 
            %makes 'FREpochOne','MatrixEpoch1'             
            disp('decode_spikedensity_events_OneLap_Epoch1')
            decode_spikedensity_events_OneLap_Epoch1_overlap(params.ident)                                                
        
            %changed numshuffles from 500 to 5000 % havn't changed these to
            %Tm, so making very few shuffles
            disp('Starting shuffle_wc_replayscore_Epoch1')
            shuffle_wc_replayscore_Epoch1_overlap(params.ident)

            disp('Starting shuffle_wc_replayscore')
            shuffle_wc_replayscore(params.ident)     
        end
        disp('doneday')
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

%% then plots figures

clearvars
basedir = 'E:\test\';
load([basedir 'dirs_linear_OneLap'],'dirs')

plot_paper_figures = true;
plot_other_figures = true;
RunReplay_ABL_OneLap_Figures(basedir,plot_paper_figures,plot_other_figures)
