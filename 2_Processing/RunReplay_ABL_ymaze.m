%%%%% Lap by lap decoding of replay for the ymaze

clearvars
basedir = 'E:\test\';
load([basedir 'dirs_ymaze_new.mat'],'dirs')
cd(dirs.spikedatadir)
d = dir;
sz = extractfield(d,'bytes')>0;
ftype = contains(extractfield(d,'name'),'.mat');
d2 = d(sz&ftype);
clear sz ftype d

for id = 1:size(d2,1)
    
    try 
        tic

        cd('D:\test\spikedata\ymaze_new')
        load(d2(id).name,'params')
        disp([ id params.Novel params.Track_Type])


        %%%% need to make laps all the same (make lap indecies into midpoint to match tings)
        get_midpoint(d2(id).name)        
        %%%% saves out 'MidTime'
        disp('Done with get_midpoint')    

        
        %%%%% saves out 'CandSeq','OutFR','InFR','Spike','Index',
        %%%%% 'OutMatrix','InMatrix','Cell_Number'
        decode_spikedensity_events_lapbylap(d2(id).name,0)
        disp('Done with decode_spikedensity_events')
    
        load(d2(id).name,'OutMatrix')
        if size(OutMatrix,3)>1
            load(d2(id).name,'InMatrix','Index','CandSeq','params','replayparams')        
            [Index,replayarms,CandSeq,InMatrix,OutMatrix,armpos] = ...
                convert_joint_to_single_runforarm(OutMatrix,InMatrix,Index,CandSeq,params,replayparams);
            save(d2(id).name,'Index','replayarms','CandSeq','InMatrix','OutMatrix','armpos','-append')
        end  
        
        %%%%% saves out CandCorr and CandDis
        wc_only(d2(id).name)
        disp('Done with WC')
    
        
    
        get_gamma_step(d2(id).name,dirs)
        %%%%% saves out 'StepGamma' and 'SeqGamma' (don't think I used this
        %%%%% in the end)
        disp('Done with get_gamma_step')
        
%         get_gamma_step_pfc(d2(id).name,dirs)
%         disp('Done with get_gamma_step_pfc')
        
        get_ripple_events(d2(id).name,dirs)
        disp('Done with get_ripple_events')

        get_cellactivitybylap(d2(id).name)
        %%%%% saves out 'SpikePerEvent','ActRatio','SpikePerAct','FRPerCell'
        disp('Done with get_cellactivitybylap')


        get_cellactivitybylap_interneurons(d2(id).name)
        disp('Done with get_cellactivitybylap_interneurons')

        get_gammaprop(d2(id).name,dirs)
        %%%%% saves out 'GammaProp','PolySlope'
        disp('Done with get_slopes_gammaprop')


        get_hover_jumps(d2(id).name)
        %%%%% saves out  'hover1','move1','hover2','move2','CandStepSize'
        disp('Done with get_hover_jumps')


         get_cellactivitybylap_hoverjump(d2(id).name)
        %%%%% saves out 'SpikePerEventHover','ActRatioHover',
        %%%%% 'SpikePerActHover','FRPerCellHover','SpikePerEventMove',
        %%%%% 'ActRatioMove','SpikePerActMove','FRPerCellMove,'
        disp('Done with get_cellactivitybylap_hoverjump')
     
        get_cellactivitybylap_first80ms(d2(id).name)
        disp('Done with get_cellactivitybylap_first80ms') 
        %did 60 80 100 and 120
    
    
    tt = toc;
    disp(['Done with # ' num2str(id) ' ' d2(id).name...
        ', done in ' num2str(round(tt/60,3,'significant')) ' minutes'])
     catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(id)
        disp('****************************************Error Occured')
    end
end

disp('Done')