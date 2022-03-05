%%%%% Lap by lap decoding of replay
clearvars
basedir = 'E:\test\';
load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')

cd(dirs.spikedatadir)
d2 = dir('*.mat');


cd(dirs.spikedatadir)
for id = 1:size(d2,1) 
    
    try 
        
    tic
    cd(dirs.spikedatadir)
    disp(['Starting # ' num2str(id) ' ' d2(id).name])
        
    %%%%% saves out 'MidTime'
    get_midpoint(d2(id).name)            
    disp('Done with get_midpoint')
   
    %%%%% saves out 'OutFRlaps','InFRlaps','CandSeq','OutFR',
    %%%%% 'InFR','Spike','Index','OutMatrix','InMatrix','Cell_Number'
    decode_spikedensity_events_lapbylap(d2(id).name,false)
    disp('Done with decode_spikedensity_events')
   
    %%%%% saves out 'CandCorr','PosShuffleCorr','CandDis',
    %%%%% 'PosShuffleDis','PosShift','CandRS'
    shuffle_wc_replayscore(d2(id).name)
    disp('Done with shuffle_wc_replayscore')       
       
    %%%%% saves out 'StepGamma' and 'SeqGamma'
    get_gamma_step(d2(id).name,dirs)
    disp('Done with get_gamma_step')        
    
    
    %%%%% saves out 'SpikePerEvent','ActRatio','ActRatioRate',
    %%%%% 'SpikePerAct','SpikePerActRate','FRPerCell'
    get_cellactivitybylap(d2(id).name)
    disp('Done with get_cellactivitybylap')        
    
    %%%%% saves out 'SpikePerEventIN','ActRatioIN','SpikePerActIN',
    %%%%% 'ActRatioRateIN','SpikePerActRateIN','FRPerCellIN'
    get_cellactivitybylap_interneurons(d2(id).name)
    disp('Done with get_cellactivitybylap_interneurons')
            
    %%%%% saves out 'GammaProp','GammaPropHigh'
    get_gammaprop(d2(id).name,dirs)
    disp('Done with get_gammaprop')      
    
    
    %%%%% saves out 'Ripple_CandEvents','HP_Ripple'
    get_ripple_events(d2(id).name,dirs)
    disp('Done with get_ripple_events')
    
    %%%%% saves out 'rp_freq'   
    set_ripple_freq(d2(id).name)
    disp('Done with get_ripple_freq')
    
    %%%%% saves out  'hover1','move1','hover2','move2','CandStepSize', 
    %%%%% 'hovertime1', 'movetime1', 'hovertime2', 'movetime2'
    get_hover_jumps(d2(id).name)
    disp('Done with get_hover_jumps')
   
    %%%%% saves out 'SpikePerEventHover','ActRatioHover',
    %%%%% 'SpikePerActHover','FRPerCellHover','SpikePerEventMove',
    %%%%% 'ActRatioMove','SpikePerActMove','FRPerCellMove,'
    get_cellactivitybylap_hoverjump(d2(id).name)
    disp('Done with get_cellactivitybylap_hoverjump')
    
    %%%%% saves out 'SpikePerEventFirst80','ActRatioFirst80',
    %%%%% 'SpikePerActFirst80','ActRatioRateFirst80',
    %%%%% 'SpikePerActRateFirst80','FRPerCellFirst80'
    get_cellactivitybylap_first80ms(d2(id).name)
    disp('Done with get_cellactivitybylap_first80ms') 
    %%%%% have looked at 60 80 100 and 120 (none of this analysis is in the
    %%%%% paper, nothing came of it)

    %%%%% smaller version    
    get_hover_jumps_small(d2(id).name)
    disp('Done with get_hover_jumps_small')

    %%%%% only using forward movement
    get_hover_jumps_fwd(d2(id).name)
    disp('Done with get_hover_jumps_fwd')
    
    
    load(d2(id).name,'params')
    if params.Novel==1
        %do nonoverlapping bins to look at significance over laps
        disp('Starting decode_spikedensity_events_lapbylap_NonOverLapping')
        decode_spikedensity_events_lapbylap_NonOverLapping(d2(id).name)
        disp('Done with decode_spikedensity_events_lapbylap_NonOverLapping')
    end
    tt = toc;
    
    disp(['Done with # ' num2str(id) ' ' d2(id).name ...
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