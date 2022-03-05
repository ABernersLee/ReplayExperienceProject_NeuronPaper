%%%% Whole session decoding of replay (not lap by lap)

clearvars
basedir = 'E:\test\';
load([basedir 'dirs_linear_lapbylap_addedIN_wholesession.mat'],'dirs')
cd(dirs.spikedatadirOrig)
d2 = dir('*.mat');

for  id = 1:size(d2,1) 
    try 
        
    tic
        
    cd(dirs.spikedatadir)
      
    %%%%% saves out 'OutFRlaps','InFRlaps','CandSeq','OutFR','InFR',
    %%%%% 'Spike','Index','OutMatrix','InMatrix','Cell_Number'    
    decode_spikedensity_events(d2(id).name,false)
    disp('Done with decode_spikedensity_events')
    
    
    %%%%% saves out 'CandCorr','PosShuffleCorr','CandDis',
    %%%%% 'PosShuffleDis','PosShift','CandRS'
    cd(dirs.spikedatadir)
    shuffle_wc_replayscore(d2(id).name)
    disp('Done with shuffle_wc_replayscore')
    
    
    %%%%% saves out 'StepGamma' and 'SeqGamma'
    cd(dirs.spikedatadir)
    get_gamma_step(d2(id).name,dirs)
    disp('Done with get_gamma_step')
    
    
    %%%%% saves out  'hover1','move1','hover2','move2','CandStepSize', 
    %%%%% 'hovertime1', 'movetime1', 'hovertime2', 'movetime2'
    cd(dirs.spikedatadir)
    get_hover_jumps(d2(id).name)
    disp('Done with get_hover_jumps')
    
%     cd(dirs.spikedatadir)
%     get_hover_jumps_small(d2(id).name)
%     disp('Done with get_hover_jumps_small')
%     
%     cd(dirs.spikedatadir)
%     get_hover_jumps_fwd(d2(id).name)
%     disp('Done with get_hover_jumps_fwd')
    
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
