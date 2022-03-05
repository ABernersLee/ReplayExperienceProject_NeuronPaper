%%%%%% This is a series of GUIs to make and correct the position
%%%%%% and behavior data which is then saved out along with figures

load('D:\test\dirs_linear.mat','dirs')

cd(dirs.spikedatadir)
d = dir;
sz = extractfield(d,'bytes')>0;
ftype = contains(extractfield(d,'name'),'.mat');
d2 = d(sz&ftype);
clear sz ftype d

for isession = 1:size(d2,1)
    load([d2(isession).name],'params','rawpos','rawspikedata')    
    
%%%% adds arms and armslength to params and makes pos that has had 
%%%% incorrect points deleted and linearly interpolated points replacing
%%%% them
    [pos,params,vel] = GUI_make_arms_cleanup_posdata(rawpos,params);
    save([dirs.spikedatadir params.ident],'pos','params','vel','-append')
    

%%%% linearizes position data and assigns arms
    [armpos,linposcat,linposnorm,linposcatnan,dirdat,cm_conv] = ...
        GUI_make_linear_position(pos,params);
    save([dirs.spikedatadir params.ident],'armpos','dirdat','linposcat',...
        'linposnorm','linposcatnan','cm_conv','-append')
    
%%%% makes behavior epoch marticies (middle, stem, neck, lick area)
    [behavior, behave_change_log, behave_ind] = ...
        GUI_make_behavior(linposcat,armpos,linposcatnan,params);
    save([dirs.spikedatadir params.ident],'behavior',...
        'behave_change_log','behave_ind','-append')
        
%%%% makes laps, a few different ways for track type 2, see script
    [laps_coverspace,laps_twoarms,laps_singlepass,headingarm,...
        error_correct] = GUI_makelaps(behavior, behave_change_log, pos, ...
        armpos, linposcat,params);
    save([dirs.spikedatadir params.ident],'laps_coverspace',...
        'laps_twoarms','laps_singlepass','headingarm',...
        'error_correct','-append')
    
    
%%%% indexes the position bins of each spike, seperates types of neurons
    [spikedata,hp_cells,other_cells,hpinterneurons] = ...
        GUI_identify_interneurons(rawspikedata,pos,params);
    save([dirs.spikedatadir params.ident],'spikedata',...
        'hp_cells','other_cells','hpinterneurons','-append')

    clearvars -except dirs d2 isession
    disp(num2str(isession))
end

cd ../