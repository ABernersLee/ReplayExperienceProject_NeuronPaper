% Run_ReplayExperience


%%%% Update: 3/5/2022
%%%% ABL
%%%% I commented some of the code a bit to put it online and took out some
%%%% (but not all) of the parts of the code that have analysis/figures that
%%%% didn't end up in the Neuron paper.

%%%% Update: 1/26/2022
%%%% ABL
%%%% I changed the units in the slopes from spatial bins/temporal bins to
%%%% meters/sec, which I had acciedentally not done. I changed these
%%%% scripts: 
% - [x] Plot_Main_HoverJump_Figures
% - [x] plot_ripple_events_over_laps_replays
% - [x] RunReplay_ABL_OneLap_Figures
% - [x] RunReplay_ABL_remote_Figures_combineLocal
% - [x] RunReplay_ABL_RemoteSleep_Figures
%%%% and updated Figure 2e and Sup Figs 1-4.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%%%%% Used these scripts to extract raw data. Data came from lots of
%%%%% different people's experiments on different harddrives
    ExtractRawData_Nlynx_xclust_ABL20180304
    ExtractRawData_ABL_RemoteSleep

%%%%% Used this script to do early processing of data including checking
%%%%% that position data is good and assigning type of neuron
    MakeProcessedData_ABL20180304

basedir = 'E:\test\'; %could make this global
% If working with processed data change this to the folder the data is in.
% If using data from the raw form these directory files would have been 
% created when you extracted them.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data is in a few different places depending on what experiment it is.
% Below are the directories of the data along with the scripts that run on
% that data.


%%%%% Novel and Familiar - lapbylap
%%%%% Figures:

% Sessions with novel and familiar linear tracks where we decode replays
% using lap by lap fields (instead of averaged across the whole session)
% use the following directories and are analyzed using 
RunReplay_ABL_lapbylap

dirs.spikedatadir = [basedir 'spikedata\final_linear_lapbylap_addedIN\'];
dirs.cscdatadir = [basedir 'cscdata\'];
dirs.paramdir = [basedir 'params\'];
dirs.homedir = basedir;
dirs.figdir = [basedir 'Figures\'];
dirs.spikedatadirOrig = dirs.spikedatadir;
save([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')




%%%%% Novel and Familiar - whole session
%%%%% Figures:

% Sessions with novel and familiar linear tracks where we decode replays
% using fields averaged across the whole session
% use the following directories and are analyzed using 
RunReplay_ABL_wholesession

dirs.spikedatadir= [basedir 'spikedata\final_linear_lapbylap_addedIN_wholesession\'];
dirs.paramdir = [basedir 'params\'];
dirs.homedir = [basedir];
dirs.figdir = [basedir 'Figures\'];
dirs.cscdatadir = [basedir 'cscdata\'];
dirs.spikedatadirOrig = [basedir 'spikedata\final_linear_lapbylap_addedIN'];
save([basedir 'dirs_linear_lapbylap_addedIN_wholesession.mat'],'dirs')




%%%%% Remote sleep replay
%%%%% Figures:

% Sessions with two novel tracks, clustered across those two tracks so that
% the same cells indecies are used in both sessions, only sessions where
% there is a rest period in between. These are analyzed using
RunReplay_ABL_RemoteSleep

dirs.spikedatadir= [basedir 'RemoteSleepData\'];
dirs.homedir= basedir;
dirs.spikedatadirOrig = [basedir 'RemoteSleepData\rawdata\'];
dirs.figdir= [basedir 'Figures\RemoteSleep'];
save([basedir 'dirs_remotesleep'],'dirs')



%%%%% Remote replay
%%%%% Figures:

% Sessions with two novel tracks, clustered across those two tracks so that
% the same cells indecies are used in both sessions. These are analyzed using
RunReplay_ABL_remote

dirs.spikedatadir =  [basedir 'RemoteClusteredData\'];
dirs.paramdir= [basedir 'params\'];
dirs.homedir = basedir;
dirs.cscdatadir = [basedir 'cscdata\'];
dirs.spikedatadirOrig= [basedir 'spikedata\final_linear_lapbylap_addedIN_wholesession\'];
dirs.remotedatadir= [basedir 'RemoteClusteredData\RawData\'];
dirs.remotesleepload = [basedir 'RemoteSleepData\'];
save([basedir 'dirs_linear_remote'],'dirs')



%%%%% OneLap
%%%%% Figures:

%Sessions where rats ran 1-2 passes on a novel track, then were allowed to
%rest in the sleep box, then ran a full session on the same track. These
%are analyzed using
RunReplay_ABL_OneLap 
% and then 
% RunReplay_ABL_OneLap_Figures(basedir,plot_paper_figures,plot_other_figures)

dirs.spikedatadir = [basedir 'OneLapData\'];
dirs.homedir = [basedir];
dirs.figdir = [basedir 'Figures\OneLap\'];
dirs.spikedatadirOrig = [basedir 'OneLapData\original_data\'];
save([basedir 'dirs_linear_OneLap'],'dirs')




%%%%% Grosmark et al Data
%%%%% Figures:

%Grosmark and Buzsaki data to have more interneurons in novel enviornments.
%These data are analyzed using
RunReplay_ABL_grosmark 
%and RunReplay_ABL_grosmark_Figures

dirs.spikedatadir = [basedir 'GrosmarkData'];
dirs.homedir = [basedir];
dirs.spikedatadirOrig = [basedir 'GrosmarkData\RawData'];
dirs.figdir = [basedir 'Figures\grosmark\'];
save([basedir 'dirs_linear_grosmark'],'dirs')




%%%%% HP and PFC data
%%%%% Figures:

% Data collected by Xiaojing with simultaneously recorded HP and PFC
% neurons. These data are analyzed using
RunReplay_ABL_ymaze

dirs.spikedatadir = [basedir 'spikedata\ymaze_new'];
dirs.paramdir = [basedir 'params\'];
dirs.homedir = [basedir];
dirs.figdir = [basedir 'Figures\PFC\'];
dirs.cscdatadir = [basedir 'cscdata\'];
dirs.spikedatadirOrig = [basedir 'spikedata\workingon_linear_add'];
save([basedir 'dirs_ymaze_new.mat'],'dirs')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Then these are the scripts used to generate the figures %%%%%%%%%%

RunReplay_ABL_OneLap_Figures(basedir,1,0)
ExampleReplays_Position(basedir)
Plot_Main_HoverJump_Figures(basedir,0)
RunReplay_ABL_remote_Figures_combineLocal(basedir) 
RunReplay_ABL_RemoteSleep_Figures(basedir,1,0)
ExampleReplays_gamma(basedir)
brad_effect(basedir,1,0)
schematic(basedir)
Hovers_Moving_laps(basedir,1,0)
RunReplay_ABL_grosmark_Figures(basedir,1,0)
INT_overlaps(basedir,1,0)
PFC_tolongreplays(basedir,1)


plot_ripple_events_over_laps_replays(basedir)
Hovers_Moving_FlatPosterior(basedir,1,0)
Hovers_Moving_Shuffle(basedir,1,0)
Plot_Hist_CellID_And_Flat(basedir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Then these are the scripts for each figure (repeated scripts) %%%%%%%%%

%%%%% Figure 1 
RunReplay_ABL_OneLap_Figures(basedir,1,0)
ExampleReplays_Position(basedir)

%%%%% Figure 2 
ExampleReplays_Position(basedir) 
Plot_Main_HoverJump_Figures(basedir,0)

%%%%% Figure 3 
ExampleReplays_Position(basedir)
Plot_Main_HoverJump_Figures(basedir,0)
RunReplay_ABL_remote_Figures_combineLocal(basedir) 

%%%%% Figure 4
RunReplay_ABL_OneLap_Figures(basedir,1,0)
RunReplay_ABL_RemoteSleep_Figures(basedir,1,0)

%%%%% Figure 5
ExampleReplays_gamma(basedir)
brad_effect(basedir,1,0)
Plot_Main_HoverJump_Figures(basedir,0)
schematic(basedir)

%%%%% Figure 6
ExampleReplays_Position(basedir)
Hovers_Moving_laps(basedir,1,0)

%%%%% Figure 7
RunReplay_ABL_grosmark_Figures(basedir,1,0)
INT_overlaps(basedir,1,0)
PFC_tolongreplays(basedir,1)


%%%%% Supplementary Figure 1
Plot_Main_HoverJump_Figures(basedir,0)
Plot_Main_HoverJump_Figures(basedir,0.4)
plot_ripple_events_over_laps_replays(basedir)

%%%%% Supplementary Figure 2
Plot_Main_HoverJump_Figures(basedir,0)
plot_ripple_events_over_laps_replays(basedir)

%%%%% Supplementary Figure 3
Plot_Main_HoverJump_Figures(basedir,0)

%%%%% Supplementary Figure 4
Plot_Main_HoverJump_Figures(basedir,0)
RunReplay_ABL_remote_Figures_combineLocal(basedir)
RunReplay_ABL_RemoteSleep_Figures(basedir,1,0)

%%%%% Supplementary Figure 5
Plot_Main_HoverJump_Figures(basedir,0)

%%%%% Supplementary Figure 6
Hovers_Moving_FlatPosterior(basedir,1,0)
Hovers_Moving_Shuffle(basedir,1,0)
Plot_Hist_CellID_And_Flat(basedir)

%%%%% Supplementary Figure 7
% none