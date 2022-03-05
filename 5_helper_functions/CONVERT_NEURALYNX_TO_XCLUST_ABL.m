function CONVERT_NEURALYNX_TO_XCLUST_ABL(Tetrodes_To_Convert,Video_String,Start_Timepoint,End_Timepoint)

%%%%% This is an edited version of Brad Pfieffer's code to convert
%%%%% neurolynx to xclust

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% 
% This program will load the Neuralynx tetrode spike data and video
% position data and convert this into a file that Xclust can load.
% 
% Tetrodes_To_Convert is a numerical list of the tetrode numbers.
% Video_String is a string that names the raw position data from Neuralynx.
% Start_Time and End_Time are single numbers that mark the start and end
% timepoints (in seconds) to pull the spikes from.
% 
% Example:  CONVERT_NEURALYNX_TO_XCLUST(1:10,'VT1.nvt',5600,8745);
%    This would convert tetrodes 1 through 10 using the VT1.nvt file for
%    position data, and would only convert spikes that happen between times
%    5600 and 8745 (the times are in seconds even though Neuralynx records
%    in microseconds -- this will be accounted for later in the program).
% 
% To convert all of the timepoints, use the following line:
%          CONVERT_NEURALYNX_TO_XCLUST(1:10,'VT1.nvt','start','end');
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% If the user did not fully define the input variables, the most likely
% values are either searched for or automatically provided.
if nargin<1,
    Tetrodes_To_Convert=1:40;
    Video_String='VT1.nvt';
    Start_Timepoint='start';
    End_Timepoint='end';
elseif nargin<2,
    Video_String='VT1.nvt';
    Start_Timepoint='start';
    End_Timepoint='end';
elseif nargin<3,
    Start_Timepoint='start';
    End_Timepoint='end';
elseif nargin<4,
    End_Timepoint='end';
end

% This section loads the position information so that this informaton can
% be attributed to the spikes prior to clustering.  If the
% Load_Raw_Position_Information function has been run, there will be a
% variable called "Position_Data", otherwise, the program will load the raw
% position information and process it a little before attributing these
% values to the spikes.
if exist('Original_Position_Data.mat')==2,
    load('Original_Position_Data','Position_Times','X_Positions','Y_Positions')
else
    [Position_Times,X_Positions,Y_Positions]=Nlx2MatVT(Video_String,[1 1 1 0 0 0],0,1);
end
Position_Times=Position_Times/1000000;        % Time is now in seconds
Good_Position_Data_Points=find(Position_Times>0 & X_Positions>0 & Y_Positions>0);
Corrected_Position_Times=Position_Times(Good_Position_Data_Points);
Corrected_X_Positions=X_Positions(Good_Position_Data_Points);
Corrected_Y_Positions=Y_Positions(Good_Position_Data_Points);
Position_Data=[Corrected_Position_Times',Corrected_X_Positions',Corrected_Y_Positions'];
% This section smooths the X- and Y-position information to correct for the
% inherent 'jitter' in how Neuralynx identifies position.
[FilterA,FilterB]=butter(2,0.05);
Position_Data(:,2)=filtfilt(FilterA,FilterB,Position_Data(:,2));
Position_Data(:,3)=filtfilt(FilterA,FilterB,Position_Data(:,3));
Position_Data=sortrows(Position_Data,1);
clear FilterA;
clear FilterB;
clear Frame_Rate;
clear Single_Dislocations;
clear Multiple_Dislocations;
clear String;
clear P;
clear Dislocated_Position_Times;
clear Position_Times;
clear X_Positions;
clear Y_Positions;
clear FilterA;
clear FilterB;
clear Dislocated_Position_Times;


for Tetrode_Number=Tetrodes_To_Convert,
    
    % This next section is just setting up the names of files and directories
    % to write the data to once it's been converted.  It also makes sure that
    % the tetrode file exists and has some real data in it.
    Tetrode_Name=sprintf('TT%d.ntt',Tetrode_Number);
    Tetrode_Directory=sprintf('Converted/tt%d',Tetrode_Number);
    Tetrode_File_Name=sprintf('tt%d.xclust',Tetrode_Number);
    if ~exist(Tetrode_Name,'file')==2;
        continue
    end
    Directory_List=dir;
    
    %This section finds the file that matches the name of the tetrode to load
    for T=1:size(Directory_List,1),
        if strcmp(Directory_List(T).name,Tetrode_Name),
            break
        end
    end
    
    %If the tetrode file is 20kb or smaller, don't bother loading it
    %because it might be completely empty and this will break the program
    if Directory_List(T).bytes<20000,
        continue
    end
    clear Directory_List;
    clear T;
    
    % This next section pulls up the Neuralynx header information to determine
    % the parameters under which the data was stored.  This will give a value
    % termed "Bit_Conversion."  If you divide the amplitude of any response by
    % this value, you will get the actual uV of the response.
    [Header]=Nlx2MatSpike(Tetrode_Name,[0 0 0 0 0],1,3,1);
    
    InputRow = strfind(Header,'-InputRange');
    InputRowConv = ~cellfun(@isempty,InputRow);
    Input_Line=cell2mat(Header(InputRowConv));
    eval(sprintf('Bit_Conversion=32768/%s;',Input_Line(13:16)));
    clear Header;
    clear Input_Line;    
    
    % This next section loads up the spike data using a program that Neuralynx
    % made called 'Nlx2MatSpike.'  It loads the spike times and the spike
    % parameters (eight values per spike that represent the max height and
    % minimum valley for each of the four tetrodes during that particular
    % spike).  You can read more about the Nlx2MatVT program in its help file.
    % The beginning and end of the spike times are "corrected" so that there
    % are no timepoints before the first position timepoint or after the
    % last position timepoint.
    %
    % This is now done below as the program loads the spike waveforms.
    %
    % This next section loads the spike waveforms using the same file as above,
    % Nlx2MatSpike.  The spike waveforms are each written as only 32 data
    % points per second (for a sampling rate of only 32 Hz rather than the 32
    % kHz that Neuralynx actually stores the data as).  These are only useful
    % for visualizing the data, probably not so useful for actual analysis.
    
    
    [Times]=Nlx2MatSpike(Tetrode_Name,[1 0 0 0 0],0,1);
    if strcmp(Start_Timepoint,'start'),
        Start_Time=min(Times);
    else
        Start_Time=Start_Timepoint*1000000;  %This is necessary to convert the Start_Time from seconds to microseconds, which is how Neuralynx stores the time
    end
    if strcmp(End_Timepoint,'end'),
        End_Time=max(Times);
    else
        End_Time=End_Timepoint*1000000;  %This is necessary to convert the End_Time from seconds to microseconds, which is how Neuralynx stores the time
    end
    clear Spike_Times;
    clear Max_Height;
    clear Max_Width;
    clear Spike_Parameters;
    for N=Start_Time:100000000:End_Time, %I can't load too much information at once or MatLab will crash, so I load the spike data in several large groups
        if ~isempty(find(Times>=N & Times<=N+99999999)),
            %Temp_Spike_Times is the timestamp of each spike, Temp_Spike_Parameters lists the peaks and minima of each spike on all four tetrodes 
            [Temp_Spike_Times,Temp_Spike_Parameters]=Nlx2MatSpike(Tetrode_Name,[1 0 0 1 0],0,4,[N min(End_Time,N+99999999)]);
            if exist('Spike_Times','var')
                Spike_Times=[Spike_Times,Temp_Spike_Times];
                Spike_Parameters=[Spike_Parameters,Temp_Spike_Parameters];
            else
                Spike_Times=Temp_Spike_Times;
                Spike_Parameters=Temp_Spike_Parameters;
            end
            if ~isempty(find((Spike_Times>=N)&(Spike_Times<=N+99999999))),
                if N<End_Time,
                    [Spike_Waves]=Nlx2MatSpike(Tetrode_Name,[0 0 0 0 1],0,4,[N min(End_Time,N+99999999)]);
                    
                    % This next bit of code determines the height and width of the tallest and
                    % widest spike of the four tetrodes at each spike.  The height is
                    % determined by the total vertical distance (microvolts) from the peak to
                    % the valley of each spike (and then picking the largest of those four for
                    % each tetrode).  The width is determined by by total horizontal distance
                    % (milliseconds) from the peak to the valley of each spike (and then
                    % picking the largest of those four for each tetrode).
                    
                    clear Temp_Max_Height;
                    clear Temp_Max_Width;
                    clear Max_Value;
                    clear Max_Index;
                    clear Min_Value;
                    clear Min_Index;
                    Temp_Max_Height=reshape(max(max(Spike_Waves)-min(Spike_Waves)),size(Spike_Waves,3),1);
                    [Max_Value,Max_Index]=max(Spike_Waves);
                    [Min_Value,Min_Index]=min(Spike_Waves);
                    Temp_Max_Width=reshape(max(Min_Index-Max_Index),size(Max_Index,3),1);
                    if(~exist('Max_Height'))
                        Max_Height=Temp_Max_Height;
                    else
                        Max_Height=[Max_Height;Temp_Max_Height];
                    end
                    if(~exist('Max_Width'))
                        Max_Width=Temp_Max_Width;
                    else
                        Max_Width=[Max_Width;Temp_Max_Width];
                    end
                end
            end
        end
    end
    clear Times;
    clear N;
    Spike_Times=Spike_Times/1000000;  %This converts the timestamps of each spike from microseconds to seconds
    
    % The following section looks for any "dislocated" spikes. A dislocated
    % spike is one that has a timestamp that appears "out of order"
    % relative to the other spikes. I don't delete these data points.  This
    % happens occasionally when the hard drive is too full or too much data
    % is being written at once.
    Dislocated_Spike_Times=find(diff(Spike_Times)<=0);
%     Single_Dislocations=0; %ABL changed from Brad's
%     Multiple_Dislocations=0;
    if ~isempty(Dislocated_Spike_Times),
        disp(sprintf('Found %d spike time dislocations in the spike data',length(Dislocated_Spike_Times)))
    end
    clear Single_Dislocations;
    clear Multiple_Dislocations;
    clear N;
    clear Temp_Max_Width;
    clear Temp_Max_Height;
    clear Max_Value;
    clear Max_Index;
    clear Min_Value;
    clear Min_Index;
    clear Temp_Spike_Times;
    clear Temp_Spike_Parameters;
    clear String;
    clear Dislocated_Spike_Times;
    
    
    % The following section determines the X- and Y-position information for
    % each spike by determining the closest timestamp in the position data and
    % pulling out the X- and Y-position data based on that timestamp.  This is
    % only used in the xclust2 program -- not in the actual analysis to
    % determine place fields.  This will just help you to see where your
    % cluster fires on the track, but this data won't be used later on in
    % subsequent analysis programs to actually calculate place fields.  I use a
    % more accurate algorithm for that.
    
    [~,Spike_Times_In_Position_Time]=histc(Spike_Times(1,:),Position_Data(:,1)); %ABL changed 1/27/2020
%     Histogram=histc(Position_Data(:,1),Spike_Times(1,:));
%     Spike_Times_In_Position_Time=cumsum(Histogram);
    Spike_Times_In_Position_Time(find(Spike_Times_In_Position_Time==0))=1;
    Matched_Spike_Times=Position_Data(Spike_Times_In_Position_Time,1);
    Spike_X_Positions=Position_Data(Spike_Times_In_Position_Time,2);
    Spike_Y_Positions=Position_Data(Spike_Times_In_Position_Time,3);
    clear Spike_Times_In_Position_Times;
    clear Matched_Spike_Times;
    
    % This section writes the information to the file that Xclust can read.
    Current_Directory=pwd;
    if(~exist(Tetrode_Directory,'dir'))
        mkdir(Tetrode_Directory);
    end
    cd(Tetrode_Directory);
    File_ID=fopen(Tetrode_File_Name,'w');
    File_String=sprintf('%%%%BEGINHEADER\n%% File type:    Binary\n%% Fields:     id,3,4,1     t_px,2,2,1     t_py,2,2,1     t_pa,2,2,1     t_pb,2,2,1     t_maxwd,2,2,1     t_maxht,2,2,1     time,5,8,1     pos_x,2,2,1     pos_y,2,2,1     velocity,4,4,1\n%%%%ENDHEADER\n');
    Count=fwrite(File_ID,File_String);
    for I=1:size(Spike_Times,2)
        Count=fwrite(File_ID,I,'uint32');
        Count=fwrite(File_ID,Spike_Parameters(1,I),'uint16');
        Count=fwrite(File_ID,Spike_Parameters(2,I),'uint16');
        Count=fwrite(File_ID,Spike_Parameters(3,I),'uint16');
        Count=fwrite(File_ID,Spike_Parameters(4,I),'uint16');
        Count=fwrite(File_ID,Max_Width(I),'uint16');
        Count=fwrite(File_ID,Max_Height(I),'uint16');
        Count=fwrite(File_ID,(Spike_Times(1,I)),'float64');
        Count=fwrite(File_ID,Spike_X_Positions(I),'uint16');
        Count=fwrite(File_ID,Spike_Y_Positions(I),'uint16');
        Count=fwrite(File_ID,0,'float32');
    end
    fclose(File_ID);
    cd(Current_Directory)
    disp(sprintf('Finished tetrode %d of %d',Tetrode_Number,max(Tetrodes_To_Convert)));
end

end







