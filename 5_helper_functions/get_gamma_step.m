function get_gamma_step(iident,dirs)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% this is ting's original script, inspired by Brad's gamma work
% %% this is to replicate his result: step size phase locked to low gamma
% %% test it on different subgroup of data: trajectory events, correlated
% %% events, noise, low spike events
% %% and potentially test the strength of gamma associated with increased
% %% experience

% %% but over here this is just to add phases for each jump of the following
% %% Calculate the gamma phase of each decoded bin - in this case
% %% we use the average phase of the corresponding tetrode of each cell of
% %% the two consecutive decoding bins (as in Brad's paper)

%%%%% Used in RunReplay_ABL_wholesession.m, RunReplay_ABL_lapbylap.m, and
%%%%% RunReplay_ABL_ymaze.m

cd(dirs.spikedatadir)
load(iident,'CandSeq','hp_cells','hpinterneurons','rawspikedata','Spike','Index','OutMatrix','InMatrix')

% %% gamma phase for the overlapping of each pair of consective bins
cd(dirs.cscdatadir)

cellstouse = setdiff(hp_cells,hpinterneurons);
Tetrode = unique(rawspikedata(ismember(rawspikedata(:,2),cellstouse),3));


Spike = Spike(ismember(Spike(:,1),cellstouse),:);    

Orig_Spike= Spike;
EstBin=0.02;

Cand=CandSeq;
N=round(diff(CandSeq')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];  

[~,I]=histc(Spike(:,2),sortrows(Cand(:))); 
B=cumsum([0 diff(Cand')]);
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1);
end
Spike(~mod(I,2),:)=[];

[~,II]=histc(Spike(:,2),0:EstBin/4:Index(end)*EstBin); 

Spike=Orig_Spike;
Spike(~mod(I,2),:)=[];

SpikeGamma=Spike(:,1); % gamma phase of each spike   
% LFP_Frequency=3255;
EEGList=dir([iident(1:6) 'tt*' iident(6:end) '.mat']);
% Tetrode

for No=1:length(Tetrode)
%     if Tetrode(No)<10
%         EEGname=[params.ident '_csc_tt0' num2str(Tetrode(No))];
%     else
     if contains(iident,'cpp') || contains(iident,'sal')
         if Tetrode(No)<10
             EEGname=[iident(1:6) 'tt0' num2str(Tetrode(No)) iident(6:end)];
         else
             EEGname=[iident(1:6) 'tt' num2str(Tetrode(No)) iident(6:end)];
         end 
         load([EEGname],'csctimes','cscsamples','Fs','-mat')         

         if exist('Times','var')
             Fs = 3255;
             csctimes = Times; clear Times
             cscsamples = LFP_Samples; clear LFP_Samples
             save(EEGname,'Fs','csctimes','cscsamples')
         end
     elseif contains(iident,'Janni_2010408_Run') && str2double(iident(end-4))>1
        EEGname=['Janni_tt' num2str(Tetrode(No)) iident(6:end-4)]; 
         load([EEGname '.mat'],'csctimes','cscsamples','Fs')
     else
        EEGname=[iident(1:end-4) '_csc_tt' num2str(Tetrode(No))];    
        if ~isfile([EEGname '.mat'])
            dEEG = dir([iident(1:end-4) '*_csc_tt' num2str(Tetrode(No)) '.mat']);
            if size(dEEG,1)>1
                csctimes_temp = []; cscsamples_temp = [];
                for de = 1:size(dEEG,1)
                   load(dEEG(de).name,'csctimes','cscsamples','Fs','-mat')   
                   csctimes_temp = cat(1,csctimes_temp,csctimes);
                   cscsamples_temp = cat(1,cscsamples_temp,cscsamples);
                end
                csctimes = csctimes_temp; cscsamples = cscsamples_temp;
            else
                 load(dEEG(1).name,'csctimes','cscsamples','Fs','-mat')  
            end
            save([EEGname '.mat'],'csctimes','cscsamples','Fs')
        else
            load([EEGname '.mat'],'csctimes','cscsamples','Fs')
        end
     end
%     end
    
    LFP_Frequency=Fs;
    LFP=cscsamples;
    LFPTime=csctimes;        

    [b,a]=butter(3,[25/LFP_Frequency*2 50/LFP_Frequency*2],'bandpass');
    LFP=filtfilt(b,a,LFP);
    [~,locs]=findpeaks(-LFP); % recording is the reverse LFP
    LFPPeakTime=LFPTime(locs); % so peak gamma is phase 0/360, same as Brad    

    Cell=unique(rawspikedata(rawspikedata(:,3)==Tetrode(No),2));
    t=find(ismember(Spike(:,1),Cell));        
    [~,I]=histc(Spike(t,2),LFPPeakTime);
    SpikeGamma(t(I~=0))=(Spike(t(I~=0),2)-LFPPeakTime(I(I~=0)))./(LFPPeakTime(I(I~=0)+1)-LFPPeakTime(I(I~=0)))*pi*2; % linear interpolation of phase
    SpikeGamma(t(I==0))=99;      

    clear EEGname LFP LFPTime Times LFP_Samples       

end

[~,I]=max(OutMatrix'+InMatrix');

Mat=OutMatrix'+InMatrix';
I2=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);

StepGamma=zeros(sum(diff(Index)*4-4),4);% gamma phase, step size - peak decoded position,step size - weighted pos, sequence number
count=0;
SeqGamma = NaN(size(CandSeq,1),2);
for SeqNo=1:size(CandSeq,1)
    
    t=find(II>=Index(SeqNo)*4+1 & II<=Index(SeqNo+1)*4-1 & SpikeGamma~=99);    
    if ~isempty(t)                            
        SeqGamma(SeqNo,1)=circ_mean(SpikeGamma(t));  
        SeqGamma(SeqNo,2)=circ_r(SpikeGamma(t));  
    else
        SeqGamma(SeqNo,1:2)=[NaN NaN];
    end
        
    for i=Index(SeqNo)*4+2:Index(SeqNo+1)*4-3 
        count=count+1;
        t=find(II>=i & II<=i+2 & SpikeGamma~=99);                      
        if ~isempty(t)                            
            StepGamma(count,1)=circ_mean(SpikeGamma(t));  
        else
            StepGamma(count,1)=NaN;
        end
        StepGamma(count,2)=abs(I(i)-I(i-1)); % note you will need to divided by length of track
        StepGamma(count,3)=abs(I2(i)-I2(i-1)); % note you will need to divided by length of track
        StepGamma(count,4)=SeqNo;
    end     
    
end

StepGamma(StepGamma(:,1)==pi,1)=-pi; % because cicr_mean give (-pi pi], I changed it to [-pi pi) so it is correct for me to use histc    
SeqGamma(SeqGamma(:,1)==pi,1) = pi;
% Spike=Orig_Spike;

clear N EstBin A B I I2 Mat count SeqNo Orig_Spike II i t CellNo Cell Number No Tetrode a b pks locs LFPPeakTime Cand LFP_Frequency

cd(dirs.spikedatadir)
save(iident,'StepGamma','SeqGamma','SpikeGamma','-append')






