function get_gammaprop(iident,dirs)

%%%%%
%%%%% Gets out aspects of gamma (freq, power, zscored power)
%%%%%

cd(dirs.spikedatadir)
load(iident,'params','CandSeq','OutMatrix','InMatrix')
cd(dirs.cscdatadir)

if contains(iident,'cpp') || contains(iident,'sal')
    EEGList=dir([iident(1:6) 'tt*' iident(6:end) '*']);
else         
    EEGList=dir([params.ident '_csc_tt*']);    
end

% this is just for the days when data aquired in different channel have different sizes
% so basically even I want to average everyting together, the data size
% is not exatcly the same so I have to just use the subset; so here I
% am using the sizes that are most common (even though the size difference were only in a few earlier animals)
sizes=zeros(size(EEGList,1),1);
for EEGno=1:size(EEGList,1)
    load(EEGList(EEGno).name,'csctimes','-mat');
    sizes(EEGno)=length(csctimes);       
end
SampleSize=[unique(sizes) zeros(length(unique(sizes)),1)];
for i=1:size(SampleSize,1)
    SampleSize(i,2)=length(find(sizes==SampleSize(i,1)));        
end
[A,I]=max(SampleSize(:,2));
t=find(sizes==SampleSize(I(1),1));
EEGList=EEGList(t);

LFP=zeros(SampleSize(I(1),1),1);

for EEGno=1:size(EEGList,1)
    if EEGno==1
        load(EEGList(EEGno).name,'csctimes','cscsamples','Fs','-mat')
        LFPTime=csctimes;
    else
        load(EEGList(EEGno).name,'cscsamples','Fs','-mat')
    end        
    LFP=LFP+cscsamples; 
end  
LFP=LFP/size(EEGList,1);

clear EEGList t A I sizes SampleSize i EEGno cscsamples Times

LFP_Frequency= Fs;
Filter=fspecial('gaussian',[1662 1],277); % brad's paper, gaussion SD=85ms

LFPsave = LFP;
[b,a]=butter(3,[25/LFP_Frequency*2 50/LFP_Frequency*2],'bandpass');
LFP=filtfilt(b,a,LFP);
[~,locs]=findpeaks(-LFP); % recording is the reverse LFP
LFPPeakTime=LFPTime(locs); % so peak gamma is phase 0/360, same as Brad     
GammaPower=filtfilt(Filter,1,abs(hilbert(LFP)));
GammaPowerZ = zscore(GammaPower);

[b,a]=butter(3,[60/LFP_Frequency*2 100/LFP_Frequency*2],'bandpass');
LFPHigh=filtfilt(b,a,LFPsave);
[~,locs]=findpeaks(-LFPHigh); % recording is the reverse LFP
LFPPeakTimeHigh=LFPTime(locs); % so peak gamma is phase 0/360, same as Brad     
GammaPowerHigh=filtfilt(Filter,1,abs(hilbert(LFPHigh)));    
GammaPowerHighZ = zscore(GammaPowerHigh);
[~,I]=histc(LFPPeakTime,sortrows([CandSeq(:,1);CandSeq(:,2)]));
[~,IH]=histc(LFPPeakTimeHigh,sortrows([CandSeq(:,1);CandSeq(:,2)]));
[~,I2]=histc(LFPTime,sortrows([CandSeq(:,1);CandSeq(:,2)]));

GammaProp=NaN*ones(size(CandSeq,1),3); % frequency, power, zscored power
GammaPropHigh = GammaProp;
% PolySlope=NaN*ones(size(CandSeq,1),2); % poly fit of peak decoded positions
for SeqNo=1:size(CandSeq,1)
    s=find(I==2*SeqNo-1,1,'first');
    e=find(I==2*SeqNo-1,1,'last');
    if ~isempty(s) && ~isempty(e) && s>1 && e<length(LFPPeakTime)              
        GammaProp(SeqNo,1)=(e-s)/(LFPPeakTime(e)-LFPPeakTime(s));  
    else
%         disp('lfp times dont match')
    end

    t=find(I2==2*SeqNo-1);
    if ~isempty(t)
        GammaProp(SeqNo,2)=mean(GammaPower(t));
        GammaProp(SeqNo,3)=mean(GammaPowerZ(t));
    end

    s=find(IH==2*SeqNo-1,1,'first');
    e=find(IH==2*SeqNo-1,1,'last');
    if ~isempty(s) && ~isempty(e) && s>1 && e<length(LFPPeakTime)              
        GammaPropHigh(SeqNo,1)=(e-s)/(LFPPeakTime(e)-LFPPeakTime(s));  
    else
        error('lfp times dont match')
    end

    t=find(I2==2*SeqNo-1);
    if ~isempty(t)
        GammaPropHigh(SeqNo,2)=mean(GammaPowerHigh(t));  
        GammaPropHigh(SeqNo,3)=mean(GammaPowerHighZ(t));  
    end

%     % peak decoded positions
%     Mat=OutMatrix(Index(SeqNo)+1:Index(SeqNo+1),:)+InMatrix(Index(SeqNo)+1:Index(SeqNo+1),:);
%     [~,II]=max(Mat');        
%     p=polyfit(1:length(II),II,1); % p(1) is the slope
%     PolySlope(SeqNo,1)=p(1);
% 
%     % weighted decoded positions
%     Mat=Mat';
%     II=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
%     p=polyfit(1:length(II),II,1); % p(1) is the slope
%     PolySlope(SeqNo,2)=p(1);        

end

clear LFPTime Filter LFP_Frequency a b LFP pks locs LFPPeakTime GammaPower A II I B I2 SeqNo s e t Mat II p 

cd(dirs.spikedatadir)
save(iident,'GammaProp','GammaPropHigh','-append')

