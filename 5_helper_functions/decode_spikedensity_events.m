function decode_spikedensity_events(iident,downsample)
load(iident,'params','hp_cells','hpinterneurons','pos','vel','dirdat','spikedata','cm_conv','linposcat')


%dirdat 1 is up/out, 0 is down/in
Bin_Size=2.5;
Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
disp(['Number of cells: ' num2str(length(Cell_Number))])
VelThresh=5;


if params.Track_Type==1 
        
    
%make directional place fields - linear
    Pos = linposcat;
    Spike=[spikedata(:,2) spikedata(:,1) Pos(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))];

    Number_Of_Bins=ceil(range(Pos)/Bin_Size);
    OutFR=zeros(length(Cell_Number),Number_Of_Bins);
    InFR=zeros(length(Cell_Number),Number_Of_Bins);
    
    Filter=fspecial('gaussian',[1 20],2); % ref Davidson et.al. 2006
     for i=1:length(Cell_Number)
        t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==1);
        if length(t)~=1
            FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==1),min(Pos):Bin_Size:max(Pos));   
        else
            FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh  & dirdat==1),min(Pos):Bin_Size:max(Pos));   
        end
        FR(isnan(FR) | isinf(FR))=0;
        FR=FR(1:Number_Of_Bins);
        OutFR(i,:)=filter2(Filter,FR');  

        t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==0);
        if length(t)~=1
            FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==0),min(Pos):Bin_Size:max(Pos));   
        else
            FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh & dirdat==0),min(Pos):Bin_Size:max(Pos));   
        end        
        FR(isnan(FR) | isinf(FR))=0;
        FR=FR(1:Number_Of_Bins);
        InFR(i,:)=filter2(Filter,FR');
    end
    clear Pos
elseif params.Track_Type==2
    %make directional place fields - open
    Pos1 = pos(:,2:3)*cm_conv;    
    
    minsave = min(Pos1);
    Pos2 = ceil((Pos1-minsave+.001)/Bin_Size);
    clear Pos1
    nb=max(Pos2);
    OutFR=zeros(length(Cell_Number),nb(1),nb(2));
    InFR=zeros(length(Cell_Number),nb(1),nb(2));
    Number_Of_Bins = nb(1)*nb(2);
    PosLong = sub2ind(nb,Pos2(:,1),Pos2(:,2));
    clear Pos2
    Spike=[spikedata(:,2) spikedata(:,1) PosLong(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))]; 

    Filter=fspecial('gaussian',[20 20],2); 
     for i=1:length(Cell_Number)
        t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==1);
        if length(t)~=1
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)./histc(PosLong(vel>VelThresh & dirdat==1),1:Number_Of_Bins);   
        else
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)'./histc(PosLong(vel>VelThresh & dirdat==1),1:Number_Of_Bins);   
        end
        FR(isnan(FR) | isinf(FR))=0;        
        twod = reshape(FR,[nb(1) nb(2)]);
        OutFR(i,:,:)=filter2(Filter,twod);  

        t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==0);
        if length(t)~=1
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)./histc(PosLong(vel>VelThresh & dirdat==0),1:Number_Of_Bins);   
        else
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)'./histc(PosLong(vel>VelThresh & dirdat==0),1:Number_Of_Bins);   
        end    
        FR(isnan(FR) | isinf(FR))=0;        
        twod = reshape(FR,[nb(1) nb(2)]);
        InFR(i,:,:)=filter2(Filter,twod);  
        clear twod FR t 
    end
    clear PosLong
end

clear Number_Of_Bins i Filter  dirdat hp_cells hpinterneurons

% get spike density events

% FOR DOWNSAMPLING ANALYSIS
if ~downsample
    % define candidate events
    Time_Bin_Size=0.001;
    Time_Bins=min(Spike(:,2)):Time_Bin_Size:min(Spike(:,2))+Time_Bin_Size*ceil(range(Spike(:,2))/Time_Bin_Size);
    [A,~]=histc(Spike(:,2),Time_Bins');
    Spike_Density=zeros(length(Time_Bins)-1,3);
    Spike_Density(:,2)=Time_Bins(1:end-1)'+Time_Bin_Size/2;
    Filter=fspecial('gaussian',[100 1],10); % ref Brad Pfeiffer 2013
    Spike_Density(:,1)=filter2(Filter,A(1:end-1));
    [~,I]=histc(Spike_Density(:,2),[min([pos(1,1) Time_Bins(1)]); pos(1:end-1,1)+diff(pos(:,1))/2 ; max([pos(end,1) Time_Bins(end)])]);
    Spike_Density(:,3)=vel(I);

    sMean=mean(Spike_Density(abs(Spike_Density(:,3))<VelThresh,1));
    sPeak=sMean+3*std(Spike_Density(abs(Spike_Density(:,3))<VelThresh,1));
    Check=diff([0 ; Spike_Density(:,1)>sMean & abs(Spike_Density(:,3))<VelThresh ; 0]);
    % Check=diff([0 ; Spike_Density(:,1)>sMean ; 0]);
    clear VelTresh

    DurationBound=[0.1 .5];
    CandSeq=[];
    SS=find(Check==1);
    EE=find(Check==-1);
    target=find(EE-SS>DurationBound(1)/Time_Bin_Size & EE-SS<DurationBound(2)/Time_Bin_Size);
    if ~isempty(target)
        CandSeq1=[Spike_Density(SS(target),2) Spike_Density(EE(target)-1,2)];
        delete=zeros(length(target),1);
        for i=1:length(target)
            if Spike_Density(SS(target(i)):EE(target(i))-1,1)<sPeak
                delete(i)=1;
            end    
        end
        target=find(delete==1);
        if ~isempty(target)
            CandSeq1(target,:)=[];
        end
    end
    CandSeq=[CandSeq;CandSeq1]; 
elseif downsample
    load(iident,'CandSeq')
end
clear A CandSeq1 Check delete DurationBound EE Filter i I sMean sPeak Spike_Density SS target Time_Bin_Size Time_Bins VelThresh spikedata vel
disp(['Number of Candidate Events: ' num2str(size(CandSeq,1))])

% decode events


EstBin=0.02;

Cand=CandSeq;
N=round(diff(CandSeq')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
clear N t

[~,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);

Orig_Spike=Spike;
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); % HERE
end
Spike(~mod(I,2),:)=[];
clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin)*4;
% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for CellID=1:length(Cell_Number)    
    for i=1:4
        c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
%         c=histcounts(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
        binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
    end                      
end

clear Start_Time End_Time CellID c i 
if params.Track_Type == 1
%     excl = max([InFR OutFR],[],2)<1;
    InFR(sum(InFR~=0,2)==0,:) = 1;
    OutFR(sum(OutFR~=0,2)==0,:) = 1;
    
elseif params.Track_Type == 2
    InFR(sum(sum(InFR~=0,2),3)==0,:,:) = 1;
    OutFR(sum(sum(OutFR~=0,2),3)==0,:,:) = 1;
end
InFR(InFR==0) = .0001;
OutFR(OutFR==0) = .0001;

% find(sum(binspike,2)==0)
if params.Track_Type == 1

    Number_Of_Bins=size(OutFR,2);
    term2=zeros(Number_Of_Bins,TimeBins);
    term4=term2;

    if Number_Of_Bins>TimeBins    
        for TBin=1:TimeBins
            term2(:,TBin)=prod(OutFR.^binspike(:,TBin),1);    
            term4(:,TBin)=prod(InFR.^binspike(:,TBin),1);
        end
    else
        for PosBin=1:Number_Of_Bins   
            term2(PosBin,:)=prod((OutFR(:,PosBin)*ones(1,TimeBins)).^binspike,1);    
            term4(PosBin,:)=prod((InFR(:,PosBin)*ones(1,TimeBins)).^binspike,1);
        end
    end

    term3=exp(-EstBin*sum(OutFR,1)')*ones(1,TimeBins);

    term5=exp(-EstBin*sum(InFR,1)')*ones(1,TimeBins);

    Mat=term2.*term3+term4.*term5;

    OutMat1=term2.*term3;
    OutMatrix=(OutMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';

    InMat1=term4.*term5;
    InMatrix=(InMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';    

    minsave = NaN;
    clear term2 term3 term4 term5 OutMat1 InMat1
    
elseif params.Track_Type==2

    Field_Data = permute(OutFR,[2 3 1]);
    OutMat=zeros(size(Field_Data,1),size(Field_Data,2),TimeBins);
    temp1=sum(Field_Data,3);
    temp(:,:)=temp1(:,:,1);

    for i=1:size(Field_Data,1)
      for j=1:size(Field_Data,2)
        CurrentFiringRate(:,1)=Field_Data(i,j,:);
        OutMat(i,j,:)=prod((CurrentFiringRate*ones(1,TimeBins)).^binspike,1).*temp(i,j);
      end
    end

    Field_Data = permute(InFR,[2 3 1]);
    InMat=zeros(size(Field_Data,1),size(Field_Data,2),TimeBins);
    temp1=sum(Field_Data,3);
    temp(:,:)=temp1(:,:,1);

    for i=1:size(Field_Data,1)
      for j=1:size(Field_Data,2)
        CurrentFiringRate(:,1)=Field_Data(i,j,:);
        InMat(i,j,:)=prod((CurrentFiringRate*ones(1,TimeBins)).^binspike,1).*temp(i,j);
      end
    end
    
    SumMat = InMat + OutMat;

    InMatrix = zeros(size(InMat));
    OutMatrix = zeros(size(OutMat));
    for i= 1:TimeBins
        InMatrix(:,:,i) = InMat(:,:,i)./sum(sum(SumMat(:,:,i)));
        OutMatrix(:,:,i) = OutMat(:,:,i)./sum(sum(SumMat(:,:,i)));
    end
        
    clear Field_Data CurrentFIringRate i j temp temp1
    

end

replayparams.Bin_Size = Bin_Size;
replayparams.minsave = minsave;
replayparams.cm_conv = cm_conv;
Spike=Orig_Spike;
Index=round(B/(EstBin));

save(iident,'replayparams','CandSeq','OutFR','InFR','Spike','Index','OutMatrix','InMatrix','Cell_Number','-append')
