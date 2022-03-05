function decode_spikedensity_events_lapbylap(iident,downsample)
load(iident,'MidTime','params','hp_cells','hpinterneurons','pos','vel','dirdat','spikedata','cm_conv','linposcat')

%%%%%
%%%%% Decode with lap by lap fields, overlapping bins. Used in
%%%%% RunReplay_ABL_lapbylap.m amd RunReplay_ABL_ymaze.m.
%%%%%

%dirdat 1 is up/out, 0 is down/in
Bin_Size=2.5;
Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
disp(['Number of cells: ' num2str(length(Cell_Number))])
VelThresh=5;


% %was used to get Midpoint
numlaps = length(MidTime)-1;
[~,I4]=histc(pos(:,1),[MidTime(1:end)]);
Pos = linposcat;
% [pks,locs]=findpeaks(-abs(Pos-max(Pos)/2-min(Pos)/2),'MINPEAKHEIGHT',-3,'MINPEAKDISTANCE',60);

cutoff1 = range(Pos)*(9/10)+min(Pos); %top
cutoff2 = (range(Pos)*(1/10))+min(Pos); % bottom
LapTime = min(pos(:,1));
Ilap = false(size(pos,1),numlaps);
I = Ilap;
Iup = NaN(size(pos,1),numlaps);
for ilap = 1:numlaps
      
    if sum(Pos(I4==ilap)<=cutoff2)>0
        %rat just ran down to the bottom on this lap (down)
        endlap = mean(pos(I4==ilap & Pos<=cutoff2));
        Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
        Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 0;
    elseif sum(Pos(I4==ilap)>=cutoff1)>0
        %rat just ran up to the top on this lap (up)
        endlap = mean(pos(I4==ilap & Pos>=cutoff1));
        Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
        Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 1;
    elseif sum(Pos(I4==ilap)>=cutoff1)==0 && sum(Pos(I4==ilap)<=cutoff2)==0
        cutoff11 = range(Pos)*(8/10)+min(Pos); %top
        cutoff22 = (range(Pos)*(2/10))+min(Pos); % bottom
         if sum(Pos(I4==ilap)<=cutoff22)>0
            %rat just ran down to the bottom on this lap (down)
            endlap = mean(pos(I4==ilap & Pos<=cutoff22));
            Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
            Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 0;
         elseif sum(Pos(I4==ilap)>=cutoff11)>0
            %rat just ran up to the top on this lap (up)
            endlap = mean(pos(I4==ilap & Pos>=cutoff11));
            Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
            Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 1;
         elseif sum(Pos(I4==ilap)<=cutoff22)==0 && sum(Pos(I4==ilap)>=cutoff11)==0
             badlap = 1; n1 = 8; n2 = 2;
             while badlap == 1
                 n1 = n1-1; n2 = n2+1;
                 cutoff11 = range(Pos)*(n1/10)+min(Pos); %top
                 cutoff22 = (range(Pos)*(n2/10))+min(Pos); % bottom
                  if sum(Pos(I4==ilap)<=cutoff22)>0
                    %rat just ran down to the bottom on this lap (down)
                    endlap = mean(pos(I4==ilap & Pos<=cutoff22));
                    Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
                    Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 0;
                    badlap = false;
                 elseif sum(Pos(I4==ilap)>=cutoff11)>0
                    %rat just ran up to the top on this lap (up)
                    endlap = mean(pos(I4==ilap & Pos>=cutoff11));
                    Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
                    Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 1;
                    badlap = false;
                  end
                  disp(['Bad Lap ' num2str(ilap)])
             end
             
         end        
    end
    LapTime = cat(1,LapTime,endlap);
    if ilap==1
        I(:,ilap) = Ilap(:,ilap);
    else
        I(:,ilap) = Ilap(:,ilap) | Ilap(:,ilap-1);
    end
end
   
% % figure; plot(pos(:,1),linposcat,'k'); hold on; plot(MidTime(:),200,'r*')
% figure; hold on; for ilap = 1:numlaps;  plot(pos(Ilap(:,ilap),1),Pos(Ilap(:,ilap),1)); end
% for ilap = 1:numlaps;  figure; hold on; plot(pos(:,1),Pos,'k'); plot(pos(I(:,ilap),1),Pos(I(:,ilap),1)); end

if params.Track_Type==1 
        
    
%make directional place fields - linear
    
    
    Spike=[spikedata(:,2) spikedata(:,1) Pos(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))];
    Spikelap = I(spikedata(:,3),:);

    Number_Of_Bins=ceil(range(Pos)/Bin_Size);
    OutFR=NaN(length(Cell_Number),Number_Of_Bins,numlaps);
    InFR=OutFR;
    InFRint = NaN(length(hpinterneurons),Number_Of_Bins,numlaps); 
    OutFRint = InFRint;
    Filter=fspecial('gaussian',[1 20],2); % ref Davidson et.al. 2006
        
    for ilap = 1:numlaps
        
         for i=1:length(Cell_Number)                     
            t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==1 & Spikelap(:,ilap)==1);
            if length(t)~=1
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==1 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            else
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh  & dirdat==1 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            end
            FR(isnan(FR) | isinf(FR))=0;
            FR=FR(1:Number_Of_Bins);
            OutFR(i,:,ilap)=filter2(Filter,FR');  

            t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==0 & Spikelap(:,ilap)==1);
            if length(t)~=1
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==0 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            else
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh & dirdat==0 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            end        
            FR(isnan(FR) | isinf(FR))=0;
            FR=FR(1:Number_Of_Bins);
            InFR(i,:,ilap)=filter2(Filter,FR');
         end
    
         %ABL added 8/23/2019
         for i=1:length(hpinterneurons)                     
            t=find(Spike(:,4)>VelThresh & Spike(:,1)==hpinterneurons(i) & Spike(:,5)==1 & Spikelap(:,ilap)==1);
            if length(t)~=1
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==1 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            else
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh  & dirdat==1 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            end
            FR(isnan(FR) | isinf(FR))=0;
            FR=FR(1:Number_Of_Bins);
            OutFRint(i,:,ilap)=filter2(Filter,FR');  

            t=find(Spike(:,4)>VelThresh & Spike(:,1)==hpinterneurons(i) & Spike(:,5)==0 & Spikelap(:,ilap)==1);
            if length(t)~=1
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==0 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            else
                FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh & dirdat==0 & I(:,ilap)==1),min(Pos):Bin_Size:max(Pos));   
            end        
            FR(isnan(FR) | isinf(FR))=0;
            FR=FR(1:Number_Of_Bins);
            InFRint(i,:,ilap)=filter2(Filter,FR');
         end
         
    end
    
     
    
    
    
    
elseif params.Track_Type==2
%     error('Shouldn''nt be running this session!')
    %make directional place fields - open
    Pos1 = pos(:,2:3)*cm_conv;    
    
    minsave = min(Pos1);
    Pos2 = ceil((Pos1-minsave+.001)/Bin_Size);
    clear Pos1
    nb=max(Pos2);
    OutFR=NaN(length(Cell_Number),nb(1),nb(2));
    InFR=OutFR;
    InFRint = NaN(length(hpinterneurons),nb(1),nb(2));
    OutFRint = InFRint;
    
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
    
     for i=1:length(hpinterneurons)
        t=find(Spike(:,4)>VelThresh & Spike(:,1)==hpinterneurons(i) & Spike(:,5)==1);
        if length(t)~=1
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)./histc(PosLong(vel>VelThresh & dirdat==1),1:Number_Of_Bins);   
        else
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)'./histc(PosLong(vel>VelThresh & dirdat==1),1:Number_Of_Bins);   
        end
        FR(isnan(FR) | isinf(FR))=0;        
        twod = reshape(FR,[nb(1) nb(2)]);
        OutFRint(i,:,:)=filter2(Filter,twod);  

        t=find(Spike(:,4)>VelThresh & Spike(:,1)==hpinterneurons(i) & Spike(:,5)==0);
        if length(t)~=1
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)./histc(PosLong(vel>VelThresh & dirdat==0),1:Number_Of_Bins);   
        else
            FR=60*histc(Spike(t,3),1:Number_Of_Bins)'./histc(PosLong(vel>VelThresh & dirdat==0),1:Number_Of_Bins);   
        end    
        FR(isnan(FR) | isinf(FR))=0;        
        twod = reshape(FR,[nb(1) nb(2)]);
        InFRint(i,:,:)=filter2(Filter,twod);  
        clear twod FR t 
    end
    clear PosLong
end

InFRlapsINT = InFRint;
OutFRlapsINT = OutFRint;


OutFRlaps = OutFR; %temp 
InFRlaps = InFR; %temp

clear OutFRint InFRint
save(iident,'OutFRlapsINT','InFRlapsINT','InFRlaps','OutFRlaps','-append')

if 1 %temp

clear Number_Of_Bins i Filter  dirdat hp_cells hpinterneurons

% get spike density events

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

%Actually doing by stopping period, look back two laps
% [~,I2]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(2:2:end);max(Spike(:,2))]); %starting with the second stopping period
% excl = I2 == 0 | I2>numlaps;
% CandSeq(excl,:) = [];
% disp(['Number of Candidate Events excluded: ' num2str(sum(excl))])
% I2(excl) = [];

[~,I2]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime;max(Spike(:,2))]);
% decode events
% [~,III]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[pos(:,1)]);
% figure; hold on; plot(pos(:,1),linposcat(:,1),'k'); plot(pos(III,1),linposcat(III,1),'r*')

if sum(I2 == 0)>0
    disp([ num2str(sum(I2 == 0)) ' Events before starting running'])
    if params.Novel==1
        disp(['Novel: ' num2str(sum(I2==0))])
%         I2(I2==0) = [];
%         CandSeq(I2==0,:) = [];
    end    
    I2(I2 == 0) = 1;
end
if sum(I2>numlaps)>0
    I2(I2>numlaps) = numlaps;    
end
EstBin=0.02;

Cand=CandSeq;
N=round(diff(CandSeq')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
clear N t

[A,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);
Index=round(B/(EstBin));


Orig_Spike=Spike;
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); 
end
Spike(~mod(I,2),:)=[];
clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin)*4;


% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for CellID=1:length(Cell_Number)    
    for i=1:4
        c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
%         c=histcounts(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
        binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
    end                      
end

InFRsave = InFR;
OutFRsave = OutFR;

clear Start_Time End_Time CellID c i 
if params.Track_Type == 1
%     excl = max([InFR OutFR],[],2)<1;
for ilap = 1:numlaps
    InFR(sum(InFR(:,:,ilap)~=0,2)==0,:,ilap) = 1;
    OutFR(sum(OutFR(:,:,ilap)~=0,2)==0,:,ilap) = 1;
end
    
elseif params.Track_Type == 2
%     error('Not supposed to use this session rn')
    InFR(sum(sum(InFR~=0,2),3)==0,:,:) = 1;
    OutFR(sum(sum(OutFR~=0,2),3)==0,:,:) = 1;
end
InFR(InFR==0) = .0001;
OutFR(OutFR==0) = .0001;

% find(sum(binspike,2)==0)
if params.Track_Type == 1

    Number_Of_Bins=size(OutFR,2);
    InMatrix = NaN(TimeBins,Number_Of_Bins); OutMatrix = InMatrix;
    for ilap = 1:numlaps
        if sum(I2==ilap)==0
            continue
        end
        timeindex = 4*Index(find(I2==ilap,1,'first'))+1:4*Index(find(I2==ilap,1,'last')+1);
        binspike1 = binspike(:,timeindex);
        OutFR1 = OutFR(:,:,ilap); InFR1 = InFR(:,:,ilap);
        TimeBins1 = size(binspike1,2);
        term2=zeros(Number_Of_Bins,TimeBins1);
        term4=term2;

        if Number_Of_Bins>TimeBins1    
            for TBin=1:TimeBins1
                term2(:,TBin)=prod(OutFR1.^binspike1(:,TBin),1);    
                term4(:,TBin)=prod(InFR1.^binspike1(:,TBin),1);
            end
        else
            for PosBin=1:Number_Of_Bins   
                term2(PosBin,:)=prod((OutFR1(:,PosBin)*ones(1,TimeBins1)).^binspike1,1);    
                term4(PosBin,:)=prod((InFR1(:,PosBin)*ones(1,TimeBins1)).^binspike1,1);
            end
        end

        term3=exp(-EstBin*sum(OutFR1,1)')*ones(1,TimeBins1);

        term5=exp(-EstBin*sum(InFR1,1)')*ones(1,TimeBins1);
    
        Mat=term2.*term3+term4.*term5;

        OutMat1=term2.*term3;
        OutMatrix(timeindex,:)=(OutMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';

        InMat1=term4.*term5;
        InMatrix(timeindex,:)=(InMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';    

        minsave = NaN;
        clear term2 term3 term4 term5 OutMat1 InMat1
    end
    
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
OutFRlaps = OutFRsave;
InFRlaps = InFRsave;
OutFR = nanmean(OutFRlaps,3);
InFR = nanmean(InFRlaps,3);

save(iident,'replayparams','CandSeq','OutFRlapsINT','InFRlapsINT','OutFR','InFR','OutFRlaps','InFRlaps','Spike','Index','Cell_Number','-append')

sz = whos('InMatrix');
if sz.bytes>1.7*10   
    OutMatrixSAVE2 = OutMatrix;
    InMatrixSAVE2 = InMatrix;
    clearvars -except iident OutMatrixSAVE2 InMatrixSAVE2
    wh = whos(matfile(iident));
    load(iident,wh.name)
    OutMatrix = OutMatrixSAVE2; InMatrix = InMatrixSAVE2;
    clear OutMatrixSAVE2 InMatrixSAVE2
    save(iident,'-v7.3')
%     save(iident,'InMatrix','OutMatrix','-append','-v7.3')
else
    save(iident,'OutMatrix','InMatrix','-append')
end
end
