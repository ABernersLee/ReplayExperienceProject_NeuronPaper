function decode_spikedensity_events_lapbylap_NonOverLapping(iident)

%%%%%
%%%%% Decode with lap by lap fields, non-overlapping bins. Used in
%%%%% RunReplay_ABL_lapbylap.m
%%%%%

load(iident,'MidTime','params','hp_cells','hpinterneurons','vel','dirdat','spikedata','linposcat')
Pos = linposcat;
%dirdat 1 is up/out, 0 is down/in
Spike=[spikedata(:,2) spikedata(:,1) Pos(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))];
Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
disp(['Number of cells: ' num2str(length(Cell_Number))])
numlaps = length(MidTime)-1;
% get spike density events
load(iident,'CandSeq','InFRlaps','OutFRlaps')
InFR = InFRlaps;
OutFR = OutFRlaps;
[~,I2]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime;max(Spike(:,2))]);

if sum(I2 == 0)>0
    disp([ num2str(sum(I2 == 0)) ' Events before starting running'])
    if params.Novel==1
        disp(['Novel: ' num2str(sum(I2==0))])
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

[~,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);
Index=round(B/(EstBin));

    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); 
end
Spike(~mod(I,2),:)=[];
clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin);


% Bayesian Decoding - 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for CellID=1:length(Cell_Number)
    binspike(CellID,:)=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time:EstBin:Start_Time+EstBin*(TimeBins-1));       
end


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


Number_Of_Bins=size(OutFR,2);
InMatrix = NaN(TimeBins,Number_Of_Bins); OutMatrix = InMatrix;
for ilap = 1:numlaps
    if sum(I2==ilap)==0
        continue
    end
    timeindex = Index(find(I2==ilap,1,'first'))+1:Index(find(I2==ilap,1,'last')+1);
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

    clear term2 term3 term4 term5 OutMat1 InMat1
end
    



D = InMatrix'+OutMatrix';

[CandCorrNOL,TmShuffleCorrNOL,CandDisNOL,TmShuffleDisNOL] = shuffle_wc_nonoverlap(D,Index,CandSeq);



save(iident,'CandCorrNOL','TmShuffleCorrNOL','CandDisNOL','TmShuffleDisNOL','-append')

