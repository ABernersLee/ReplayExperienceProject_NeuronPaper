function [OutMatrixShuff,InMatrixShuff] =  decode_spikedensity_events_shuffle(iident,numshuff)
load(iident,'params','hp_cells','hpinterneurons','vel','dirdat','spikedata','linposcat')

%%%%%
%%%%%
%%%%% Make cellID shuffle posteriors. Used in Hovers_Moving_laps.m for
%%%%% supplemental figure 6
%%%%%
%%%%%

%dirdat 1 is up/out, 0 is down/in

Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
disp(['Starting ' iident])
Pos = linposcat;
Spike=[spikedata(:,2) spikedata(:,1) Pos(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))];    

load(iident,'CandSeq','InFR','OutFR')


% decode events


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
Number_Of_Bins=size(OutFR,2);

% shuffle
InMatrixShuff = NaN(TimeBins,Number_Of_Bins,numshuff);
OutMatrixShuff = InMatrixShuff;
for ishuff = 1:numshuff
    ss = randperm(length(Cell_Number));
    OutFRs = OutFR(ss,:);
    InFRs = InFR(ss,:);
    
    term2=zeros(Number_Of_Bins,TimeBins);
    term4=term2;

    if Number_Of_Bins>TimeBins    
        for TBin=1:TimeBins
            term2(:,TBin)=prod(OutFRs.^binspike(:,TBin),1);    
            term4(:,TBin)=prod(InFRs.^binspike(:,TBin),1);
        end
    else
        for PosBin=1:Number_Of_Bins   
            term2(PosBin,:)=prod((OutFRs(:,PosBin)*ones(1,TimeBins)).^binspike,1);    
            term4(PosBin,:)=prod((InFRs(:,PosBin)*ones(1,TimeBins)).^binspike,1);
        end
    end

    term3=exp(-EstBin*sum(OutFRs,1)')*ones(1,TimeBins);

    term5=exp(-EstBin*sum(InFRs,1)')*ones(1,TimeBins);

    Mat=term2.*term3+term4.*term5;

    OutMat1=term2.*term3;
    OutMatrix=(OutMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';

    InMat1=term4.*term5;
    InMatrix=(InMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';    
    
    clear term2 term3 term4 term5 OutMat1 InMat1
    
    nospks = sum(binspike>0)==0;
    OutMatrix(nospks,:) = NaN;
    InMatrix(nospks,:) = NaN;
    
    OutMatrixShuff(:,:,ishuff) = OutMatrix;
    InMatrixShuff(:,:,ishuff) = InMatrix;
    clear OutMatrix InMatrix
    if rem(ishuff,10)==0
        disp(['Done ' num2str(ishuff) ' shuffles'])
    end
end




% save(iident,'OutMatrixShuff','InMatrixShuff','-append')
