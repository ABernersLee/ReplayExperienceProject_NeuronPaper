function decode_spikedensity_events_OneLap_Epoch1(iident)
load(iident,'params','hp_cells','hpinterneurons','pos','vel','dirdat','spikedata','CandSeq','linposcat')

%%%%%
%%%%% Non-overlapping, decoding making the spikedensity sperately for each
%%%%% epoch using first run (run1) fields. Used in RunReplay_ABL_OneLap.m
%%%%%

%dirdat 1 is up/out, 0 is down/in
Bin_Size=2.5;
Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
disp(['Number of cells: ' num2str(length(Cell_Number))])
VelThresh=5;

Pos = linposcat;    
Pos(pos(:,1)<=params.Run_Times(1,1) | pos(:,1)>=params.Run_Times(1,2),:) = NaN;    
posbins = min(Pos):Bin_Size:max(Pos);

spikedata1 = spikedata(~isnan(spikedata(:,3)),:);
spikedata2 = spikedata(isnan(spikedata(:,3)),:);
SpikeSave = [[spikedata1(:,2) spikedata1(:,1) Pos(spikedata1(:,3)) vel(spikedata1(:,3)) dirdat(spikedata1(:,3))]; [spikedata2(:,2) spikedata2(:,1) NaN(size(spikedata2,1),3)]];
[~,ind] = sort(SpikeSave(:,2),'ascend');
SpikeSave = SpikeSave(ind,:);
clear spikedata 
Spike=[spikedata1(:,2) spikedata1(:,1) Pos(spikedata1(:,3)) vel(spikedata1(:,3)) dirdat(spikedata1(:,3))];
% Spike=[spikedata(:,2) spikedata(:,1) Pos(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))];
% SpikeSave = Spike;
Spike = Spike(Spike(:,2)>params.Run_Times(1,1) & Spike(:,2)>pos(1,1) & Spike(:,2)<params.Run_Times(1,2),:);
Number_Of_Bins=length(posbins);
FREpochOne=zeros(length(Cell_Number),Number_Of_Bins);    
Filter=fspecial('gaussian',[1 20],2); % changed from 4 to 2  11/21/2019 and again (because it was 4) 2/25/20
 for i=1:length(Cell_Number)
    t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i));

    if length(t)~=1
        FR=60*histc(Spike(t,3),posbins)./histc(Pos(vel>VelThresh),posbins);   
    else
        FR=60*histc(Spike(t,3),posbins)'./histc(Pos(vel>VelThresh),posbins);   
    end
    FR(isnan(FR) | isinf(FR))=0;
    FR=FR(1:Number_Of_Bins);
    FREpochOne(i,:)=filter2(Filter,FR');  
 end

clear Number_Of_Bins i Filter  dirdat hp_cells hpinterneurons
% Pos = linposcat;
Spike = SpikeSave;
% vel = velsave;


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

% Orig_Spike=Spike;
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); % HERE
end
Spike(~mod(I,2),:)=[];
clear Cand I i

Start_Time=0;
End_Time=B(end);
% TimeBins=round((End_Time-Start_Time)/EstBin)*4;
TimeBins=round((End_Time-Start_Time)/EstBin);
% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);

% for CellID=1:length(Cell_Number)    
% %     for i=1:4
%         c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
% %         c=histcounts(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
%         binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
%     end                      
% end

% Bayesian Decoding - 20ms window estimate, nonoverlapping
for CellID=1:length(Cell_Number)
    binspike(CellID,:)=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time:EstBin:Start_Time+EstBin*(TimeBins-1));       
end  

clear Start_Time End_Time CellID c i 

FREpochOne(sum(FREpochOne~=0,2)==0,:) = 1;
FREpochOne(FREpochOne==0) = .0001;




Number_Of_Bins=size(FREpochOne,2);
term2=zeros(Number_Of_Bins,TimeBins);    

if Number_Of_Bins>TimeBins    
    for TBin=1:TimeBins
        term2(:,TBin)=prod(FREpochOne.^binspike(:,TBin),1);                
    end
else
    for PosBin=1:Number_Of_Bins   
        term2(PosBin,:)=prod((FREpochOne(:,PosBin)*ones(1,TimeBins)).^binspike,1);                
    end
end

term3=exp(-EstBin*sum(FREpochOne,1)')*ones(1,TimeBins);


Mat=term2.*term3;

Mat1=term2.*term3;
MatrixEpoch1=(Mat1./(ones(size(Mat,1),1)*sum(Mat,1)))';
clear term2 term3 Mat1
    
save(iident,'FREpochOne','MatrixEpoch1','-append')
