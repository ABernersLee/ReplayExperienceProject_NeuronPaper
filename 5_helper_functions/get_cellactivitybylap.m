function get_cellactivitybylap(iident)

%%%%%
%%%%%
%%%%% Gets out aspects of how cells are firing and participating in
%%%%% candidate events. This is used in the big figure in 
%%%%% Plot_Main_HoverJump_Figures.m but isn't really used in the paper. 
%%%%%
%%%%%

load(iident)

Cell_Number=size(OutFR,1);

CN = hp_cells(~ismember(hp_cells,hpinterneurons));

CandSpike=zeros(Cell_Number,size(CandSeq,1));    
for CellNo=1:Cell_Number
    c=histc(Spike(Spike(:,1)==CN(CellNo),2),sortrows([CandSeq(:,1);CandSeq(:,2)]));
    CandSpike(CellNo,:)=c(1:2:end);
end
CandSpike2 = CandSpike; CandSpike2(CandSpike==0) = NaN;

% MaxLap=ceil(length(MidTime)/2);
MaxLap=length(MidTime)-1;
[~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);    

SpikePerEvent=NaN*ones(Cell_Number,MaxLap,5);
ActRatio=NaN*ones(Cell_Number,MaxLap,5);
SpikePerAct=NaN*ones(Cell_Number,MaxLap,5);
ActRatioRate=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActRate=NaN*ones(Cell_Number,MaxLap,5);
FRPerCell=NaN*ones(Cell_Number,MaxLap,5);
FRall = SpikePerEvent;

for Quan=1:5 % I,II,III,IV quadrant, counter clockwise; and II+III+IV together

    for LapNo=1:MaxLap           

        if Quan==1
            t=find(I==LapNo & abs(CandCorr)>0.6 & CandDis/size(OutFR,2)<0.4);
        elseif Quan==2
            t=find(I==LapNo & abs(CandCorr)>0.6 & CandDis/size(OutFR,2)>=0.4);                
        elseif Quan==3
            t=find(I==LapNo & abs(CandCorr)<=0.6 & CandDis/size(OutFR,2)>=0.4);
        elseif Quan==4
            t=find(I==LapNo & abs(CandCorr)<=0.6 & CandDis/size(OutFR,2)<0.4);
        else
            t=find(I==LapNo & (abs(CandCorr)<=0.6 | CandDis/size(OutFR,2)>=0.4) );
        end

        if length(t)>1
%             SpikePerEvent(:,LapNo,Quan)=sum(CandSpike(:,t)')/length(t);  
%             FRPerCell(:,LapNo,Quan)=sum(CandSpike(:,t)')/sum(CandSeq(t,2)-CandSeq(t,1)); 
%             ActRatio(:,LapNo,Quan)=sum(CandSpike(:,t)'>0)/length(t);
%             ActRatioRate(:,LapNo,Quan)=(sum(CandSpike(:,t)'>0))/sum(CandSeq(t,2)-CandSeq(t,1)); 
%             SpikePerAct(:,LapNo,Quan)=sum(CandSpike(:,t)')./sum(CandSpike(:,t)'>0);
%             SpikePerActRate(:,LapNo,Quan)=(sum(CandSpike(:,t)')/sum(CandSeq(t,2)-CandSeq(t,1)))./sum(CandSpike(:,t)'>0);
            SpikePerEvent(:,LapNo,Quan)=nanmean(CandSpike(:,t),2);  
            FRPerCell(:,LapNo,Quan)=mean(CandSpike(:,t)./(CandSeq(t,2)-CandSeq(t,1))',2);
            ActRatio(:,LapNo,Quan)=nanmean(~isnan(CandSpike2(:,t)),2);            
            ActRatioRate(:,LapNo,Quan)=mean((double(CandSpike(:,t)>0)./(CandSeq(t,2)-CandSeq(t,1))'),2)';
            SpikePerAct(:,LapNo,Quan)= nanmean(CandSpike2(:,t),2);
            SpikePerActRate(:,LapNo,Quan)=nanmean(CandSpike2(:,t)./(CandSeq(t,2)-CandSeq(t,1))',2);             
            FRall(:,LapNo,Quan)=sum(CandSpike(:,t)')/sum(CandSeq(t,2)-CandSeq(t,1)); 
        elseif length(t)==1
            SpikePerEvent(:,LapNo,Quan)=CandSpike(:,t); 
            FRPerCell(:,LapNo,Quan)=CandSpike(:,t)/(CandSeq(t,2)-CandSeq(t,1));
            ActRatio(:,LapNo,Quan)=CandSpike(:,t)>0;
            FRall(:,LapNo,Quan)=CandSpike(:,t)/sum(CandSeq(t,2)-CandSeq(t,1)); 
%             ActRatioRate(:,LapNo,Quan)=CandSpike(:,t)>0/(CandSeq(t,2)-CandSeq(t,1)); 
%             SpikePerAct(:,LapNo,Quan)=CandSpike(:,t)./(CandSpike(:,t)>0);
%             SpikePerActRate(:,LapNo,Quan)=(CandSpike(:,t)/(CandSeq(t,2)-CandSeq(t,1)))./(CandSpike(:,t)>0); %
            ActRatioRate(:,LapNo,Quan)=double(CandSpike(:,t)>0)./(CandSeq(t,2)-CandSeq(t,1)); 
            SpikePerAct(:,LapNo,Quan)= CandSpike2(:,t);
            SpikePerActRate(:,LapNo,Quan)=CandSpike2(:,t)./(CandSeq(t,2)-CandSeq(t,1))';
             
        end 
    end
end     

clear MaxLap A I Quan LapNo t Cell_Number CandSpike CellNo c


save(iident,'SpikePerEvent','ActRatio','SpikePerAct','ActRatioRate','SpikePerActRate','FRPerCell','FRall','-append')

    
    