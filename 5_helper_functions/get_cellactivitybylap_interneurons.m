function get_cellactivitybylap_interneurons(iident)


%%%%%
%%%%%
%%%%% Gets out aspects of how interneurons are firing and participating in
%%%%% candidate events. This is used in the big figure in 
%%%%% Plot_Main_HoverJump_Figures.m but isn't really used in the paper. 
%%%%%
%%%%%

load(iident)


CN = hpinterneurons;
Cell_Number = length(CN);

CandSpike=zeros(Cell_Number,size(CandSeq,1));    
for CellNo=1:Cell_Number
    c=histc(Spike(Spike(:,1)==CN(CellNo),2),sortrows([CandSeq(:,1);CandSeq(:,2)]));
    CandSpike(CellNo,:)=c(1:2:end);
end
CandSpike2 = CandSpike; CandSpike2(CandSpike==0) = NaN;
% MaxLap=ceil(length(MidTime)/2);
MaxLap=length(MidTime)-1;
[~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);    

SpikePerEventIN=NaN*ones(Cell_Number,MaxLap,5);
ActRatioIN=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActIN=NaN*ones(Cell_Number,MaxLap,5);
ActRatioRateIN=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActRateIN=NaN*ones(Cell_Number,MaxLap,5);
FRPerCellIN=NaN*ones(Cell_Number,MaxLap,5);
FRallIN = SpikePerEventIN;
if ~isempty(CN)
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
            SpikePerEventIN(:,LapNo,Quan)=nanmean(CandSpike(:,t),2);  
            FRPerCellIN(:,LapNo,Quan)=mean(CandSpike(:,t)./(CandSeq(t,2)-CandSeq(t,1))',2);
            ActRatioIN(:,LapNo,Quan)=nanmean(~isnan(CandSpike2(:,t)),2);            
            ActRatioRateIN(:,LapNo,Quan)=mean((double(CandSpike(:,t)>0)./(CandSeq(t,2)-CandSeq(t,1))'),2)';
            SpikePerActIN(:,LapNo,Quan)= nanmean(CandSpike2(:,t),2);
            SpikePerActRateIN(:,LapNo,Quan)=nanmean(CandSpike2(:,t)./(CandSeq(t,2)-CandSeq(t,1))',2);     
            
            FRallIN(:,LapNo,Quan)=sum(CandSpike(:,t)')/sum(CandSeq(t,2)-CandSeq(t,1)); 
        elseif length(t)==1
            SpikePerEventIN(:,LapNo,Quan)=CandSpike(:,t); 
            FRPerCellIN(:,LapNo,Quan)=CandSpike(:,t)/(CandSeq(t,2)-CandSeq(t,1)); 
            ActRatioIN(:,LapNo,Quan)=CandSpike(:,t)>0;
            
            ActRatioRateIN(:,LapNo,Quan)=double(CandSpike(:,t)>0)./(CandSeq(t,2)-CandSeq(t,1)); 
            SpikePerActIN(:,LapNo,Quan)= CandSpike2(:,t);
            SpikePerActRateIN(:,LapNo,Quan)=CandSpike2(:,t)./(CandSeq(t,2)-CandSeq(t,1))';
            
            FRallIN(:,LapNo,Quan)=CandSpike(:,t)/sum(CandSeq(t,2)-CandSeq(t,1)); 
            
        end 
    end
end     

clear MaxLap A I Quan LapNo t Cell_Number CandSpike CellNo c

end
save(iident,'SpikePerEventIN','ActRatioIN','SpikePerActIN','ActRatioRateIN','SpikePerActRateIN','FRPerCellIN','FRallIN','-append')

    
    