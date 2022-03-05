function get_cellactivitybylap_hoverjump(iident)


%%%%%
%%%%%
%%%%% Gets out aspects of how cells are firing and participating in
%%%%% hovers/moves in candidate events. This is used in the big figure in 
%%%%% Plot_Main_HoverJump_Figures.m but isn't really used in the paper. 
%%%%%
%%%%%

load(iident)

Cell_Number=size(OutFR,1);
% MaxLap=ceil(length(MidTime)/2);
MaxLap=length(MidTime)-1;

SpikePerEventHover=NaN*ones(Cell_Number,MaxLap,5);
ActRatioHover=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActHover=NaN*ones(Cell_Number,MaxLap,5);
ActRatioRateHover=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActRateHover=NaN*ones(Cell_Number,MaxLap,5);
FRPerCellHover=NaN*ones(Cell_Number,MaxLap,5);

SpikePerEventMove=NaN*ones(Cell_Number,MaxLap,5);
ActRatioMove=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActMove=NaN*ones(Cell_Number,MaxLap,5);
ActRatioRateMove=NaN*ones(Cell_Number,MaxLap,5);
SpikePerActRateMove=NaN*ones(Cell_Number,MaxLap,5);
FRPerCellMove=NaN*ones(Cell_Number,MaxLap,5);
CandDis2 = CandDis/size(OutFR,2);

if ~isempty(hovertime2)
clear hovertime movetime
hovertime(:,1) = CandSeq(hovertime2(:,3),1)+(hovertime2(:,1)-1)*.005;
hovertime(:,2) = CandSeq(hovertime2(:,3),1)+(hovertime2(:,2))*.005;
movetime(:,1) = CandSeq(movetime2(:,3),1)+(movetime2(:,1)-1)*.005;
movetime(:,2) = CandSeq(movetime2(:,3),1)+(movetime2(:,2))*.005;

CN = hp_cells(~ismember(hp_cells,hpinterneurons));

HoverSpike=zeros(Cell_Number,size(hovertime,1));    
MoveSpike=zeros(Cell_Number,size(movetime,1));    
for CellNo=1:Cell_Number
    c=histc(Spike(Spike(:,1)==CN(CellNo),2),sortrows([hovertime(:,1);hovertime(:,2)]));
    HoverSpike(CellNo,:)=c(1:2:end);
    c=histc(Spike(Spike(:,1)==CN(CellNo),2),sortrows([movetime(:,1);movetime(:,2)]));
    MoveSpike(CellNo,:)=c(1:2:end);
end
MoveSpike2 = MoveSpike; MoveSpike2(MoveSpike==0) = NaN;
HoverSpike2 = HoverSpike; HoverSpike2(HoverSpike==0) = NaN;

[~,Ihover]=histc(hovertime(:,1)/2+hovertime(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);   
[~,Imove]=histc(movetime(:,1)/2+movetime(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);   


for Quan=1:5 % I,II,III,IV quadrant, counter clockwise; and II+III+IV together

    for LapNo=1:MaxLap           

        if Quan==1
            t=find(Ihover==LapNo & abs(CandCorr(hovertime2(:,3)))>0.6 & CandDis2(hovertime2(:,3))<0.4);
        elseif Quan==2
            t=find(Ihover==LapNo & abs(CandCorr(hovertime2(:,3)))>0.6 & CandDis2(hovertime2(:,3))>=0.4);                
        elseif Quan==3
            t=find(Ihover==LapNo & abs(CandCorr(hovertime2(:,3)))<=0.6 & CandDis2(hovertime2(:,3))>=0.4);
        elseif Quan==4
            t=find(Ihover==LapNo & abs(CandCorr(hovertime2(:,3)))<=0.6 & CandDis2(hovertime2(:,3))<0.4);
        else
            t=find(Ihover==LapNo & (abs(CandCorr(hovertime2(:,3)))<=0.6 | CandDis2(hovertime2(:,3))>=0.4) );
        end

        if length(t)>1
            SpikePerEventHover(:,LapNo,Quan)=sum(HoverSpike(:,t)')/length(t); 
            ActRatioHover(:,LapNo,Quan)=sum(HoverSpike(:,t)'>0)/length(t);
            ActRatioRateHover(:,LapNo,Quan)=(sum(HoverSpike(:,t)'>0))/sum(hovertime(t,2)-hovertime(t,1));% (sum(HoverSpike(:,t)'>0)/length(t))/sum(hovertime(t,2)-hovertime(t,1));
            SpikePerActHover(:,LapNo,Quan)=sum(HoverSpike(:,t)')./sum(HoverSpike(:,t)'>0);
            SpikePerActRateHover(:,LapNo,Quan)=(sum(HoverSpike(:,t)')/sum(hovertime(t,2)-hovertime(t,1)))./sum(HoverSpike(:,t)'>0);
            FRPerCellHover(:,LapNo,Quan)=sum(HoverSpike(:,t)')/sum(hovertime(t,2)-hovertime(t,1)); 
            
            SpikePerEventHover(:,LapNo,Quan)=nanmean(HoverSpike(:,t),2);  
            FRPerCellHover(:,LapNo,Quan)=mean(HoverSpike(:,t)./(hovertime(t,2)-hovertime(t,1))',2);
            ActRatioHover(:,LapNo,Quan)=nanmean(~isnan(HoverSpike2(:,t)),2);            
            ActRatioRateHover(:,LapNo,Quan)=mean((double(HoverSpike(:,t)>0)./(hovertime(t,2)-hovertime(t,1))'),2)';
            SpikePerActHover(:,LapNo,Quan)= nanmean(HoverSpike2(:,t),2);
            SpikePerActRateHover(:,LapNo,Quan)=nanmean(HoverSpike2(:,t)./(hovertime(t,2)-hovertime(t,1))',2);             
        elseif length(t)==1
            SpikePerEventHover(:,LapNo,Quan)=HoverSpike(:,t); 
            FRPerCellHover(:,LapNo,Quan)=HoverSpike(:,t)/(hovertime(t,2)-hovertime(t,1)); 
            ActRatioHover(:,LapNo,Quan)=HoverSpike(:,t)>0;
            
            ActRatioRateHover(:,LapNo,Quan)=HoverSpike(:,t)>0/sum(hovertime(t,2)-hovertime(t,1));
            SpikePerActHover(:,LapNo,Quan)=HoverSpike(:,t)./(HoverSpike(:,t)>0);
            SpikePerActRateHover(:,LapNo,Quan)=(HoverSpike(:,t)/(hovertime(t,2)-hovertime(t,1)))./(HoverSpike(:,t)>0);
            
            ActRatioRateHover(:,LapNo,Quan)=double(HoverSpike(:,t)>0)./(hovertime(t,2)-hovertime(t,1)); 
            SpikePerActHover(:,LapNo,Quan)= HoverSpike2(:,t);
            SpikePerActRateHover(:,LapNo,Quan)=HoverSpike2(:,t)./(hovertime(t,2)-hovertime(t,1))';
        end 
        
        if Quan==1
            t=find(Imove==LapNo & abs(CandCorr(movetime2(:,3)))>0.6 & CandDis2(movetime2(:,3))<0.4);
        elseif Quan==2
            t=find(Imove==LapNo & abs(CandCorr(movetime2(:,3)))>0.6 & CandDis2(movetime2(:,3))>=0.4);                
        elseif Quan==3
            t=find(Imove==LapNo & abs(CandCorr(movetime2(:,3)))<=0.6 & CandDis2(movetime2(:,3))>=0.4);
        elseif Quan==4
            t=find(Imove==LapNo & abs(CandCorr(movetime2(:,3)))<=0.6 & CandDis2(movetime2(:,3))<0.4);
        else
            t=find(Imove==LapNo & (abs(CandCorr(movetime2(:,3)))<=0.6 | CandDis2(movetime2(:,3))>=0.4) );
        end

        if length(t)>1
%             SpikePerEventMove(:,LapNo,Quan)=sum(MoveSpike(:,t)')/length(t); 
%             ActRatioMove(:,LapNo,Quan)=sum(MoveSpike(:,t)'>0)/length(t);
%             ActRatioRateMove(:,LapNo,Quan)=(sum(MoveSpike(:,t)'>0))/sum(movetime(t,2)-movetime(t,1)); %(sum(MoveSpike(:,t)'>0)/length(t))/sum(movetime(t,2)-movetime(t,1));
%             SpikePerActMove(:,LapNo,Quan)=sum(MoveSpike(:,t)')./sum(MoveSpike(:,t)'>0);
%             SpikePerActRateMove(:,LapNo,Quan)=(sum(MoveSpike(:,t)')/sum(movetime(t,2)-movetime(t,1)))./sum(MoveSpike(:,t)'>0);
%             FRPerCellMove(:,LapNo,Quan)=sum(MoveSpike(:,t)')/sum(movetime(t,2)-movetime(t,1)); 
            
            SpikePerEventMove(:,LapNo,Quan)=nanmean(MoveSpike(:,t),2);  
            FRPerCellMove(:,LapNo,Quan)=mean(MoveSpike(:,t)./(movetime(t,2)-movetime(t,1))',2);
            ActRatioMove(:,LapNo,Quan)=nanmean(~isnan(MoveSpike2(:,t)),2);            
            ActRatioRateMove(:,LapNo,Quan)=mean((double(MoveSpike(:,t)>0)./(movetime(t,2)-movetime(t,1))'),2)';
            SpikePerActMove(:,LapNo,Quan)= nanmean(MoveSpike2(:,t),2);
            SpikePerActRateMove(:,LapNo,Quan)=nanmean(MoveSpike2(:,t)./(movetime(t,2)-movetime(t,1))',2);             
        elseif length(t)==1
            SpikePerEventMove(:,LapNo,Quan)=MoveSpike(:,t); 
            FRPerCellMove(:,LapNo,Quan)=MoveSpike(:,t)/(movetime(t,2)-movetime(t,1)); 
            ActRatioMove(:,LapNo,Quan)=MoveSpike(:,t)>0;
            
%             ActRatioRateMove(:,LapNo,Quan)=MoveSpike(:,t)>0/(movetime(t,2)-movetime(t,1));
%             SpikePerActMove(:,LapNo,Quan)=MoveSpike(:,t)./(MoveSpike(:,t)>0);
%             SpikePerActRateMove(:,LapNo,Quan)=(MoveSpike(:,t)/(movetime(t,2)-movetime(t,1)))./(MoveSpike(:,t)>0);
            
            ActRatioRateMove(:,LapNo,Quan)=double(MoveSpike(:,t)>0)./(movetime(t,2)-movetime(t,1)); 
            SpikePerActMove(:,LapNo,Quan)= MoveSpike2(:,t);
            SpikePerActRateMove(:,LapNo,Quan)=MoveSpike2(:,t)./(movetime(t,2)-movetime(t,1))';
        end 
        
        
    end
end     

clear MaxLap A I Quan LapNo t Cell_Number MoveSpike HoverSpike CellNo c
else
    disp('no hovertime2')
end

save(iident,'SpikePerEventHover','ActRatioHover','SpikePerActHover','ActRatioRateHover','SpikePerActRateHover','FRPerCellHover',...
    'SpikePerEventMove','ActRatioMove','SpikePerActMove','ActRatioRateMove','SpikePerActRateMove','FRPerCellMove','-append')

    
    