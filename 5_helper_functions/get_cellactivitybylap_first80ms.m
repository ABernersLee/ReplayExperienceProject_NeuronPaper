function get_cellactivitybylap_first80ms(iident)


%%%%%
%%%%%
%%%%% Gets out aspects of how cells are firing and participating in
%%%%% candidate events. I was specifically looking to see if anything about 
%%%%% the beggining of the event could predict how short/long the replay
%%%%% would become.
%%%%% This is used in the big figure in 
%%%%% Plot_Main_HoverJump_Figures.m but isn't used in the paper. 
%%%%%
%%%%%

load(iident)
addon = .08;
% excl = diff(CandSeq,[],2)<= addon;
CN = hp_cells(~ismember(hp_cells,hpinterneurons));
Cell_Number=size(OutFR,1);
CandSpike1=zeros(Cell_Number,size(CandSeq,1));    CandSpike2 = CandSpike1;
for CellNo=1:Cell_Number
    c=histc(Spike(Spike(:,1)==CN(CellNo),2),sortrows([CandSeq(:,1);CandSeq(:,1)+addon]));
    CandSpike1(CellNo,:)=c(1:2:end);
    c=histc(Spike(Spike(:,1)==CN(CellNo),2),sortrows([CandSeq(:,1)+addon;CandSeq(:,2)]));
    CandSpike2(CellNo,:)=c(1:2:end);
end

% MaxLap=ceil(length(MidTime)/2);
MaxLap=length(MidTime)-1;
[~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);    

SpikePerEventFirst80=NaN*ones(Cell_Number,MaxLap,5,2);
ActRatioFirst80=NaN*ones(Cell_Number,MaxLap,5,2);
SpikePerActFirst80=NaN*ones(Cell_Number,MaxLap,5,2);
ActRatioRateFirst80=NaN*ones(Cell_Number,MaxLap,5,2);
SpikePerActRateFirst80=NaN*ones(Cell_Number,MaxLap,5,2);
FRPerCellFirst80=NaN*ones(Cell_Number,MaxLap,5,2);
for iff = 1:2
    if iff == 1
        CandSpike = CandSpike1; CandSeq1 = [CandSeq(:,1) CandSeq(:,1)+addon]; 
    elseif iff ==2  
        CandSpike = CandSpike2; CandSeq1 = [CandSeq(:,1)+addon CandSeq(:,2)]; 
    end
    CandSpike2B = CandSpike; CandSpike2B(CandSpike==0) = NaN;
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
%                 SpikePerEventFirst80(:,LapNo,Quan,iff)=sum(CandSpike(:,t)')/length(t);             
                SpikePerEventFirst80(:,LapNo,Quan)=nanmean(CandSpike(:,t),2);  
%                 ActRatioFirst80(:,LapNo,Quan,iff)=sum(CandSpike(:,t)'>0)/length(t);
                ActRatioFirst80(:,LapNo,Quan)=nanmean(~isnan(CandSpike2B(:,t)),2);
%                 ActRatioRateFirst80(:,LapNo,Quan,iff)=(sum(CandSpike(:,t)'>0))/sum(CandSeq1(t,2)-CandSeq1(t,1)); %(sum(CandSpike(:,t)'>0)/length(t))/sum(CandSeq1(t,2)-CandSeq1(t,1));
                ActRatioRateFirst80(:,LapNo,Quan)=mean((double(CandSpike(:,t)>0)./(CandSeq(t,2)-CandSeq(t,1))'),2)';
%                 SpikePerActFirst80(:,LapNo,Quan,iff)=sum(CandSpike(:,t)')./sum(CandSpike(:,t)'>0);
                SpikePerActFirst80(:,LapNo,Quan)= nanmean(CandSpike2B(:,t),2);
%                 SpikePerActRateFirst80(:,LapNo,Quan,iff)=(sum(CandSpike(:,t)')/sum(CandSeq1(t,2)-CandSeq1(t,1)))./sum(CandSpike(:,t)'>0);
                SpikePerActRateFirst80(:,LapNo,Quan)=nanmean(CandSpike2B(:,t)./(CandSeq(t,2)-CandSeq(t,1))',2);
%                 FRPerCellFirst80(:,LapNo,Quan,iff)=sum(CandSpike(:,t)')/sum(CandSeq1(t,2)-CandSeq1(t,1)); 
                FRPerCellFirst80(:,LapNo,Quan)=mean(CandSpike(:,t)./(CandSeq(t,2)-CandSeq(t,1))',2);
            elseif length(t)==1
                SpikePerEventFirst80(:,LapNo,Quan,iff)=CandSpike(:,t); 
                ActRatioFirst80(:,LapNo,Quan,iff)=CandSpike(:,t)>0;
%                 ActRatioRateFirst80(:,LapNo,Quan,iff)=CandSpike(:,t)>0/(CandSeq1(t,2)-CandSeq1(t,1));
                ActRatioRateFirst80(:,LapNo,Quan)=double(CandSpike(:,t)>0)./(CandSeq(t,2)-CandSeq(t,1)); 
%                 SpikePerActFirst80(:,LapNo,Quan,iff)=CandSpike(:,t)./(CandSpike(:,t)>0);
                SpikePerActFirst80(:,LapNo,Quan)= CandSpike2B(:,t);
%                 SpikePerActRateFirst80(:,LapNo,Quan,iff)=(CandSpike(:,t)/(CandSeq1(t,2)-CandSeq1(t,1)))./(CandSpike(:,t)>0);
                SpikePerActRateFirst80(:,LapNo,Quan)=CandSpike2B(:,t)./(CandSeq(t,2)-CandSeq(t,1))';
                FRPerCellFirst80(:,LapNo,Quan,iff)=CandSpike(:,t)/(CandSeq1(t,2)-CandSeq1(t,1)); 
            end 
        end
    end     
end
clear MaxLap A I Quan LapNo t Cell_Number CandSpike CellNo c


save(iident,'SpikePerEventFirst80','ActRatioFirst80','SpikePerActFirst80','ActRatioRateFirst80','SpikePerActRateFirst80','FRPerCellFirst80','-append')
% save(iident,'ActRatioRateFirst80','SpikePerActRateFirst80','FRPerCellFirst80','-append')

    
    