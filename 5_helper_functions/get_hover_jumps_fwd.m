function get_hover_jumps_fwd(iident)
             
%%%%%
%%%%%
%%%%% This script gets hovers out of replays. It was written by Ting Feng
%%%%% and edited by ABL. SImilar to get_hover_jumps but only fwd moevement
%%%%%
%%%%%

load(iident)

% CandStepSize: hover step number; jump step number; predict step size
% based on range (earlier version); predict step size based on sum, the
% range the replay covers in space, proportional to the length of the track
CandStepSize=zeros(size(CandSeq,1),5,2); % 2 stands for peakpos and weight pos

% duration; total travel (sum, earlier version); total travel (range); SeqNo
hover_fwd1=[];
move_fwd1=[];
hover_fwd2=[];
move_fwd2=[];
hovertime_fwd1=[];
movetime_fwd1=[];
hovertime_fwd2=[];
movetime_fwd2=[];

% Ok Index is the same, 
% but need to make sure (moving decoded window) for each decoding only the
% 4*Index(SeqNo)+1:Index(SeqNo+1)*4-3

[A,I]=max(OutMatrix'+InMatrix');
% disp(num2str(size(InMatrix,2)*replayparams.Bin_Size))
Mat=OutMatrix'+InMatrix';
I2=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);

for SeqNo=1:size(CandSeq,1)        
    for method=1:2 % peak decoded positon versus weighted position
        if method==1
            Loc=I(Index(SeqNo)*4+1:Index(SeqNo+1)*4-3);
%             Loc=I(Index(SeqNo)*4+1:Index(SeqNo+1)*4);
        else
            Loc=I2(Index(SeqNo)*4+1:Index(SeqNo+1)*4-3);
%             Loc=I2(Index(SeqNo)*4+1:Index(SeqNo+1)*4);
        end
%         DiffLoc=abs(diff(Loc));
        DiffLoc=diff(Loc);
        if CandCorr(SeqNo)<0
            DiffLoc = -DiffLoc;
        end
        DiffLoc(DiffLoc<0) = 0; %only taking the forward movement
        Mark=DiffLoc>1; % Ting used 1, ABL trying more to be more robust (1 = 2.5)
        Check=diff([0 Mark 0]);
        SS=find(Check==1);
        EE=find(Check==-1);
        if length(SS)==0
            move=[];
            hover=[length(Mark) sum(DiffLoc) range(Loc) SeqNo];       
            movetime = []; hovertime = [];
        else
            move=zeros(length(SS),4); % duration (time bins), total travel pos bins(sum), total travel pos bins(range); SeqNo
            hover=zeros(length(SS),4);
            move(:,1)=EE-SS;
            hover(:,1)=[SS(2:end) length(Mark)+1 ]-EE;
            move(:,4)=SeqNo;
            hover(:,4)=SeqNo;
            for i=1:length(SS)
                move(i,2)=sum(DiffLoc(SS(i):EE(i)-1));
                move(i,3)=range(Loc(SS(i):EE(i)));
                if i~=length(SS)
                    hover(i,2)=sum(DiffLoc(EE(i):SS(i+1)-1));
                    hover(i,3)=range(Loc(EE(i):SS(i+1)));
                else
                    if EE(end)>length(Mark)
                        hover(i,:)=[];
                    else
                        hover(i,2)=sum(DiffLoc(EE(i):length(Mark)));
                        hover(i,3)=range(Loc(EE(i):length(Mark)+1));
                    end
                end
            end  
            movetime =[SS;EE;SeqNo*ones(1,length(SS))]';
            hovertime = [EE;SS(2:end) length(Mark)+1;SeqNo*ones(1,length(SS))]';
            
            if SS(1)>1
                hover=[[SS(1)-1 sum(DiffLoc(1:SS(1)-1)) range(Loc(1:SS(1))) SeqNo];hover];    
                hovertime = [[1 SS(1) SeqNo];hovertime];
            end
            

        end
        CandStepSize(SeqNo,:,method)=[size(hover,1) size(move,1) range(Loc)/length(Loc) mean(DiffLoc) range(Loc)/size(InMatrix,2)];    
        eval(['hover_fwd' num2str(method) '=[hover_fwd' num2str(method) ';hover];']);
        eval(['move_fwd' num2str(method) '=[move_fwd' num2str(method) ';move];']);    
        eval(['hovertime_fwd' num2str(method) '=[hovertime_fwd' num2str(method) ';hovertime];']);
        eval(['movetime_fwd' num2str(method) '=[movetime_fwd' num2str(method) ';movetime];']);    
    end        
end   
CandStepSize_fwd = CandStepSize;
clear hover move hovertime movetime
save(iident,'hover_fwd*','move_fwd*','hovertime_fwd*','movetime_fwd*','CandStepSize_fwd','-append')
