function get_hover_jumps_small(iident)
             
    
load(iident)

% CandStepSize: hover step number; jump step number; predict step size
% based on range (earlier version); predict step size based on sum, the
% range the replay covers in space, proportional to the length of the track
CandStepSize=zeros(size(CandSeq,1),5,2); % 2 stands for peakpos and weight pos

% duration; total travel (sum, earlier version); total travel (range); SeqNo
hover_sm1=[];
move_sm1=[];
hover_sm2=[];
move_sm2=[];
hovertime_sm1=[];
movetime_sm1=[];
hovertime_sm2=[];
movetime_sm2=[];

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
        DiffLoc=abs(diff(Loc));                   
        Mark=DiffLoc>1; % 1 instead of 2 here (1 = 2.5cm)
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
        eval(['hover_sm' num2str(method) '=[hover_sm' num2str(method) ';hover];']);
        eval(['move_sm' num2str(method) '=[move_sm' num2str(method) ';move];']);    
        eval(['hovertime_sm' num2str(method) '=[hovertime_sm' num2str(method) ';hovertime];']);
        eval(['movetime_sm' num2str(method) '=[movetime_sm' num2str(method) ';movetime];']);    
    end        
end   
CandStepSize_sm = CandStepSize;
clear hover move hovertime movetime
save(iident,'hover_sm*','move_sm*','hovertime_sm*','movetime_sm*','CandStepSize_sm','-append')
