function [hover2,move2,hovertime2,movetime2,CandStepSize] = get_hover_jumps_remotesleep(iident,irun)
             
%%%%% edited to go through the runs of the posterior and not save out but
%%%%% output
    
load(iident)

% CandStepSize: hover step number; jump step number; predict step size
% based on range (earlier version); predict step size based on sum, the
% range the replay covers in space, proportional to the length of the track
CandStepSize=zeros(size(CandSeq,1),5,2); % 2 stands for peakpos and weight pos

% duration; total travel (sum, earlier version); total travel (range); SeqNo
hover1=[];
move1=[];
hover2=[];
move2=[];
hovertime1=[];
movetime1=[];
hovertime2=[];
movetime2=[];

% Ok Index is the same, 
% but need to make sure (moving decoded window) for each decoding only the
% 4*Index(SeqNo)+1:Index(SeqNo+1)*4-3
OutMatrix = OutMatrix(:,:,irun);
InMatrix = InMatrix(:,:,irun);
    
[~,I]=max(OutMatrix'+InMatrix');
% disp(num2str(size(InMatrix,2)*replayparams.Bin_Size))
Mat=OutMatrix'+InMatrix';
I2=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);

for SeqNo=1:size(CandSeq,1)        
    for method=1:2 % peak decoded positon versus weighted position
        if method==1
            Loc=I(Index(SeqNo)*4+1:Index(SeqNo+1)*4-3);
        else
            Loc=I2(Index(SeqNo)*4+1:Index(SeqNo+1)*4-3);
        end
        DiffLoc=abs(diff(Loc));                   
        Mark=DiffLoc>2; % Ting used 1, ABL trying more to be more robust (1 = 2.5)
        Check=diff([0 Mark 0]);
        SS=find(Check==1);        
        EE=find(Check==-1);
        
        if isempty(SS)
            move=[];
            hover=[length(Mark) sum(DiffLoc) range(Loc) SeqNo 0];       
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
                    hover(i,2)=sum(DiffLoc(EE(i):length(Mark)));
                    hover(i,3)=range(Loc(EE(i):length(Mark)+1));
                end
            end              
            movetime =[SS;EE;SeqNo*ones(1,length(SS))]';
            hovertime = [EE;SS(2:end) length(Mark)+1;SeqNo*ones(1,length(SS))]';
            hover(:,5) = 1; % 1 is in the middle, 0 is on the ends
            if SS(1)>1
                hover=[[SS(1)-1 sum(DiffLoc(1:SS(1)-1)) range(Loc(1:SS(1))) SeqNo 0];hover];    
                hovertime = [[1 SS(1) SeqNo];hovertime];
            end
            if hovertime(end,2)==length(Loc)
                hover(end,5) = 0;
            end
            
        end
        CandStepSize(SeqNo,:,method)=[size(hover,1) size(move,1) range(Loc)/length(Loc) mean(DiffLoc) range(Loc)/size(InMatrix,2)];    
        eval(['hover' num2str(method) '=[hover' num2str(method) ';hover];']);
        eval(['move' num2str(method) '=[move' num2str(method) ';move];']);    
        eval(['hovertime' num2str(method) '=[hovertime' num2str(method) ';hovertime];']);
        eval(['movetime' num2str(method) '=[movetime' num2str(method) ';movetime];']);    
    end        
end   
clear hover move hovertime movetime
% save(iident,'hover*','move*','hovertime*','movetime*','CandStepSize','-append')
CandStepSize = CandStepSize(:,:,2);