function [CandCorr,CandDis,CandRS] = shuffle_wc_replayscore_remotesleep(iident,irun)

%%%%%
%%%%% Calculate position (column-cycle) shuffle and weighted correlation
%%%%% and replay score for replays with overlapping bins.
%%%%% Same as shuffle_wc_replayscore.m but goes through runs in posterior
%%%%% Used in RunReplay_ABL_RemoteSleep.m
%%%%%

    load(iident,'OutMatrix','InMatrix','Index','CandSeq')
    OutMatrix = OutMatrix(:,:,irun);
    InMatrix = InMatrix(:,:,irun);
    
    if size(InMatrix,2)>size(InMatrix,1)
        OutMatrix = OutMatrix';
        InMatrix = InMatrix';
        save(iident,'OutMatrix','InMatrix','-append')
    end
    
    D=OutMatrix'+InMatrix';   

   
    YY=size(Index,2)-1;
    Z=NaN(size(D,1), max(diff(Index))*4, YY);
    for c=1:size(Index,2)-1
        Z(:,1:((4*Index(c+1)-3)-(4*Index(c))),c)=D(:,4*Index(c)+1:4*Index(c+1)-3); 
%         Z(:,1:((4*Index(c+1))-(4*Index(c))),c)=D(:,4*Index(c)+1:4*Index(c+1)); 
    end

    %nonoverlapping
%     Z=zeros(size(D,1), max(diff(Index)), YY);
%     for c=1:size(Index,2)-1
%        Z(:,1:(Index(c+1)-Index(c)),c)=D(:,Index(c)+1:Index(c+1)); 
%     end

    denom=nansum(nansum(Z,1),2); 
    Tx=repmat(ones(size(Z,1),1)*[1:size(Z,2)],[1 1 size(Z,3)]);

    mxw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,2)]),1),2)./denom; 
    myw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,1)]'),1),2)./denom; 
    mxyw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,1)]'*[1:size(Z,2)]),1),2)./denom; 
    mxxw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,2)].^2),1),2)./denom; 
    myyw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,1)]'.^2),1),2)./denom;

    wc=(mxyw-mxw.*myw)./(sqrt(mxxw-mxw.*mxw).*sqrt(myyw-myw.*myw));
%     ZZ = Z(:,4*Index(c)+1:4*Index(c+1),1);
%     r = weighted_corr([1:size(ZZ,2)],[1:size(ZZ,1)]',ZZ); %testing
    out=reshape(wc,size(wc,3),1);
    
    
    %shuffles
    n_shuffles=500; 

    PosShuffleCorr=zeros(size(CandSeq,1),n_shuffles);
    
    P=size(Z,1);
    S=size(Z,2);
    PosShift=zeros(S,n_shuffles);
    
    [nanI,I2]=max(Z,[],1);    
%     I2=sum([1:size(Z,1)]'*ones(1,size(Z,2)).*Z); %ABL changed from max to weighted 6/18
    A=zeros(S,size(I2,3));
    A(:,:)=I2(1,:,:);
    A(isnan(nanI)) = NaN;
    

    %examples
%     for SeqNo = [4 48 68]
%         x = find(~isnan(II(:,SeqNo)));
%         y = II(~isnan(II(:,SeqNo)),SeqNo);
%         [p1,S1] = polyfit(x,y,1);
%         R2_1 = 1 - (S1.normr/norm(y - mean(y)))^2;
%         R2_11 = 1-( (1- R2_1)*(length(x)-1)./(length(x)-1-1));
%         [p,S] = polyfit(x,y,3);
%         R2_3 = 1 - (S.normr/norm(y - mean(y)))^2;
%         R2_33 = 1-( (1- R2_3)*(length(x)-1)./(length(x)-3-1));
%         [SeqNo R2_3-R2_1 R2_3/R2_1 R2_33-R2_11 R2_33/R2_11]    
%         y1 = polyval(p,x);
%         figure; hold on; plot(x,y,'k.','MarkerSize',20); 
%         plot(x,y1,'r')
%         y1 = polyval(p1,x);
%         plot(x,y1,'b')
%         text(1,20,['Radj = ' num2str(round(R2_11,2,'significant'))],'Color','b')
%         text(1,15,['Radj = ' num2str(round(R2_33,2,'significant'))],'Color','r')
%         text(1,10,['Diff = ' num2str(round(R2_33-R2_11,2,'significant'))],'Color','k')
%         yl = get(gca,'ylim');
%         set(gca,'ylim',[0 yl(2)])
%     %     p
%     end

    PeakProbPos=zeros(S,YY,n_shuffles+1,'single');
    PeakProbPos(:,:,end)=A;

    for s=1:n_shuffles
        
        a=rand(S,1);
        PosShift(:,s)=floor(P*a);
        PeakProbPos(:,:,s)=mod(A+floor(P*a)*ones(1,YY),P)+1;
        
        no=mod(floor(P*a)*ones(1,P)+ones(S,1)*[1:P],P)'+1;               

        myw=nansum(nansum(bsxfun(@times,Z,no),1),2)./denom; 
        mxyw=nansum(nansum(bsxfun(@times,Z,no.*Tx(:,:,1)),1),2)./denom;
        myyw=nansum(nansum(bsxfun(@times,Z,no.^2),1),2)./denom;

        wc=(mxyw-mxw.*myw)./(sqrt(mxxw-mxw.*mxw).*sqrt(myyw-myw.*myw));
        
        Sout=reshape(wc,size(wc,3),1);
        
        PosShuffleCorr(:,s)=Sout;
     
    end
    
    CandCorr=out; 
    
    clear c
    PosShuffleDis=zeros(size(CandSeq,1),n_shuffles);
    CandDis=zeros(size(CandSeq,1),1);
    Mat=diff(PeakProbPos);
    for SeqNo=1:size(CandSeq,1)
        A=max(abs(Mat(1:((4*Index(SeqNo+1)-3)-(4*Index(SeqNo)))-1,SeqNo,:)),[],1);                    
        CandDis(SeqNo)=A(1,1,end);
        PosShuffleDis(SeqNo,:)=A(1,1,1:n_shuffles);       
    end
    clear c s S no Ty myw Txy mxyw myyw wc Sout out Z Tx mxw mxxw PeakProbPos a Mat SeqNo A I 
    clear n_shuffles
        
    disp(['Number Replays with abs(wc)>.6 & jd<.4: ' num2str(sum(abs(CandCorr)>.6 & CandDis/size(OutMatrix,2)<.4))])
 %    global i AA BB    
%     AA = 851*4; %big
%     BB = 101; %big
%     AA = 851;

    %%%overlapping
    AA = 150;
    BB = 100;
    Z=zeros(AA,BB, YY);
    for c=1:YY
        len = (Index(c+1)-Index(c))*4-3;
        leftover = BB-len;
        Z(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((leftover/2))+1:floor((BB/2-len/2))+len,c)=D(:,4*Index(c)+1:4*Index(c+1)-3); 
    end    

    InterS=[-10:.2:10]; %slopes
    InterP=[-35:35]; %intersepts    
    RS=NaN(length(InterS),length(InterP),YY);
    d = 10;
    tic
    for i=1:length(InterS)
        for j=1:length(InterP)
            
            k = 0:BB-1;
            traj = (k-BB/2)*InterS(i)+InterP(j)+AA/2;
%             offtrack = or(traj<1,traj>AA);
            pos_lim = [(floor(traj-d))' (ceil(traj+d))'];
            pos_lim(pos_lim<1) = 1;
            pos_lim(pos_lim>AA) = AA;
            idx_lo = sub2ind([AA BB],pos_lim(:,1)',k+1);
            idx_hi = sub2ind([AA BB],pos_lim(:,2)',k+1);   
%             temp = zeros(length(k),1);
            Mark=zeros(AA,BB); 
            for x = 1:length(k)
                Mark(idx_lo(x):idx_hi(x)) = 1;
            end
    %         temp(offtrack) = median(decoded_matrix(:,offtrack),1);
            
            RS(i,j,:)=nansum(nansum(bsxfun(@times,Z,Mark),1),2)./denom;                        
    %             ApplyInd(@getRS,InterP); for arrayfun insead of loop
        end
    end
    
    [A,I]=nanmax(RS);
    [B,I2]=nanmax(A);
    
    tt = toc;
    disp([num2str(round(tt/60,3,'significant')) ' minutes for making replay score'])
    
    
    CandRS=zeros(size(CandSeq,1),4); % Replay Score : RS | Slope | InterSec | RSorder 
    CandRS(:,1:3)=[B(:) InterS(I(I2(:)+([1:size(I,3)]'-1)*size(I,2)))' InterP(I2(:))'];    
%     figure; imagesc(Mark); axis xy
    clear Z
    
    D=OutMatrix'; 
    Z1=zeros(AA, BB, YY);
    for c=1:size(Index,2)-1            
        len = (Index(c+1)-Index(c))*4-3;
        leftover = BB-len;
        Z1(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((leftover/2))+1:floor((BB/2-len/2))+len,c)=D(:,4*Index(c)+1:4*Index(c+1)-3);       
    end

    clear D
    D=InMatrix'; 
    Z2=zeros(AA, BB, YY);
    for c=1:size(Index,2)-1            
         len = (Index(c+1)-Index(c))*4-3;
        leftover = BB-len;
        Z2(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((leftover/2))+1:floor((BB/2-len/2))+len,c)=D(:,4*Index(c)+1:4*Index(c+1)-3); 
    end
%     
    Mark=zeros(AA,BB,YY);
    for SeqNo=1:size(CandSeq,1)
       k = 0:BB-1;
        traj = (k-BB/2)*CandRS(SeqNo,2)+CandRS(SeqNo,3)+AA/2;
        pos_lim = [(floor(traj-d))' (ceil(traj+d))'];
        pos_lim(pos_lim<1) = 1;
        pos_lim(pos_lim>AA) = AA;
        idx_lo = sub2ind([AA BB],pos_lim(:,1)',k+1);
        idx_hi = sub2ind([AA BB],pos_lim(:,2)',k+1);   
        Mark1=zeros(AA,BB); 
        for x = 1:length(k)
            Mark1(idx_lo(x):idx_hi(x)) = 1;
        end
        Mark(:,:,SeqNo) = Mark1;
    end    
    
    CandRS(:,4)=(nansum(nansum(Z1.*Mark,1),2)-nansum(nansum(Z2.*Mark,1),2))./(nansum(nansum(Z1.*Mark,1),2)+nansum(nansum(Z2.*Mark,1),2));        

end