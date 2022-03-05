function shuffle_wc_replayscore_Epoch1_overlap(iident)

%%%%%
%%%%% Calculate position (column-cycle) shuffle and weighted correlation
%%%%% and replay score for replays with overlapping bins.
%%%%% Same as shuffle_wc_replayscore.m except uses MatrixEpoch1
%%%%% Used in RunReplay_ABL_OneLap.m
%%%%%

    load(iident,'MatrixEpoch1_overlap','Index_overlap','CandSeq_overlap')
 % position (column-cycle) shuffle and weighted correlation
    
    if size(MatrixEpoch1_overlap,2)>size(MatrixEpoch1_overlap,1)
        MatrixEpoch1_overlap = MatrixEpoch1_overlap';
        save(iident,'MatrixEpoch1_overlap','-append')
    end
    CandSeq = CandSeq_overlap;
    Index = Index_overlap;
    D=MatrixEpoch1_overlap';   

   
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
    n_shuffles=5;

    PosShuffleCorr=zeros(size(CandSeq,1),n_shuffles);
    
    P=size(Z,1);
    S=size(Z,2);
    PosShift=zeros(S,n_shuffles);
    
    [nanI,I2]=max(Z,[],1);    
%     I2=sum([1:size(Z,1)]'*ones(1,size(Z,2)).*Z); %ABL changed from max to weighted 6/18
    A=zeros(S,size(I2,3));
    A(:,:)=I2(1,:,:);
    A(isnan(nanI)) = NaN;

    %     SpearmanR = corr(squeeze(I2),[1:size(I2,2)]'*ones(1,size(I2,3)),'type','Spearman','rows','complete');
%     SpearmanR = SpearmanR(:,1);
%     KendallR = corr(squeeze(I2),[1:size(I2,2)]'*ones(1,size(I2,3)),'type','Kendall','rows','complete');
%     KendallR = KendallR(:,1);
%     II = squeeze(I2);
%         
%     Radjusted = NaN(size(II,2),3);
%     for SeqNo = 1:size(II,2)
%         x = find(~isnan(II(:,SeqNo)));
%         y = II(~isnan(II(:,SeqNo)),SeqNo);
%         for pp = 1:3
%             [~,S1] = polyfit(x,y,pp);
%             R2_1 = 1 - (S1.normr/norm(y - mean(y)))^2;
%             Radjusted(SeqNo,pp) = 1-( (1- R2_1)*(length(x)-1)./(length(x)-pp-1));
%         end
%     end


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
        A=max(abs(Mat(1:Index(SeqNo+1)-Index(SeqNo)-1,SeqNo,:)),[],1);
%         A=max(abs(Mat(1:((4*Index(SeqNo+1)-3)-(4*Index(SeqNo)))-1,SeqNo,:)),[],1);        
%         A=max(abs(Mat(1:((4*Index(SeqNo+1))-(4*Index(SeqNo)))-1,SeqNo,:)),[],1);                
        CandDis(SeqNo)=A(1,1,end);
        PosShuffleDis(SeqNo,:)=A(1,1,1:n_shuffles);       
    end
    
    clear c s S no Ty myw Txy mxyw myyw wc Sout out Z Tx mxw mxxw PeakProbPos a Mat SeqNo A I 
%     disp(['Done with ' num2str(n_shuffles) ' shuffles in ' num2str(tt/60) ' minutes'])
    clear n_shuffles
        
    disp(['Number Replays with abs(wc)>.6 & jd<.4: ' num2str(sum(abs(CandCorr)>.6 & CandDis/size(MatrixEpoch1_overlap,2)<.4))])
 %    global i AA BB    
%     AA = 851*4; %big
%     BB = 101; %big
%     AA = 851;

% %     %with overlapping 
    AA = 150;
    BB = 100;
    Z=zeros(AA,BB, YY);
    for c=1:YY
%         Z(floor((851-P)/2)+1:floor((851-P)/2)+P,floor((BB-Index(c+1)+Index(c))/2)+1:floor((BB-Index(c+1)+Index(c))/2)+Index(c+1)-Index(c),c)=D(:,Index(c)+1:Index(c+1)); 
%        Z([1:floor((851-P)/2) floor((851-P)/2)+P+1:end] ,floor((BB-Index(c+1)+Index(c))/2)+1:floor((BB-Index(c+1)+Index(c))/2)+Index(c+1)-Index(c),c)=ones(851-P,1)*median(D(:,Index(c)+1:Index(c+1)));  
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

    
%     D=OutMatrix'; 
%     Z1=zeros(AA, BB, YY);
%     for c=1:size(Index,2)-1            
% %        Z1(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((BB-(4*Index(c+1))+(Index(c)*4))/2)+1:floor((BB-(4*Index(c+1))+(Index(c)*4))/2)+(4*Index(c+1))-(Index(c)*4),c)=D(:,(Index(c)*4)+1:(4*Index(c+1))); 
% %        Z1([1:floor((AA-P)/2) floor((AA-P)/2)+P+1:end] ,floor((BB-(4*Index(c+1))+(Index(c)*4))/2)+1:floor((BB-(4*Index(c+1))+(Index(c)*4))/2)+(4*Index(c+1))-(Index(c)*4),c)=ones(AA-P,1)*median(D(:,(Index(c)*4)+1:(4*Index(c+1))));          
%         len = (Index(c+1)-Index(c))*4-3;
%         leftover = BB-len;
%         Z1(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((leftover/2))+1:floor((BB/2-len/2))+len,c)=D(:,4*Index(c)+1:4*Index(c+1)-3);       
%     end
% 
%     clear D
%     D=InMatrix'; 
%     Z2=zeros(AA, BB, YY);
%     for c=1:size(Index,2)-1            
% %        Z2(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((BB-(4*Index(c+1)-3)+(Index(c)*4))/2)+1:floor((BB-(4*Index(c+1)-3)+(Index(c)*4))/2)+(4*Index(c+1)-3)-(Index(c)*4),c)=D(:,(Index(c)*4)+1:(4*Index(c+1)-3)); 
% %        Z2([1:floor((AA-P)/2) floor((AA-P)/2)+P+1:end] ,floor((BB-(4*Index(c+1)-3)+(Index(c)*4))/2)+1:floor((BB-(4*Index(c+1)-3)+(Index(c)*4))/2)+(4*Index(c+1)-3)-(Index(c)*4),c)=ones(AA-P,1)*median(D(:,(Index(c)*4)+1:(4*Index(c+1)-3)));          
%          len = (Index(c+1)-Index(c))*4-3;
%         leftover = BB-len;
%         Z2(floor((AA-P)/2)+1:floor((AA-P)/2)+P,floor((leftover/2))+1:floor((BB/2-len/2))+len,c)=D(:,4*Index(c)+1:4*Index(c+1)-3); 
%     end
%     
%     Mark=zeros(AA,BB,YY);
%     for SeqNo=1:size(CandSeq,1)
% %         p=round(ones(13,1)*(CandRS(SeqNo,3)+CandRS(SeqNo,2)*[-(BB-1)/2:(BB-1)/2]+(AA+1)/2)+[-(BB-1)/4:(BB-1)/4]'*ones(1,BB));
% %         Mark(p+ones(13,1)*[0:BB-1]*AA+AA*BB*(SeqNo-1))=1;          
%         k = 0:BB-1;
%         traj = (k-BB/2)*CandRS(SeqNo,2)+CandRS(SeqNo,3)+AA/2;
% %             offtrack = or(traj<1,traj>AA);
%         pos_lim = [(floor(traj-d))' (ceil(traj+d))'];
%         pos_lim(pos_lim<1) = 1;
%         pos_lim(pos_lim>AA) = AA;
%         idx_lo = sub2ind([AA BB],pos_lim(:,1)',k+1);
%         idx_hi = sub2ind([AA BB],pos_lim(:,2)',k+1);   
% %             temp = zeros(length(k),1);
%         Mark1=zeros(AA,BB); 
%         for x = 1:length(k)
%             Mark1(idx_lo(x):idx_hi(x)) = 1;
%         end
%         Mark(:,:,SeqNo) = Mark1;
%     end    
%     
%     CandRS(:,4)=(nansum(nansum(Z1.*Mark,1),2)-nansum(nansum(Z2.*Mark,1),2))./(nansum(nansum(Z1.*Mark,1),2)+nansum(nansum(Z2.*Mark,1),2));        
    
 
%     CandRS(:,4)=(sum(sum(Z1.*Mark,1),2)-sum(sum(Z2.*Mark,1),2))./(sum(sum(Z1.*Mark,1),2)+sum(sum(Z2.*Mark,1),2));        
%     tt = toc;
%     disp(['Done making rest of replay score in ' num2str(round(tt/60,3,'significant')) ' minutes'])

%%% Add replay score 
%     Z=zeros(851, 25, YY);
%     for c=1:size(Index,2)-1
%        Z(floor((851-P)/2)+1:floor((851-P)/2)+P,floor((25-Index(c+1)+Index(c))/2)+1:floor((25-Index(c+1)+Index(c))/2)+Index(c+1)-Index(c),c)=D(:,Index(c)+1:Index(c+1)); 
%        Z([1:floor((851-P)/2) floor((851-P)/2)+P+1:end] ,floor((25-Index(c+1)+Index(c))/2)+1:floor((25-Index(c+1)+Index(c))/2)+Index(c+1)-Index(c),c)=ones(851-P,1)*median(D(:,Index(c)+1:Index(c+1)));          
%     end
%     
%     InterS=[-30:30];
%     InterP=[-42:42];
%     
%     RS=zeros(length(InterS),length(InterP),YY);
%     for i=1:length(InterS)
%         for j=1:length(InterP)
%             Mark=zeros(851,25);
%             p=round(ones(13,1)*(InterP(j)+InterS(i)*[-12:12]+426)+[-6:6]'*ones(1,25));
%             Mark(p+ones(13,1)*[0:24]*851)=1;  
%             RS(i,j,:)=sum(sum(bsxfun(@times,Z,Mark),1),2)./denom;
%         end
%     end
%     
%     [A,I]=max(RS);
%     [B,I2]=max(A);
%     
%     CandRS=zeros(size(CandSeq,1),4); % Replay Score : RS | Slope | InterSec | RSorder 
%     CandRS(:,1:3)=[B(:) InterS(I(I2(:)+([1:size(I,3)]'-1)*size(I,2)))' InterP(I2(:))'];
%     tt = toc;
%     disp([num2str(round(tt/60,3,'significant')) ' minutes for making replay score'])
%     clear Z
%     
%     CandRSEpoch1_overlaps = CandRS;    
%     CandCorrEpoch1_overlaps = CandCorr;
%     PosShuffleCorrEpoch1_overlaps = PosShuffleCorr; 
%     CandDisEpoch1_overlaps = CandDis;
%     PosShuffleDisEpoch1_overlaps = PosShuffleDis;
%     PosShiftEpoch1_overlaps = PosShift;
%     
%     
%     save(iident,'CandCorrEpoch1','PosShuffleCorrEpoch1','CandDisEpoch1','PosShuffleDisEpoch1','PosShiftEpoch1','-append')
    CandCorrEpoch1_overlap = CandCorr; 
    PosShuffleCorrEpoch1_overlap = PosShuffleCorr;
    CandDisEpoch1_overlap = CandDis;
    PosShuffleDisEpoch1_overlap = PosShuffleDis;
    PosShiftEpoch1_overlap = PosShift;
    CandRSEpoch1_overlap = CandRS;
    save(iident,'CandCorrEpoch1_overlap','PosShuffleCorrEpoch1_overlap','CandDisEpoch1_overlap','PosShuffleDisEpoch1_overlap','PosShiftEpoch1_overlap','CandRSEpoch1_overlap','-append')
end