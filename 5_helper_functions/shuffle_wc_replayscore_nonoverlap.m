function shuffle_wc_replayscore_nonoverlap(iident)


%%%%%
%%%%% Calculate timebin shuffle and weighted correlation
%%%%% and replay score for replays with non-overlapping bins.
%%%%% Used in RunReplay_ABL_RemoteSleep.m, and RunReplay_ABL_OneLap.m
%%%%%


%an older version with commented out parts including an old circular shift
%of position (column-cycle shuffle) called
%shuffle_wc_replayscore_nonoverlap - Old - changed July2020


    load(iident,'OutMatrix','InMatrix','Index','CandSeq')
 % position (column-cycle) shuffle and weighted correlation
    
    if size(InMatrix,2)>size(InMatrix,1)
        OutMatrix = OutMatrix';
        InMatrix = InMatrix';
        save(iident,'OutMatrix','InMatrix','-append')
    end
    
    D=OutMatrix'+InMatrix';   

   
   
    YY=size(Index,2)-1;
     %nonoverlapping
    Z=NaN(size(D,1), max(diff(Index)), YY);
    for c=1:size(Index,2)-1
       Z(:,1:(Index(c+1)-Index(c)),c)=D(:,Index(c)+1:Index(c+1)); 
    end    
    
    denom=nansum(nansum(Z,1),2); 
    mxw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,2)]),1),2)./denom; 
    myw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,1)]'),1),2)./denom; 
    mxyw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,1)]'*[1:size(Z,2)]),1),2)./denom; 
    mxxw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,2)].^2),1),2)./denom; 
    myyw=nansum(nansum(bsxfun(@times,Z,[1:size(Z,1)]'.^2),1),2)./denom;
    wc=(mxyw-mxw.*myw)./(sqrt(mxxw-mxw.*mxw).*sqrt(myyw-myw.*myw));
    CandCorr=reshape(wc,size(wc,3),1);
    denomsave = denom;
    
    %shuffles
    n_shuffles=5000;
    P=size(Z,1);
    S=size(Z,2);
    
    
    [nanI,I2]=max(Z,[],1);    
%     I2=sum([1:size(Z,1)]'*ones(1,size(Z,2)).*Z); %ABL changed from max to weighted 6/18
    A=zeros(S,size(I2,3));
    A(:,:)=I2(1,:,:);
    A(isnan(nanI)) = NaN;
    CandDis = max(abs(diff(A)),[],1);
    CandDis = CandDis';
    
    if 0
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
    end
    
    TmShuffleCorr=NaN(size(CandSeq,1),n_shuffles);
    TmShuffleDis=TmShuffleCorr;
    for SeqNo=1:size(CandSeq,1)
        len = Index(SeqNo+1)-Index(SeqNo);
        [~,RandInd] = sort(rand(n_shuffles,len),2);
        RandInd = repmat(RandInd,[1 1 size(Z,1)]);
        RandInd = permute(RandInd,[2 3 1]); 
        addon = permute(repmat(len*[1:size(Z,1)]-len,[len 1 n_shuffles]),[1 2 3])+permute(repmat([1:n_shuffles]'*(len*size(Z,1))-(len*size(Z,1)),[1 len size(Z,1)]),[2 3 1]);
        RandInd = RandInd+addon;        
        Z1 = repmat(Z(:,1:len,SeqNo),[1 1 n_shuffles]);
        Z1 = permute(Z1,[2 1 3]);
%         Z2 = cat(3,Z1(RandInd),Z1(:,:,1)); %testing
        Z2 = Z1(RandInd);
        Z2= permute(Z2,[2 1 3]);
        denom=nansum(nansum(Z2,1),2); 
        mxw=nansum(nansum(bsxfun(@times,Z2,[1:size(Z2,2)]),1),2)./denom; 
        myw=nansum(nansum(bsxfun(@times,Z2,[1:size(Z2,1)]'),1),2)./denom; 
        mxyw=nansum(nansum(bsxfun(@times,Z2,[1:size(Z2,1)]'*[1:size(Z2,2)]),1),2)./denom; 
        mxxw=nansum(nansum(bsxfun(@times,Z2,[1:size(Z2,2)].^2),1),2)./denom; 
        myyw=nansum(nansum(bsxfun(@times,Z2,[1:size(Z2,1)]'.^2),1),2)./denom;
        wc=(mxyw-mxw.*myw)./(sqrt(mxxw-mxw.*mxw).*sqrt(myyw-myw.*myw));    
        TmShuffleCorr(SeqNo,:)=reshape(wc,size(wc,3),1);
        
        [nanI,I2]=max(Z2,[],1);  
        A=zeros(size(Z2,2),size(I2,3));
        A(:,:)=I2(1,:,:);
        A(isnan(nanI)) = NaN;        
        TmShuffleDis(SeqNo,:) = max(abs(diff(A)),[],1);
        if mod(SeqNo,500)==0
            disp(['Done with ' num2str(SeqNo) ' of ' num2str(size(CandSeq,1)) ' Events'])
        end
    end
    

    sigst1 = ((sum(((TmShuffleDis))<=CandDis,2)+1)./(size(TmShuffleDis,2)+1))<.05;
    sigst2 = ((sum(abs(TmShuffleCorr)>=abs(CandCorr),2)+1)./(size(TmShuffleCorr,2)+1))<.05;
    disp([num2str(sum(sigst1 & sigst2)) ' Replays pass significance: ' num2str(round(100*(sum(sigst1 & sigst2)./size(CandCorr,1)))) '%'])
    disp(['Number Replays with abs(wc)>.6 & jd<.4: ' num2str(sum(abs(CandCorr)>.6 & (CandDis/size(InMatrix,2))<.4))])

    
%%% Add replay score   
    Z=zeros(851, 25, YY);
    for c=1:size(Index,2)-1
       Z(floor((851-P)/2)+1:floor((851-P)/2)+P,floor((25-Index(c+1)+Index(c))/2)+1:floor((25-Index(c+1)+Index(c))/2)+Index(c+1)-Index(c),c)=D(:,Index(c)+1:Index(c+1)); 
       Z([1:floor((851-P)/2) floor((851-P)/2)+P+1:end] ,floor((25-Index(c+1)+Index(c))/2)+1:floor((25-Index(c+1)+Index(c))/2)+Index(c+1)-Index(c),c)=ones(851-P,1)*median(D(:,Index(c)+1:Index(c+1)));          
    end
    
    InterS=[-30:30];
    InterP=[-42:42];
    
    RS=zeros(length(InterS),length(InterP),YY);
    for i=1:length(InterS)
        for j=1:length(InterP)
            Mark=zeros(851,25);
            p=round(ones(13,1)*(InterP(j)+InterS(i)*[-12:12]+426)+[-6:6]'*ones(1,25));
            Mark(p+ones(13,1)*[0:24]*851)=1;  
            RS(i,j,:)=sum(sum(bsxfun(@times,Z,Mark),1),2)./denomsave;
        end
    end
    
    [A,I]=max(RS);
    [B,I2]=max(A);
    
    CandRS=zeros(size(CandSeq,1),4); % Replay Score : RS | Slope | InterSec | RSorder 
    CandRS(:,1:3)=[B(:) InterS(I(I2(:)+([1:size(I,3)]'-1)*size(I,2)))' InterP(I2(:))'];
    tt = toc;
    disp([num2str(round(tt/60,3,'significant')) ' minutes for making replay score'])
    clear Z
%     disp(['Done making rest of replay score in ' num2str(round(tt/60,3,'significant')) ' minutes'])
        
%     save(iident,'Radjusted','KendallR','SpearmanR','CandCorr','PosShuffleCorr','CandDis','PosShuffleDis','PosShift','CandRS','-append')
    save(iident,'CandCorr','TmShuffleCorr','CandDis','TmShuffleDis','CandRS','-append')
    
%     save(iident,'CandCorr','PosShuffleCorr','CandDis','PosShuffleDis','PosShift','-append')
end