function wc_only(iident)

%%%%%
%%%%% Calculate position (column-cycle) shuffle and weighted correlation
%%%%% (NOT replay score) for replays with overlapping bins.
%%%%% Same as shuffle_wc_replayscore but w/o replay score.
%%%%% Used in RunReplay_ABL_ymaze.m
%%%%%

load(iident,'params')


load(iident,'OutMatrix','InMatrix','Index','CandSeq')
% position (column-cycle) shuffle and weighted correlation

if size(InMatrix,2)>size(InMatrix,1)
    OutMatrix = OutMatrix';
    InMatrix = InMatrix';
    save(iident,'OutMatrix','InMatrix','-append')
end

D=OutMatrix'+InMatrix';   


YY=size(Index,2)-1;
Z=zeros(size(D,1), max(diff(Index)), YY);
for c=1:size(Index,2)-1
    Z(:,1:((4*Index(c+1))-(4*Index(c))),c)=D(:,4*Index(c)+1:4*Index(c+1)); 
end

denom=sum(sum(Z,1),2); 
Tx=repmat(ones(size(Z,1),1)*[1:size(Z,2)],[1 1 size(Z,3)]);

mxw=sum(sum(bsxfun(@times,Z,[1:size(Z,2)]),1),2)./denom; 
myw=sum(sum(bsxfun(@times,Z,[1:size(Z,1)]'),1),2)./denom; 
mxyw=sum(sum(bsxfun(@times,Z,[1:size(Z,1)]'*[1:size(Z,2)]),1),2)./denom; 
mxxw=sum(sum(bsxfun(@times,Z,[1:size(Z,2)].^2),1),2)./denom; 
myyw=sum(sum(bsxfun(@times,Z,[1:size(Z,1)]'.^2),1),2)./denom;

wc=(mxyw-mxw.*myw)./(sqrt(mxxw-mxw.*mxw).*sqrt(myyw-myw.*myw));
out=reshape(wc,size(wc,3),1);
disp(['Number Replays with abs(wc)>.6: ' num2str(sum(abs(out)>.6))])
CandCorr=out; 

clear c
P=size(Z,1);
S=size(Z,2);
[~,I]=max(Z,[],1);    
A=zeros(S,size(I,3));
A(:,:)=I(1,:,:);
n_shuffles = 1;
PeakProbPos=zeros(S,YY,n_shuffles+1);    

PeakProbPos(:,:,end)=A;

CandDis=zeros(size(CandSeq,1),1);
Mat=diff(PeakProbPos);
for SeqNo=1:size(CandSeq,1)
    A=max(abs(Mat(1:((4*Index(SeqNo+1))-(4*Index(SeqNo)+1))-1,SeqNo,:)),[],1);        
    CandDis(SeqNo)=A(1,1,end);    
end


save(iident,'CandCorr','CandDis','-append')

end