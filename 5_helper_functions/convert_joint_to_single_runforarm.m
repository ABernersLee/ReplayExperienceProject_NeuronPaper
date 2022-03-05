function [Index,replayarms,CandSeq,InMatrix,OutMatrix,armpos] = convert_joint_to_single_runforarm(OutMatrix,InMatrix,Index,CandSeq,params,replayparams)
%goal of this is to take the decoded posteriors and break up each replay
%into which arm it is. 

%if it is two arms, split that into to single replays
%to analyze and update the Index and CandSeq to reflect that change

%also spits out linearized open field decoding (could have done other way
%around, but already had it from open field analysis) and zeros the
%posterior of the arms that it isnt



nb = [size(InMatrix,1) size(InMatrix,2)];
[x,y] = ind2sub(nb,1:(size(OutMatrix,1)*size(OutMatrix,2)));
clear nb

input = (([x' y']*replayparams.Bin_Size)+replayparams.minsave)/replayparams.cm_conv;
[Dlinpos,armind] = helper_projectarms(input,params,replayparams.Bin_Size); %p is x,y, po is    
clear input

armpos = NaN(max(Dlinpos),1);
for iarm = 1:length(params.armslength)
    armpos(unique(Dlinpos(armind==iarm))) = iarm;
end
clear armind

InMatrix2 = zeros(max(Dlinpos),size(InMatrix,3));
OutMatrix2 = zeros(max(Dlinpos),size(InMatrix,3));
for ii = 1:size(Dlinpos,1)
    InMatrix2(Dlinpos(ii),:) = squeeze(squeeze(InMatrix(x(ii),y(ii),:)))' + squeeze(InMatrix2(Dlinpos(ii),:));
    OutMatrix2(Dlinpos(ii),:) = squeeze(squeeze(OutMatrix(x(ii),y(ii),:)))' + squeeze(OutMatrix2(Dlinpos(ii),:));                
end
clear x y
InMatrix = InMatrix2;
OutMatrix = OutMatrix2;
clear InMatrix2 OutMatrix2
D = InMatrix+OutMatrix;



YY=size(Index,2)-1;
Z=zeros(size(D,1), max(diff(Index)*4), YY);
for c=1:size(Index,2)-1
    Z(:,1:((4*Index(c+1)-3)-(4*Index(c))),c)=D(:,4*Index(c)+1:4*Index(c+1)-3); 
end

MAP_arms = NaN(size(Z,2),size(Z,3),max(armpos));         
MAPraw = MAP_arms;

gauss_x=-20:1:20;
gauss_y=normpdf(gauss_x,0,1);
gaussFilter=gauss_y/sum(gauss_y);


for iarm = 1:max(armpos)    
    segment = Z(armpos==iarm,:,:); 
    MAP1 = squeeze(max(segment,[],1));    
    MAP = conv2(gaussFilter,ones*size(gaussFilter),MAP1,'same');
%     MAP=withwings(((size(gauss_x,2)-1)/2)+1:end-(size(gauss_x,2)-1)/2,:);
    MAP_arms(:,:,iarm) = MAP; 
    MAPraw(:,:,iarm) = MAP1;
end    

%find which has the max at each time point and then see how long each one stays a certain arm
MAPlong = reshape(MAP_arms,[size(MAP_arms,1)*size(MAP_arms,2) size(MAP_arms,3)]);
[~,y] = max(MAPlong,[],2);
whicharm = reshape(y,[size(MAP_arms,1) size(MAP_arms,2)]);
whicharm(sum(MAPraw==0,3)==size(MAPraw,3)) = 0;

contigabl = @(x) max([find(x==-1)-find(x==1);0]);
ApplyToGivenCol = @(func, matrix) @(col) func(matrix(:,col));
ApplyToCols = @(func, matrix) arrayfun(ApplyToGivenCol(func, matrix), 1:size(matrix,2));

numcontig = NaN(size(MAP_arms,2),size(MAP_arms,3));
for iarm = 1:size(MAP_arms,3)
   thisarm = whicharm==iarm;
   cutadd = [zeros(1,size(thisarm,2),size(thisarm,3));thisarm;zeros(1,size(thisarm,2),size(thisarm,3))];
   chn = diff(cutadd);
   numcontig(:,iarm) = ApplyToCols(contigabl,chn);
end

isrep = numcontig>10;
[~,lonrep] = max(numcontig,[],2);
%if more than 2 have map then only take the longest two
for iarm = 1:size(MAP_arms,3)
    isrep(y==iarm,iarm) = 0;
end

%this is how many arms each replay has
numrep = sum(isrep,2);

lonrep(numrep==2) = NaN;
lonrep2 = [lonrep NaN(size(lonrep,1),1)];
%for events that have two arms, split it in between them into two
Ji = find(numrep==2);
newIndex = NaN(length(Ji),1);
newSeq = NaN(length(Ji),2);
% CandSeq2 = CandSeq;
for ij = 1:length(Ji)
    arms = find(isrep(Ji(ij),:));
    whichbit = false(size(whicharm,1),2);
    for iarm = 1:2
        whichbit(:,iarm) = bwlabel(whicharm(:,Ji(ij))==arms(iarm))==find(histc(bwlabel(whicharm(:,Ji(ij))==arms(iarm)),setdiff(unique(bwlabel(whicharm(:,Ji(ij))==arms(iarm))),0))>10,1,'first');
    end
    [x,y] = find(whichbit==1);
    firstarm = y(min(x)==x);    
    nd = find(whichbit(:,firstarm)==1,1,'last');
    st = find(whichbit(:,setdiff(1:2,firstarm))==1,1,'first');
    ind2 = (st+(nd-st)/2);
    newIndex(ij) = round(ind2/4)+Index(Ji(ij));
    newtime = CandSeq(Ji(ij),1)+ind2*.005;
    oldtime = CandSeq(Ji(ij),2);
    CandSeq(Ji(ij),2) = newtime;
    newSeq(ij,:) = [newtime oldtime];
    lonrep2(Ji(ij),:) = [arms(firstarm) arms(setdiff(1:2,firstarm))];
end

Index = sort([Index newIndex']);
CandSeq = sort([CandSeq; newSeq],1);
replayarms = lonrep2(:);
replayarms(isnan(replayarms)) = [];

II = InMatrix;
OO = OutMatrix;
for c=1:size(Index,2)-1
    II(armpos~=replayarms(c),4*Index(c)+1:4*Index(c+1)) = 0;
    OO(armpos~=replayarms(c),4*Index(c)+1:4*Index(c+1)) = 0;
end

InMatrix = II'; OutMatrix = OO';
% if alreadydone
%     save(iident,'replayarms','-append')
% else
%     save(iident,'Index','replayarms','CandSeq','InMatrix','OutMatrix','-append')
% end

