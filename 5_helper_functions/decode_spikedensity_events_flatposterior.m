function [OutMatrixShuff,InMatrixShuff] =  decode_spikedensity_events_flatposterior(dirs,iident,numshuff,toplot,final)

load(iident,'params','hp_cells','hpinterneurons','pos','vel','dirdat','spikedata','cm_conv','linposcat')

%%%%%
%%%%% Make decoding with flattened FR distribution in a similar way to 
%%%%% replays to compare. 
%%%%% Used in Hovers_Moving_laps.m for Supplemental figure 6. 
%%%%%

%dirdat 1 is up/out, 0 is down/in
Bin_Size=2.5;
Cell_Number = hp_cells(~ismember(hp_cells,hpinterneurons));
disp(['Number of cells: ' num2str(length(Cell_Number))])
VelThresh=5;



%make directional place fields - linear
Pos = linposcat;
Spike=[spikedata(:,2) spikedata(:,1) Pos(spikedata(:,3)) vel(spikedata(:,3)) dirdat(spikedata(:,3))];

Number_Of_Bins=ceil(range(Pos)/Bin_Size);
OutFR=zeros(length(Cell_Number),Number_Of_Bins);
InFR=zeros(length(Cell_Number),Number_Of_Bins);

Filter=fspecial('gaussian',[1 20],2); % ref Davidson et.al. 2006
 for i=1:length(Cell_Number)
    t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==1);
    if length(t)~=1
        FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==1),min(Pos):Bin_Size:max(Pos));   
    else
        FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh  & dirdat==1),min(Pos):Bin_Size:max(Pos));   
    end
    FR(isnan(FR) | isinf(FR))=0;
    FR=FR(1:Number_Of_Bins);
    OutFR(i,:)=filter2(Filter,FR');  
%     OutFR(i,:)=FR;  

    t=find(Spike(:,4)>VelThresh & Spike(:,1)==Cell_Number(i) & Spike(:,5)==0);
    if length(t)~=1
        FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))./histc(Pos(vel>VelThresh & dirdat==0),min(Pos):Bin_Size:max(Pos));   
    else
        FR=60*histc(Spike(t,3),min(Pos):Bin_Size:max(Pos))'./histc(Pos(vel>VelThresh & dirdat==0),min(Pos):Bin_Size:max(Pos));   
    end        
    FR(isnan(FR) | isinf(FR))=0;
    FR=FR(1:Number_Of_Bins);
    InFR(i,:)=filter2(Filter,FR');
%     InFR(i,:)=FR;
end
clear Pos

% add noise in a biased mway that makes the density of fields flat for both
% directions

OutFRsave = OutFR;
InFRsave = InFR;

% [~,m] = max(OutFR');
% [~,o1] = sort(m);
% figure; imagesc(OutFR(o,:)./max(OutFR(o,:),[],2))
% num = 0;
% while any(mean(OutFR)~=mean(mean(OutFR)))    
%     num = num+1;
%       bias = repmat(mean(OutFR)-mean(mean(OutFR)),[size(OutFR,1) 1]);
%       OutFR = OutFR+((-bias).*rand(size(OutFR)));
%       OutFR(OutFR<0) = .00001;
% end
% % [~,m] = max(OutFR');
% % [~,o] = sort(m);
% % figure; imagesc(OutFR(o,:)./max(OutFR(o,:),[],2))
% % figure; imagesc(OutFR(o1,:)./max(OutFR(o1,:),[],2))
% %       
% num2 = 0;
% while any(mean(InFR)~=mean(mean(InFR)))    
%     num2 = num2+1;
%       bias = repmat(mean(InFR)-mean(mean(InFR)),[size(InFR,1) 1]);
%       InFR = InFR+((-bias).*rand(size(InFR)));
%       InFR(InFR<0) = .00001;
% end

% OutFRv1 = OutFR;
% InFRv1 = InFR;
                                    % OutFR = OutFRsave; InFR = InFRsave;
clear Number_Of_Bins i Filter  dirdat hp_cells hpinterneurons

% get spike density events


load(iident,'CandSeq')



% decode events


EstBin=0.02;

Cand=CandSeq;
N=round(diff(CandSeq')/EstBin)';
t=find(mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
t=find(~mod(N,2));
Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
clear N t

[A,I]=histc(Spike(:,2),sortrows(Cand(:)));

B=cumsum([0 diff(Cand')]);

Orig_Spike=Spike;
    
for i=1:2:max(I)
    Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); % HERE
end
Spike(~mod(I,2),:)=[];
clear Cand I i

Start_Time=0;
End_Time=B(end);
TimeBins=round((End_Time-Start_Time)/EstBin)*4;
% Cell_Number=size(OutFR,1);

% Bayesian Decoding - 5ms moving step, 20ms window estimate
binspike=zeros(length(Cell_Number),TimeBins);
for CellID=1:length(Cell_Number)    
    for i=1:4
        c=histc(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
%         c=histcounts(Spike(Spike(:,1)==Cell_Number(CellID),2),Start_Time+(i-1)*EstBin/4:EstBin:Start_Time+EstBin/4*TimeBins);
        binspike(CellID,i:4:TimeBins)=c(1:TimeBins/4);
    end                      
end

Number_Of_Bins=size(OutFR,2);
clear Start_Time End_Time CellID c i 
% shuffle
InMatrixShuff = NaN(TimeBins,Number_Of_Bins,numshuff);
OutMatrixShuff = InMatrixShuff;        
    
for ishuff = 1:numshuff

    InFR(sum(InFR~=0,2)==0,:) = 1;
    OutFR(sum(OutFR~=0,2)==0,:) = 1;    
    InFR(InFR==0) = .0001;
    OutFR(OutFR==0) = .0001;


    
%     OutFRs = OutFRsave;
%     InFRs = InFRsave;
    
    %faster but produces negative spots in the fields
%     tic
% %     bias = repmat(mean(mean(OutFRs))-mean(OutFRs),[size(OutFRs,1) 1]);
%     bias = repmat(sum(sum(OutFRs))-sum(OutFRs),[size(OutFRs,1) 1]);
%     r = rand(size(OutFRs));
%     rz = [r./sum(r)].*(bias.*size(OutFRs,1));
%     OutFRs = OutFRs+rz;
%     
%     
%     bias = repmat(mean(mean(InFRs))-mean(InFRs),[size(InFRs,1) 1]);
%     r = rand(size(InFRs));
%     rz = [r./sum(r)].*(bias.*size(InFRs,1));
%     InFRs = InFRs+rz;
%     t = toc
%     
%     clear r rz bias
    
%     slower but fields dont go under 0
    OutFRs = OutFRsave;
    InFRs = InFRsave;
%     tic
    while max(abs(mean(OutFRs)-mean(mean(OutFRs))))>.00001
          bias = repmat(mean(OutFRs)-mean(mean(OutFRs)),[size(OutFRs,1) 1]);
          OutFRs = OutFRs+((-bias).*rand(size(OutFRs)));
          OutFRs(OutFRs<0) = .00001;
    end

    while max(abs(mean(InFRs)-mean(mean(InFRs))))>.00001
          bias = repmat(mean(InFRs)-mean(mean(InFRs)),[size(InFRs,1) 1]);
          InFRs = InFRs+((-bias).*rand(size(InFRs)));
          InFRs(InFRs<0) = .00001;
    end
%     t2 = toc

    
    if toplot && ishuff==1 %length(Cell_Number)>100 
       
        figure; hold on;
        [~,m] = max(OutFRsave');
        [~,o1] = sort(m);
        subplot(3,3,1); hold on
        imagesc(OutFRsave(o1,:)./max(OutFRsave(o1,:),[],2))
        xlabel('Position'); ylabel('Cells')
        axis tight; colormap hot
        
        subplot(3,3,2); hold on
        imagesc(OutFRs(o1,:)./max(OutFRs(o1,:),[],2))
        xlabel('Position (Flattened Field Distribution)'); ylabel('Cells (in same order)')
        axis tight;
        
        subplot(3,3,3); hold on
        oind = o1(randperm(length(o1)));
        imagesc(OutFRsave(oind,:)./max(OutFRsave(oind,:),[],2))
        xlabel('Position (Cell ID Shuffle)'); ylabel('Cells (in same order)')
        axis tight;
        
        subplot(3,3,4); hold on;
        plot(sum(OutFRs),'r');   plot(sum(OutFRsave),'k');   
        axis tight;        
        xlabel('Position'); ylabel('Sum Place Fields')
        
        subplot(3,3,5); hold on;
        plot(mean(OutFRs),'r'); plot(mean(OutFRsave),'k')
        axis tight;
        
        xlabel('Position'); ylabel('Mean Place Fields')
        for icell = 1:3
            subplot(3,3,6+icell); hold on;        
            plot(OutFRsave(icell,:),'k'); plot(OutFRs(icell,:),'r'); axis tight
            xlabel('Position'); ylabel(['Cell ' num2str(icell) ' FR'])
        end
        
        suptitle([iident(1:end-4) ', ' num2str(length(Cell_Number)) ' Cells'])
        set(gcf,'Position',[  689    89   782   707])
        set(gcf, 'Renderer', 'painters');        
        if final            
            helper_savefig([dirs.figdir '\Final\' iident(1:end-4) '_ExampleFlat'])
        end
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' iident(1:end-4) '_ExampleFlat'])
    end
    term2=zeros(Number_Of_Bins,TimeBins);
    term4=term2;

    if Number_Of_Bins>TimeBins    
        for TBin=1:TimeBins
            term2(:,TBin)=prod(OutFRs.^binspike(:,TBin),1);    
            term4(:,TBin)=prod(InFRs.^binspike(:,TBin),1);
        end
    else
        for PosBin=1:Number_Of_Bins   
            term2(PosBin,:)=prod((OutFRs(:,PosBin)*ones(1,TimeBins)).^binspike,1);    
            term4(PosBin,:)=prod((InFRs(:,PosBin)*ones(1,TimeBins)).^binspike,1);
        end
    end

    term3=exp(-EstBin*sum(OutFRs,1)')*ones(1,TimeBins);

    term5=exp(-EstBin*sum(InFRs,1)')*ones(1,TimeBins);

    Mat=term2.*term3+term4.*term5;

    OutMat1=term2.*term3;
    OutMatrix=(OutMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';

    InMat1=term4.*term5;
    InMatrix=(InMat1./(ones(size(Mat,1),1)*sum(Mat,1)))';    
    
    clear term2 term3 term4 term5 OutMat1 InMat1
    
    nospks = sum(binspike>0)==0;
    OutMatrix(nospks,:) = NaN;
    InMatrix(nospks,:) = NaN;
    OutMatrixShuff(:,:,ishuff) = OutMatrix;
    InMatrixShuff(:,:,ishuff) = InMatrix;
    clear OutMatrix InMatrix
    if rem(ishuff,10)==0
        disp(['Done ' num2str(ishuff) ' shuffles'])
    end
end





% save(iident,'OutMatrixShuff','InMatrixShuff','-append')
