function [r,p] = plotsubplot_TFAlldat_phase(f,col,dat,label,yl,plotcol,isNovel,isDValue,TrackType,Run,sigonly)

if isempty(plotcol)
    plotcol = col;
end
set(groot,'CurrentFigure',f);
gca=subplot(dat.sz(1),dat.sz(2),plotcol-1);
hold on;

r = ones(2,1); p = r;

for iSig = 0:1
    if sigonly && iSig==0
        continue
    end
    A=NaN(dat.Lap,2);
    
    %significant replays or non significant replays
    if iSig==1
        AllDat1 = dat.AllSig;
    elseif iSig==0
        AllDat1 = dat.AllnonSig;
    else
        AllDat1 = [dat.AllSig; dat.AllnonSig];
    end
    
    %whether the session is novel, linear/not, or has a reward change
    if isNovel==1 
        ind1 = AllDat1(:,end-1)==1;
    elseif isNovel==0
        ind1 = AllDat1(:,end-1)==0; 
    else
       ind1 = true(size(AllDat1,1),1);
    end        
    
    if TrackType>0
        ind2 = AllDat1(:,end-2)==TrackType;
    else
        ind2 = true(size(AllDat1,1),1);
    end
    
    if isDValue==1
       ind3 = AllDat1(:,end)~=0;
    elseif isDValue==0
       ind3 = AllDat1(:,end)==0;
    elseif isDValue==2
        ind3 = AllDat1(:,end)==1;
    elseif isDValue==3
        ind3 = AllDat1(:,end)==-1;
    else
       ind3 = true(size(AllDat1,1),1);
    end
    
    if Run>0
        ind4 = AllDat1(:,end-3)==Run;
    else
        ind4 = true(size(AllDat1,1),1);
    end

    AllDat = AllDat1(ind1 & ind2 & ind3 & ind4,:);
    if isempty(AllDat)
        continue
    end

    for j=1:dat.Lapskip:dat.Lap    
        t=AllDat(:,1)==j | AllDat(:,1)==j+dat.Lapskip-1;
        if sum(t)>0
            [x,~,y1,y2] = abl_circ_mean_ci(AllDat(t,col),[],[],[],[],0);
            A(j,1)= x; 
            A(j,2)= rad2deg(min(abs(x-y1),abs(x-y2)));    
        end
    end    
    A(sum(~isnan(A),2)==0,:) = [];
    A(:,1) = rad2deg(A(:,1));
    
    if iSig
        errorbar(A(:,1),A(:,2),'k','LineWidth',2);
    else
        errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
    end
    
    xlabel('Lap number')
    set(gca,'XTick',[0:5:dat.Lap])
    t=find(AllDat(:,1)>0 & AllDat(:,1)<=dat.Lap);
    [r(iSig+1),p(iSig+1)]=corr([1:size(A,1)]',A(:,1),'rows','complete');         
end

ylabel(label)
xlim([0 (dat.Lap/dat.Lapskip)+1])
if sum(yl==0)~=2
    ylim(yl)
end
if ~sigonly

    title({['\color{gray}r=' num2str(round(r(1),2,'significant')) ...
             ', p=' num2str(round(p(1),2,'significant'))  ';     \color{black}r=' num2str(round(r(2),2,'significant')) ...
             ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
elseif sigonly
        title({['\color{black}r=' num2str(round(r(2),2,'significant')) ...
                 ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
end
