function [r,p] = plotsubplot_TFAlldat_SN(f,col,dat,label,yl,plotcol,isNovel,isDValue,TrackType,Run,sigonly,finalfigfolder)

if isempty(plotcol)
    plotcol = col;
end
set(groot,'CurrentFigure',f);
gca=subplot(dat.sz(1),dat.sz(2),plotcol-1);
hold on;

r = ones(2,1); p = r;
n = NaN(2,1);
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
        t=find(AllDat(:,1)==j | AllDat(:,1)==j+dat.Lapskip-1);
        A(j,1)=nanmean(AllDat(t,col));
        A(j,2)=nanstd(AllDat(t,col))/sqrt(sum(~isnan(AllDat(t,col))));        
    end    
    A(sum(~isnan(A),2)==0,:) = [];
    
    if iSig
        errorbar(A(:,1),A(:,2),'k','LineWidth',2);
    else
        errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
    end
    
    xlabel('Lap number')
    set(gca,'XTick',[0:5:dat.Lap])
    t=find(AllDat(:,1)>0 & AllDat(:,1)<=dat.Lap);
    [r(iSig+1),p(iSig+1)]=corr(AllDat(t,1),AllDat(t,col),'rows','complete');       
    if iSig==0
        AllDat0 = AllDat;
    end
    n(iSig+1) = length(t);
end
if ~sigonly
adat = [AllDat(AllDat(:,1)>0 & AllDat(:,1)<=dat.Lap,col); AllDat0(AllDat0(:,1)>0 & AllDat0(:,1)<=dat.Lap,col)];
dat11 = [AllDat(AllDat(:,1)>0 & AllDat(:,1)<=dat.Lap,1); AllDat0(AllDat0(:,1)>0 & AllDat0(:,1)<=dat.Lap,1)];
dat22 = [ones(size(AllDat(AllDat(:,1)>0 & AllDat(:,1)<=dat.Lap,col))); 2*ones(size(AllDat0(AllDat0(:,1)>0 & AllDat0(:,1)<=dat.Lap,col)))];
[pp,tbl] = anovan(adat,[dat11 dat22],'continuous',1,'model','interaction','display','off');
end
ylabel(label)
xlim([0 (dat.Lap/dat.Lapskip)+1])
if sum(yl==0)~=2
    ylim(yl)
end
if ~sigonly

        title({['\color{gray}n=' num2str(n(1)) ', r=' num2str(round(r(1),2,'significant')) ...
                 ', p=' num2str(round(p(1),2,'significant'))  ';\color{black}n=' num2str(n(2)) ', r=' num2str(round(r(2),2,'significant')) ...
                 ', p=' num2str(round(p(2),1,'significant')) ', i=' num2str(round(pp(3),1,'significant'))]},'FontSize',8);
         
               if col <4
    figure; hold on;
    plab = {'Group';'Lap';'Interaction'};
    for ip = 1:3
        text(0,ip/4,[plab{ip} ': F = ' num2str(tbl{ip+1,6}) ' df = ' num2str(tbl{ip+1,3}) ',' num2str(length(adat)-tbl{ip+1,3}) ' , p = ' num2str(pp(ip))])
    end
    set(gcf,'Position',[1919         133         861         403])
    axis off
    helper_saveandclosefig([finalfigfolder 'Replays_SigNonSig_AllGroupsStats_ANOVA_Laps_' label])
               end

elseif sigonly
        title({['\color{black}n=' num2str(n(2)) ', r=' num2str(round(r(2),2,'significant')) ...
                 ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
end