function [r,p] = plotsubplot_TFAlldat_phase_DayDay(f,col,dat,label,yl,plotcol,isNovel,isDValue,TrackType,sigonly,ratday)
% ratdayplay = ~isnan(ratday);
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
        ind1 = AllDat1(:,end-1)==0; % & AllDat(:,end)==0;
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
    
%     if Run>0
%         ind4 = AllDat1(:,end-3)==Run;
%     else
%         ind4 = true(size(AllDat1,1),1);
%     end
%     ind4 = ratdayplay(AllDat1(:,end-4));
    ind4 = ismember(AllDat1(:,end-4),find(~isnan(ratday(:,1))));

    AllDat = AllDat1(ind1 & ind2 & ind3 & ind4,:);
    if isempty(AllDat)
        continue
    end
    A=NaN(max(ratday(AllDat(:,end-4),1)),2);
    for j=1:max(ratday(AllDat(:,end-4),1))
        t=ratday(AllDat(:,end-4),1)==j;
        if sum(t)>0
            [x,~,y1,y2] = abl_circ_mean_ci(AllDat(t,col),[],[],[],[],0);
            A(j,1)= x; %nanmean(AllDat(t,col)); %        
            A(j,2)= rad2deg(min(abs(x-y1),abs(x-y2))); %/sqrt(length(t)); %nanstd(AllDat(t,col))/sqrt(length(t));       
        end
    end    
    A(sum(~isnan(A),2)==0,:) = [];
%     A(:,1) = rad2deg(unwrap(A(:,1)+pi));
    A(:,1) = rad2deg(A(:,1));
    
    if iSig
        errorbar(A(:,1),A(:,2),'k','LineWidth',2);
    else
        errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
    end
    
    xlabel('Day')
    set(gca,'XTick',[1:max(ratday(AllDat(:,end-4),1))])
    [r(iSig+1),p(iSig+1)]=corr([1:max(ratday(AllDat(:,end-4),1))]',A(:,1),'rows','complete');       
end

ylabel(label)
xlim([0 max(ratday(AllDat(:,end-4)))+1])
if sum(yl==0)~=2
    ylim(yl)
end
if ~sigonly
% if p(2)<.05 && p(1)>=.05
%     title({['\color{gray}r=' num2str(round(r(1),2,'significant')) ...
%              ', p=' num2str(round(p(1),2,'significant'))  ';     \color{red}r=' num2str(round(r(2),2,'significant')) ...
%              ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
% elseif p(1)<.05 && p(2)>=.05
%     title({['\color{red}r=' num2str(round(r(1),2,'significant')) ...
%              ', p=' num2str(round(p(1),2,'significant'))  ';     \color{black}r=' num2str(round(r(2),2,'significant')) ...
%              ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
% elseif p(1)<.05 && p(2)<.05
%     title({['\color{red}r=' num2str(round(r(1),2,'significant')) ...
%              ', p=' num2str(round(p(1),2,'significant'))  ';     \color{red}r=' num2str(round(r(2),2,'significant')) ...
%              ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);    
% else
    title({['\color{gray}r=' num2str(round(r(1),2,'significant')) ...
             ', p=' num2str(round(p(1),2,'significant'))  ';     \color{black}r=' num2str(round(r(2),2,'significant')) ...
             ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
% end
elseif sigonly
    title({['\color{black}r=' num2str(round(r(2),2,'significant')) ...
             ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
end