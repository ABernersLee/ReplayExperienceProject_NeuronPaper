function [r,p] = plotsubplot_TFhoverjump_SN(f,col,dat,label,yl,plotcol,datt,isNovel,isDValue,TrackType,Run,sigonly)

set(groot,'CurrentFigure',f)
gca=subplot(dat.sz(1),dat.sz(2),plotcol-1);
hold on;
r = ones(2,1); p = r;

for iSig = 0:1    
    if sigonly && iSig==0
        continue
    end
    datt1 = datt(datt(:,2)==iSig,[1 3:end]);
    

    A=NaN(dat.Lap,2);
   
    if isNovel==1
        ind1 = datt1(:,end-1)==1;
    elseif isNovel==0
        ind1 = datt1(:,end-1)==0;
    else
        ind1 = true(size(datt1,1),1);
    end
    
    if TrackType>0
        ind2 = datt1(:,end-2)==TrackType;
    else
        ind2 = true(size(datt1,1),1);
    end
    
    if isDValue==1
       ind3 = datt1(:,end)~=0;
    elseif isDValue==0
       ind3 = datt1(:,end)==0;
    elseif isDValue==2
         ind3 = datt1(:,end)==1;
    elseif isDValue==3
        ind3 = datt1(:,end)==-1;
    else
       ind3 = true(size(datt1,1),1);
    end
    
    if Run>0
        ind4 = datt1(:,end-3)==Run;
    else
        ind4 = true(size(datt1,1),1);
    end

    datt2 = datt1(ind1 & ind2 & ind3 & ind4,:);
    
    if isempty(datt2)
        continue
    end
    
    for j=1:dat.Lapskip:dat.Lap    
        t=find(datt2(:,1)==j | datt2(:,1)==j+dat.Lapskip-1);
        A(j,1)=nanmean(datt2(t,col));
        A(j,2)=nanstd(datt2(t,col))/sqrt(sum(~isnan(datt2(t,col))));        
    end    
    A(sum(~isnan(A),2)==0,:) = [];

    if iSig
        errorbar(A(:,1),A(:,2),'k','LineWidth',2);
    else
        errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
    end
    xlabel('Lap number')
    set(gca,'XTick',[0:5:dat.Lap])
    t=find(datt2(:,1)>0 & datt2(:,1)<=dat.Lap);
    if dat.Lapskip>1 % comparing for Ting's old effect
        [r(iSig+1),p(iSig+1)]=corr(ceil(datt2(t,1)./dat.Lapskip),datt2(t,col),'rows','complete');   
    else
        [r(iSig+1),p(iSig+1)]=corr(datt2(t,1),datt2(t,col),'rows','complete');
    end

    

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
    

