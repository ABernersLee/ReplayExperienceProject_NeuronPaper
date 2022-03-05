function [r,p] = plotsubplot_TFhoverjump_SN_DayDay(f,col,dat,label,yl,plotcol,datt,isNovel,isDValue,TrackType,sigonly,ratday,zscored)

set(groot,'CurrentFigure',f)
gca=subplot(dat.sz(1),dat.sz(2),plotcol-1);
hold on;
r = ones(2,1); p = r;

for iSig = 0:1    
    if sigonly && iSig==0
        continue
    end
    datt1 = datt(datt(:,2)==iSig,[1 3:end]);
    

    
   
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
    
    ind4 = ismember(datt1(:,end-4),find(~isnan(ratday(:,1))));

    datt2 = datt1(ind1 & ind2 & ind3 & ind4,:);
    
    if isempty(datt2)
        continue
    end
    
    if zscored==1
        for irat = 1:max(ratday(:,2))
            ds = ratday(ratday(:,2)==irat,1);
            if ~isempty(ds)
                rind = ismember(datt2(:,end-4),find(ratday(:,2)==irat));

                datt2(rind,col) = nanzscore(datt2(rind,col));
            end
        end
    end
    A=NaN(max(ratday(datt2(:,end-4),1)),2);
    for j=1:max(ratday(datt2(:,end-4),1))
        t=find(ratday(datt2(:,end-4),1)==j);
        A(j,1)=mean(datt2(t,col));
        A(j,2)=std(datt2(t,col))/sqrt(length(t));        
    end    
    A(sum(~isnan(A),2)==0,:) = [];

    if iSig
        errorbar(A(:,1),A(:,2),'k','LineWidth',2);
    else
        errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
    end
    xlabel('Day')
    set(gca,'XTick',[1:max(ratday(datt2(:,end-4),1))+1])
    [r(iSig+1),p(iSig+1)]=corr(ratday(datt2(:,end-4)),datt2(:,col),'rows','complete');          

end
ylabel(label)
xlim([0 max(ratday(datt2(:,end-4)))+1])
    
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
 

