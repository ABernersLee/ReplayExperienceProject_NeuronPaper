function [r,p] = plotsubplot_TF_FR_SN_DayDay(f,col,dat,label,yl,plotcol,datt,isNovel,isDValue,dattind,TrackType,allassig,sigonly,ratday,zscored)

set(groot,'CurrentFigure',f)
gca=subplot(dat.sz(1),dat.sz(2),plotcol-1);
hold on;

r = ones(2,1); p = r;


if isNovel==1
    ind1 = dattind(:,end-1)==1;
elseif isNovel==0
    ind1 = dattind(:,end-1)==0; 
else
    ind1 = true(size(dattind,1),1);
end

if isDValue==1
    ind3 = dattind(:,end)~=0;
elseif isDValue==0
    ind3 = dattind(:,end)==0;
elseif isDValue==2
    ind3 = dattind(:,end)==1;
elseif isDValue==3
    ind3 = dattind(:,end)==-1;
else
    ind3 = true(size(dattind,1),1);
end


 if TrackType>0
    ind2 = dattind(:,end-2)==TrackType;
else
    ind2 = true(size(dattind,1),1);
 end


ind4 = ismember(dattind(:,end-4),find(~isnan(ratday(:,1))));

datt = datt(ind1 & ind2 & ind3 & ind4,:,:,:);
dattind2 = dattind(ind1&ind2 &ind3 & ind4,:);
if ~isempty(datt)

    for iSig = 1:2
        
        if allassig && iSig ==1
            datt1 = cat(1,datt(:,:,:,1),datt(:,:,:,2));
        elseif allassig && iSig ==2
            continue
        else
            datt1 = datt(:,:,:,iSig);
        end
        
        if sigonly && iSig==2
            continue
        end
        
        if zscored==1
            for irat = 1:max(ratday(:,2))                
                ds = ratday(ratday(:,2)==irat,1);
                if ~isempty(ds)
                    rind = ismember(dattind2(:,end-4),find(ratday(:,2)==irat));
                    datt1(rind,:,col) = nanzscore(datt1(rind,:,col));
                end
            end
        end
        A=NaN(max(ratday(dattind(:,end-4),1)),2);
        for j=1:max(ratday(dattind(:,end-4),1))
            da = datt1(ratday(dattind2(:,end-4),1)==j,:,col);
            A(j,1)=nanmean(da(:));
            A(j,2)=nanstd(da(:))/sqrt(sum(~isnan(da(:))));        
        end    
        A(sum(~isnan(A),2)==0,:) = [];


        if iSig == 1
            errorbar(A(:,1),A(:,2),'k','LineWidth',2);
        elseif iSig == 2
            errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
        end


        FR=datt1(:,:,col); 
        A=[dattind2(:,end-3)]*ones(1,size(FR,2));
        t=find(~isnan(FR));

        [r(setdiff(1:2,iSig)),p(setdiff(1:2,iSig))]=corr(A(t),FR(t),'rows','complete');  
    end

    xlabel('Day')
    set(gca,'XTick',[1:max(ratday(dattind(:,end-4),1))+1])
    ylabel(label)
    xlim([0 max(ratday(dattind(:,end-4)))+1])
    if sum(yl==0)~=2
        ylim(yl)
    end
    if ~sigonly
%     if p(2)<.05 && p(1)>=.05
%         title({['\color{gray}r=' num2str(round(r(1),2,'significant')) ...
%                  ', p=' num2str(round(p(1),2,'significant'))  ';     \color{red}r=' num2str(round(r(2),2,'significant')) ...
%                  ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
%     elseif p(1)<.05 && p(2)>=.05
%         title({['\color{red}r=' num2str(round(r(1),2,'significant')) ...
%                  ', p=' num2str(round(p(1),2,'significant'))  ';     \color{black}r=' num2str(round(r(2),2,'significant')) ...
%                  ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
%     elseif p(1)<.05 && p(2)<.05
%         title({['\color{red}r=' num2str(round(r(1),2,'significant')) ...
%                  ', p=' num2str(round(p(1),2,'significant'))  ';     \color{red}r=' num2str(round(r(2),2,'significant')) ...
%                  ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);    
%     else
        title({['\color{gray}r=' num2str(round(r(1),2,'significant')) ...
                 ', p=' num2str(round(p(1),2,'significant'))  ';     \color{black}r=' num2str(round(r(2),2,'significant')) ...
                 ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
%     end
    elseif sigonly
        title({['\color{black}r=' num2str(round(r(2),2,'significant')) ...
             ', p=' num2str(round(p(2),1,'significant'))]},'FontSize',8);
    end
end