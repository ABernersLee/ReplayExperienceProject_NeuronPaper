function [r,p] = plotsubplot_TF_FR_SN(f,col,dat,label,yl,plotcol,datt,isNovel,isDValue,dattind,TrackType,Run,allassig,sigonly)

set(groot,'CurrentFigure',f)
gca=subplot(dat.sz(1),dat.sz(2),plotcol-1);
hold on;

r = ones(2,1); p = r;

if ~isempty(dattind)
    
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

 if Run>0
    ind4 =dattind(:,end-3)==Run;
else
    ind4 = true(size(dattind,1),1);
end


datt = datt(ind1 & ind2 & ind3 & ind4,:,:,:);
if ~isempty(datt)

    for iSig = 1:2
        A=NaN(dat.Lap*dat.Lapskip,2);
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

        for j=1:dat.Lapskip:min(dat.Lap,size(datt1,2))
            da = datt1(:,j:j+dat.Lapskip-1,col);
            A(j,1)=nanmean(da(:));
            A(j,2)=nanstd(da(:))/sqrt(sum(~isnan(da(:))));        
        end    
        A(sum(~isnan(A),2)==0,:) = [];


        if iSig == 1
            errorbar(A(:,1),A(:,2),'k','LineWidth',2);
        elseif iSig == 2
            errorbar(A(:,1),A(:,2),'Color',[.5 .5 .5],'LineWidth',2);
        end


        FR=datt1(:,1:dat.Lap,col);
        A=ones(size(FR,1),1)*[1:dat.Lap];
        t=find(~isnan(FR));
        
        if size(A(t),1)==1
            [r(setdiff(1:2,iSig)),p(setdiff(1:2,iSig))]=corr(A(t)',FR(t)','rows','complete');  
        else
            [r(setdiff(1:2,iSig)),p(setdiff(1:2,iSig))]=corr(A(t),FR(t),'rows','complete'); 
        end
    end

    xlabel('Lap number')
    set(gca,'XTick',[0:5:dat.Lap])
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
end
end