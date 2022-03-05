function Hovers_Moving_Shuffle(basedir,plot_paper_figures,plot_other_figures)

%%%%%
%%%%%
%%%%% This uses a cell ID shuffle to test whether this distribution of
%%%%% place fields would always produce bumpy posteriors or whether it
%%%%% depends on the specific firing of the HP neurons during replay.
%%%%%
%%%%%


%%%%% set up parameters, labels, and directories

wc_cutoff = 0;
jd_cutoff = 1;
coveragecutoff = 0;
load([basedir 'dirs_linear_lapbylap_addedIN_wholesession.mat'],'dirs')
cd(dirs.spikedatadir)
if ~isfolder([dirs.figdir '/HoversMovingShuffle/'])
    mkdir([dirs.figdir '/HoversMovingShuffle/'])
end
addlab = ['_wholefields_wc' num2str(wc_cutoff) '_jd' ...
    num2str(jd_cutoff) '_coverage' num2str(coveragecutoff) '_CellIDShuffle'];
List=dir('*.mat');
numses = size(List,1);
numperm = 100;    
r_all = []; cellnum = [];

%%%%% go through sessions

for ListNo= 1:numses 

    load(List(ListNo).name,...
        'Index','InMatrix','OutMatrix',...
        'Spike','CandSeq','MidTime',...
        'hp_cells','hpinterneurons','-mat') 
    
    cn = hp_cells(~ismember(hp_cells,hpinterneurons));
    cellnum = cat(1,cellnum,length(cn));
     [OutMatrixShuff,InMatrixShuff] = decode_spikedensity_events_shuffle(List(ListNo).name,numperm);
    
    [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);
    if sum(I==0)>0
        I(I==0) = NaN;
    end
    I(I>length(MidTime)) = length(MidTime);
    
    

    Mat=OutMatrix'+InMatrix';
    Mat2 = permute(OutMatrixShuff,[2 1 3])+permute(InMatrixShuff,[2 1 3]);
    clear OutMatrixShuff InMatrixShuff
    D=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
    D2=squeeze(sum([1:size(Mat2,1)]'*ones(1,size(Mat2,2)).*Mat2));
    
    isN = squeeze(isnan(Mat2(1,:,:)));
    D2(isN) = NaN;
    D(isN) = NaN;
 
    t = 1:length(CandSeq);       
    
    hovpos = []; hovpos2 = [];
    for it = 1:length(t)
       pos = D(Index(t(it))*4+1:Index(t(it)+1)*4-3);              
       DiffLoc=abs(diff(pos))';     
       Hov1 = [false(1,size(DiffLoc,2));DiffLoc<2];
       Hov2 = [DiffLoc<2;false(1,size(DiffLoc,2))];
       Hov3 = Hov1 | Hov2;
%        figure; plot(pos,'.-'); hold on; plot(find(Hov3(:,1)),pos(Hov3(:,1)),'r*')        
       pos(~Hov3) = NaN; 
       hovpos = [hovpos; pos' ones(length(pos),1)*I(t(it)) ones(length(pos),1)*t(it)];
       
       pos2 = D2(Index(t(it))*4+1:Index(t(it)+1)*4-3,:);       
       DiffLoc=abs(diff(pos2));     
       Hov1 = [false(1,size(DiffLoc,2));DiffLoc<2];
       Hov2 = [DiffLoc<2;false(1,size(DiffLoc,2))];
       Hov3 = Hov1 | Hov2;
%        figure; plot(pos2(:,1),'.-'); hold on; plot(find(Hov3(:,1)),pos2(Hov3(:,1),1),'r*')        
       pos2(~Hov2) = NaN; 
       hovpos2 = cat(1,hovpos2,cat(3,pos2,ones(size(pos2))*(I(t(it))),ones(size(pos2))*t(it)));
    end
    hovpos2 = permute(hovpos2,[1 3 2]);
    
    
    hhh = histc(hovpos(:,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))));
    hhh2 = NaN(length(hhh),size(hovpos2,3));
    r2 = NaN(size(hovpos2,3),1);
    for ishuff = 1:size(hovpos2,3)
            hhh1 = histc(hovpos2(:,1,ishuff),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))));
            r2(ishuff) = corr(hhh,hhh1);
            hhh2(:,ishuff) = hhh1;
    end
    
    clear hhh1
     
    %%%%% plot example session
    if plot_other_figures || (plot_paper_figures && ListNo==35)
        figure; hold on; 
        set(gcf,'Position',[  2098        -210         438         726])
        subplot(3,1,1); hold on
        plot(100*(hhh./sum(hhh)),'k-'); 
        plot(100*(hhh2(:,1)./sum(hhh2(:,1))),'r-'); 
        plot(100*(hhh2(:,2)./sum(hhh2(:,2))),'b-'); 
        title([List(ListNo).name(1:end-4)])
        set(gca,'FontSize',18,'FontName','Arial')
        axis tight
        subplot(3,1,2); hold on        
        histogram(r2,'FaceColor','w')
        title(['Mean Raw R = ' num2str(round(mean(r2),2))])
        set(gca,'FontSize',18,'FontName','Arial')
        
        set(gcf,'renderer','Painters')
        if (plot_paper_figures && ListNo==35)            
            helper_savefig([dirs.figdir 'Final\' List(ListNo).name(1:end-4) '_corr' addlab])
        end
        helper_saveandclosefig([dirs.figdir 'HoversMovingShuffle\' List(ListNo).name(1:end-4) '_corr' addlab])
    end
    
    
   r_all = cat(2,r_all,r2);
   
end

%%%%% only sessions with over 60 HP cells

if plot_paper_figures
    r_all2 = r_all(:,cellnum>60);

    figure; hold on
    histogram(r_all2(:),'FaceColor','w')
    title(['Over 60 Cells, Mean Raw R = ' num2str(round(mean(r_all2(:)),2))])
    set(gca,'FontSize',18,'FontName','Arial')
    yl = get(gca,'ylim');
    plot([mean(mean(r_all2)) mean(mean(r_all2))],yl,'r-','LineWidth',2)
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([dirs.figdir 'HoversMovingShuffle\corr_60cell_all' addlab])

    figure; hold on
    histogram(mean(r_all2),5,'FaceColor','w')
    title(['Over 60 Cells, Mean Raw R = ' num2str(round(mean(mean(r_all2)),2))])
    set(gca,'FontSize',18,'FontName','Arial')
    yl = get(gca,'ylim');
    plot([mean(mean(r_all2)) mean(mean(r_all2))],yl,'r-','LineWidth',2)
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([dirs.figdir 'HoversMovingShuffle\corr_60cell_mean' addlab])
end

    