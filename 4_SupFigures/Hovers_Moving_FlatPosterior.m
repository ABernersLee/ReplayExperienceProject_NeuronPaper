function Hovers_Moving_FlatPosterior(basedir,plot_paper_figures,plot_other_figures)

%%%%%
%%%%% 
%%%%% This makes a flat posterior in order to compare with a real posterior
%%%%% to test whether the bumpiness of the replays are only a function of a
%%%%% bumpy posterior. I do this by adding biased noise to each neruon such
%%%%% that each neuron's field shouldn't matter much but overall the
%%%%% distribution is flat. This is for supplemental figure 6. This also
%%%%% generates figures that aren't in the paper (plot_other_figures) and
%%%%% example sessions
%%%%%
%%%%%


%%%%% set up parameters and directories

wc_cutoff = 0;
jd_cutoff = 1;
coveragecutoff = 0; 
load([basedir 'dirs_linear_lapbylap_addedIN_wholesession.mat'],'dirs')
cd(dirs.spikedatadir)
if ~isfolder([dirs.figdir '/HoversMovingShuffle/'])
    mkdir([dirs.figdir '/HoversMovingShuffle/'])
end
addlab = ['_wholefields_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_coverage' num2str(coveragecutoff) '_Flat_FamandNov'];
shuffr = []; realr = []; 
List=dir('*.mat');
numses = size(List,1);
numshuff = 100;
r_all = []; cellnum = []; r_running_all = []; r_running_all_sm = [];


for ListNo= 1:numses 

    %%%%% get data out
    load(List(ListNo).name,'InFR',...
        'replayparams','Index','InMatrix','OutMatrix',...
        'Spike','OutFR','CandSeq','MidTime','params',...
        'hp_cells','hpinterneurons','-mat')
    
    cn = hp_cells(~ismember(hp_cells,hpinterneurons));
    cellnum = cat(1,cellnum,length(cn));
    [OutMatrixFlat,InMatrixFlat] = decode_spikedensity_events_flatposterior(dirs,List(ListNo).name,...
        numshuff,plot_other_figures || (plot_paper_figures && ListNo==35),plot_paper_figures && ListNo==35);
    
    [OutMatrixRunning,InMatrixRunning,binspike] = decode_spikedensity_events_running(List(ListNo).name);
    MatRunning = OutMatrixRunning'+InMatrixRunning';
    Drunning=sum([1:size(MatRunning,1)]'*ones(1,size(MatRunning,2)).*MatRunning);
    Drunning(sum(binspike)<2)=NaN; 
    
    [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);
    if sum(I==0)>0
        I(I==0) = NaN;
    end
    I(I>length(MidTime)) = length(MidTime);
        
    Mat=OutMatrix'+InMatrix';
    Mat2 = permute(OutMatrixFlat,[2 1 3])+permute(InMatrixFlat,[2 1 3]);
    
    D=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
    clear OutMatrix InMatrix 
    D2=squeeze(sum([1:size(Mat2,1)]'*ones(1,size(Mat2,2)).*Mat2));    
    isN = squeeze(isnan(Mat2(1,:,:)));
    D2(isN) = NaN;
    D(isN) = NaN;
    clear OutMatrixFlat InMatrixFlat 
    
     
    %%%%% get hovers out
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
    

    hhhR = histc(Drunning',floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))));   
    hhh = histc(hovpos(:,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))));    
    r3 = corr(hhh,hhhR);   
    
    hhh22 = smoothdata(hhh,'movmean',5);        
    hhhR2 = smoothdata(hhhR,'movmean',5); 
    r32 = corr(hhh22,hhhR2);  
    
    %%%%% find peaks
    hhh2 = NaN(length(hhh),size(hovpos2,3));
    r2 = NaN(size(hovpos2,3),1); hh = r2; locs1 = cell(2,1);
    for ishuff = 1:size(hovpos2,3)
            hhh1 = histc(hovpos2(:,1,ishuff),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))));
            r2(ishuff) = corr(hhh,hhh1);
            hhh2(:,ishuff) = hhh1;
            [~,locs11,~,~] = findpeaks((hhh1./sum(hhh1))*100,'MinPeakProminence',.4,'MinPeakDistance',size(hhh1,1)*.10);     
            hh(ishuff) = length(locs11);
            if ishuff<3
                locs1{ishuff} = locs11;
            end
    end
    [~,locs2,~,~] = findpeaks((hhh./sum(hhh))*100,'MinPeakProminence',.4,'MinPeakDistance',size(hhh,1)*.10);
    h = length(locs2); 
 
    if plot_other_figures || (plot_paper_figures && ListNo==35)
        
        figure; hold on; 
        set(gcf,'Position',[    2098          71         438         445])
        subplot(2,1,1); hold on
        plot(100*(hhh./sum(hhh)),'k-'); 
        plot(100*(hhh2(:,1)./sum(hhh2(:,1))),'r-'); 
        plot(100*(hhh2(:,2)./sum(hhh2(:,2))),'b-');         
        title([List(ListNo).name(1:end-4)])
        set(gca,'FontSize',18,'FontName','Arial')
        axis tight        
        subplot(2,1,2); hold on        
        histogram(r2,'FaceColor','w')
        title(['Mean Raw R = ' num2str(round(mean(r2),2))])
        set(gca,'FontSize',18,'FontName','Arial')        
        set(gcf,'renderer','Painters')
        if (plot_paper_figures && ListNo==35)                 
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_corr' addlab])
        end
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' List(ListNo).name(1:end-4) '_corr' addlab])
        
        figure; hold on; 
        set(gcf,'Position',[    2098          71         438         445])
        subplot(2,1,1); hold on
        plot(100*(hhh./sum(hhh)),'k-');       
        plot(100*(hhhR./sum(hhhR)),'g-'); 
        set(gca,'FontSize',18,'FontName','Arial')
        axis tight             
        title(['Mean Raw R = ' num2str(round(mean(r3),2))])
        subplot(2,1,2); hold on        
        title([List(ListNo).name(1:end-4)])
        set(gca,'FontSize',18,'FontName','Arial')        
        set(gcf,'renderer','Painters')
        if (plot_paper_figures && ListNo==35)                 
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_RunningHist' addlab])
        end
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' List(ListNo).name(1:end-4) '_RunningHist' addlab])
        
        figure; hold on; 
        set(gcf,'Position',[    2098          71         438         445])
        subplot(2,1,1); hold on         
        plot(100*(hhh./sum(hhh)),'k-');          
        plot(100*(hhh22./sum(hhh22)),'k-','LineWidth',3);            
        plot(100*(hhhR./sum(hhhR)),'g-');             
        plot(100*(hhhR2./sum(hhhR2)),'g-','LineWidth',3);
        set(gca,'FontSize',18,'FontName','Arial')
        axis tight                      
        title(['Mean Raw R = ' num2str(round(mean(r32),2))])
        subplot(2,1,2); hold on        
        title([List(ListNo).name(1:end-4)])
        set(gca,'FontSize',18,'FontName','Arial')        
        set(gcf,'renderer','Painters')
        if (plot_paper_figures && ListNo==35)                 
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_RunningHistSmooth' addlab])
        end
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' List(ListNo).name(1:end-4) '_RunningHistSmooth' addlab])
    end
    
    if plot_other_figures || (plot_paper_figures && ListNo==35)
        
        figure; hold on; 
        set(gcf,'Position',[  2098        -210         438         726])
        subplot(3,1,1); hold on
        plot(100*(hhh2(:,1)./sum(hhh2(:,1))),'r-'); plot(locs1{1},100*(hhh2(locs1{1},1)./sum(hhh2(:,1))),'r*')
        plot(100*(hhh2(:,2)./sum(hhh2(:,2))),'b-'); plot(locs1{2},100*(hhh2(locs1{2},2)./sum(hhh2(:,2))),'b*')
        plot(100*(hhh./sum(hhh)),'k-'); plot(locs2,100*(hhh(locs2)./sum(hhh)),'k*')
        title([List(ListNo).name(1:end-4)])
        set(gca,'FontSize',18,'FontName','Arial')
        axis tight
        subplot(3,1,3); hold on        
        histogram(hh,'FaceColor','w')
        yy = get(gca,'ylim');
        plot([h h],yy,'r-','LineWidth',3)     
        set(gca,'FontSize',18,'FontName','Arial')                  
        subplot(3,1,2); hold on        
        histogram(r2,'FaceColor','w')
        title(['Mean Raw R = ' num2str(round(mean(r2),2))])
        set(gca,'FontSize',18,'FontName','Arial')        
        set(gcf,'renderer','Painters')
        if (plot_paper_figures && ListNo==35)            
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_PeaksAndCorr' addlab])
        end            
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' List(ListNo).name(1:end-4) '_PeaksAndCorr' addlab])
        
        figure; hold on; 
        set(gcf,'Position',[  2098        -210         438         726])
        subplot(3,1,1); hold on
        yyaxis left
        plot(100*(hhh2(:,1)./sum(hhh2(:,1))),'r-'); plot(locs1{1},100*(hhh2(locs1{1},1)./sum(hhh2(:,1))),'r*')
        plot(100*(hhh2(:,2)./sum(hhh2(:,2))),'b-'); plot(locs1{2},100*(hhh2(locs1{2},2)./sum(hhh2(:,2))),'b*')
        plot(100*(hhh./sum(hhh)),'k-'); plot(locs2,100*(hhh(locs2)./sum(hhh)),'k*')
        yyaxis right    
        plot(1:size(InFR,2),mean(InFR+OutFR),'k--')         
        title([List(ListNo).name(1:end-4)])
        set(gca,'FontSize',18,'FontName','Arial')
        axis tight
        subplot(3,1,3); hold on        
        histogram(hh,'FaceColor','w')
        yy = get(gca,'ylim');
        plot([h h],yy,'r-','LineWidth',3)     
        set(gca,'FontSize',18,'FontName','Arial')                  
        subplot(3,1,2); hold on        
        histogram(r2,'FaceColor','w')
        title(['Mean Raw R = ' num2str(round(mean(r2),2))])
        set(gca,'FontSize',18,'FontName','Arial')        
        set(gcf,'renderer','Painters')
        if (plot_paper_figures && ListNo==35)            
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_PeaksAndCorr' addlab])
        end            
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' List(ListNo).name(1:end-4) '_PeaksAndCorr_PlusFields_' addlab])
        
    end
    
    if plot_other_figures
        fontS = 18;        
        num = 7;
        
        f = figure; hold on
        subplot(num,1,1); hold on
        histogram(hovpos(:,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['Posterior (all)'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)   

        subplot(num,1,2); hold on
        histogram(hovpos(hovpos(:,2)<6,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['First 5 laps'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)

        subplot(num,1,3); hold on
        histogram(hovpos(hovpos(:,2)>9 & hovpos(:,2)<15,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['Last 5 laps'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)

        subplot(num,1,4); hold on
        bar(mean(InFR+OutFR),'FaceColor','r')
        ylabel(['Fields (mean FR)'])        
        xlabel('Decoded Position Bin')
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)

        subplot(num,1,5); hold on
        histogram(hovpos2(:,1,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['Posterior (Sh)'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)   

        subplot(num,1,6); hold on
        histogram(hovpos2(hovpos(:,2,1)<6,1,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['Sh: First 5 laps'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)

        subplot(num,1,7); hold on
        histogram(hovpos(hovpos(:,2,1)>9 & hovpos(:,2,1)<15,1,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['Sh: Last 5 laps'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',fontS)

        suptitle([List(ListNo).name(1:end-4)])
        set(gcf,'Position',[ 2093        -732         690        1653])
        helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\' List(ListNo).name(1:end-4) '_Hist' addlab])
    end
    
    %%%%% add data to group data (real, shuffled, and running)
   r_all = cat(2,r_all,r2);
   r_running_all = cat(2,r_running_all,r3);
   r_running_all_sm = cat(2,r_running_all_sm,r32);
   realr = cat(2,realr,h);
   shuffr = cat(2,shuffr,hh);
end


if plot_paper_figures
    
    totp = (sum(mean(shuffr,1)>=mean(realr))+1)./(size(hovpos2,3)+1);

    figure; hold on
    r_all_mean = mean(r_all);
    [p,~,stats] = signrank(r_all_mean,r_running_all);
    h1 = histogram(r_all_mean,-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor','k','Normalization','probability');   
    h2 = histogram(r_running_all,-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor',[.5 .5 .5],'Normalization','probability');         
    yl = get(gca,'ylim');
    plot([median(r_running_all) median(r_running_all)],yl,'--','Color',[.5 .5 .5],'LineWidth',2)
    text(median(r_running_all)+.01,yl(2)*.98,num2str(round(median(r_running_all),2)),'Color',[.5 .5 .5],'FontSize',12)
    plot([median(r_all_mean) median(r_all_mean)],yl,'--','Color','k','LineWidth',2)
    text(median(r_all_mean)+.01,yl(2)*.98,num2str(round(median(r_all_mean),2)),'Color','k','FontSize',12)
    legend([h2 h1],'Running distribution','Flat field distribution','Location','northwest')
    xlabel('Correlation')
    ylabel('Proportion of sessions')
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_savefig([basedir 'Figures\HoversMovingShuffle\Hist_Running_And_Flat'])    
    helper_saveandclosefig([basedir 'Figures\Final\Hist_Running_And_Flat'])
    
    figure; hold on;
    text(mean(r_running_all_sm)+.01,.8,num2str(round(mean(r_running_all_sm),2)),'Color',[.5 .5 .5],'FontSize',12)   
    text(mean(r_running_all_sm)+.01,.5,['SEM = ' num2str(std(r_running_all_sm)./sqrt(sum(~isnan(r_running_all_sm))))],'Color',[.5 .5 .5],'FontSize',12)  
    text(mean(r_all_mean)+.01,.8,num2str(round(mean(r_all_mean),2)),'Color','k','FontSize',12)    
    text(mean(r_all_mean)+.01,.5,['SEM = ' num2str(std(r_all_mean)./sqrt(sum(~isnan(r_all_mean))))],'Color','k','FontSize',12)  
    title({['Means; N = ' num2str(length(r_running_all_sm)), ', signedrank = ' num2str(stats.signedrank)];['zval = ' num2str(stats.zval) ', p = ' num2str(p)]})    
    helper_savefig([basedir 'Figures\HoversMovingShuffle\Hist_Running_And_Flat_Stats'])
    helper_saveandclosefig([basedir 'Figures\Final\Hist_Running_And_Flat_Stats'])
    
    figure; hold on
    r_all_mean = mean(r_all);
    [p,~,stats] = signrank(r_all_mean,r_running_all_sm);
    h1 = histogram(r_all_mean,-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor','k','Normalization','probability');   
    h2 = histogram(r_running_all_sm,-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor',[.5 .5 .5],'Normalization','probability');         
    yl = get(gca,'ylim');
    plot([median(r_running_all_sm) median(r_running_all_sm)],yl,'--','Color',[.5 .5 .5],'LineWidth',2)
    text(median(r_running_all_sm)+.01,yl(2)*.98,num2str(round(median(r_running_all_sm),2)),'Color',[.5 .5 .5],'FontSize',12)
    plot([median(r_all_mean) median(r_all_mean)],yl,'--','Color','k','LineWidth',2)
    text(median(r_all_mean)+.01,yl(2)*.98,num2str(round(median(r_all_mean),2)),'Color','k','FontSize',12)
    legend([h2 h1],'Running distribution','Flat field distribution','Location','northwest')
    xlabel('Correlation')
    ylabel('Proportion of sessions')
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_savefig([basedir 'Figures\HoversMovingShuffle\Hist_Running_And_Flat_Sm'])    
    helper_saveandclosefig([basedir 'Figures\Final\Hist_Running_And_Flat_Sm'])
    
    figure; hold on;
    text(mean(r_running_all_sm)+.01,.8,num2str(round(mean(r_running_all_sm),2)),'Color',[.5 .5 .5],'FontSize',12)   
    text(mean(r_running_all_sm)+.01,.5,['SEM = ' num2str(std(r_running_all_sm)./sqrt(sum(~isnan(r_running_all_sm))))],'Color',[.5 .5 .5],'FontSize',12)  
    text(mean(r_all_mean)+.01,.8,num2str(round(mean(r_all_mean),2)),'Color','k','FontSize',12)    
    text(mean(r_all_mean)+.01,.5,['SEM = ' num2str(std(r_all_mean)./sqrt(sum(~isnan(r_all_mean))))],'Color','k','FontSize',12)  
    title({['Means; N = ' num2str(length(r_running_all_sm)), ', signedrank = ' num2str(stats.signedrank)];['zval = ' num2str(stats.zval) ', p = ' num2str(p)]})    
    helper_savefig([basedir 'Figures\HoversMovingShuffle\Hist_Running_And_Flat_Stats_Sm'])
    helper_saveandclosefig([basedir 'Figures\Final\Hist_Running_And_Flat_Stats_Sm'])

    figure; hold on; 
    histogram(mean(shuffr),'FaceColor','w','LineWidth',2); 
    yl = get(gca,'ylim'); 
    plot([nanmean(realr) nanmean(realr)],yl,'k','LineWidth',3)
    title(['p = ' num2str(round(totp,2,'significant'))])
    ylabel('Shuffle count')
    xlabel('Number of Peaks')
    set(gca,'FontSize',18)
    set(gcf,'Position',[  2044          58         874         595])
    set(gcf,'renderer','Painters')
    helper_savefig([dirs.figdir '\HoversMovingShuffle\Shuff_HistPeaks' addlab '_numshuff' num2str(size(hovpos2,3))])
    helper_saveandclosefig([dirs.figdir '\Final\MovingShuff_HistPeaks' addlab '_numshuff' num2str(size(hovpos2,3))])
end

if plot_other_figures
    figure; hold on; 
    plot(realr-nanmean(shuffr,1),'ok','MarkerSize',10)
    plot(1:length(realr),zeros(size(realr)),'r--','LineWidth',3)
    errorbar(length(realr)/2,mean(realr-nanmean(shuffr,1)),std(realr-nanmean(shuffr,1))./sqrt(length(realr)),'r','LineWidth',3)
    p = signrank(realr,nanmean(shuffr,1));
    title(['Signrank p = ' num2str(p)])
    xlabel('Session')
    ylabel('Number of peaks minus average shuffle peaks')
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    set(gcf,'Position',[680   386   760   592])
    helper_savefig([dirs.figdir '\HoversMovingShuffle\Shuff_PeakSessions' addlab '_numshuff' num2str(size(hovpos2,3))])
    helper_saveandclosefig([dirs.figdir '\Final\MovingShuff_PeakSessions' addlab '_numshuff' num2str(size(hovpos2,3))])

    figure; hold on
    histogram(r_all(:),'FaceColor','w')
    title(['Mean Raw R = ' num2str(round(mean(r_all(:)),2))])
    set(gca,'FontSize',18,'FontName','Arial')
    yl = get(gca,'ylim');
    plot([mean(mean(r_all)) mean(mean(r_all))],yl,'r-','LineWidth',2)
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\corr_all' addlab])

    figure; hold on
    histogram(mean(r_all),5,'FaceColor','w')
    title(['Mean Raw R = ' num2str(round(mean(mean(r_all)),2))])
    set(gca,'FontSize',18,'FontName','Arial')
    yl = get(gca,'ylim');
    plot([mean(mean(r_all)) mean(mean(r_all))],yl,'r-','LineWidth',2)
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\corr_mean' addlab])

    figure; hold on;
    plot(cellnum,mean(r_all),'ko','MarkerSize',15)
    [r,p] = corr(cellnum,mean(r_all)');
    yl = get(gca,'ylim'); xl = get(gca,'xlim');
    if p<.05
        text(mean(xl),mean(yl),['r = ' num2str(r) ', p = ' num2str(p)],'Color','r')
    else
        text(mean(xl),mean(yl),['r = ' num2str(r) ', p = ' num2str(p)],'Color','k')
    end
    xlabel('Number of Cells')
    ylabel('Posterior Correlation')
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\corr_vs_cells' addlab])
end


%%%%% using only sessions with more than 60 HP cells
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
    helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\corr_60cell_all' addlab])

    figure; hold on
    histogram(mean(r_all2),5,'FaceColor','w')
    title(['Over 60 Cells, Mean Raw R = ' num2str(round(mean(mean(r_all2)),2))])
    set(gca,'FontSize',18,'FontName','Arial')
    yl = get(gca,'ylim');
    plot([mean(mean(r_all2)) mean(mean(r_all2))],yl,'r-','LineWidth',2)
    set(gca,'FontSize',18)
    set(gcf,'renderer','Painters')
    helper_saveandclosefig([dirs.figdir '\HoversMovingShuffle\corr_60cell_mean' addlab])
end
