function PFC_tolongreplays(basedir,plot_paper_figures)
%%%%% looking at whether the pfc population activity change 
%%%%% stays elevated for longer when there are longer replays

%%%%% set up params and dirs
load([basedir 'dirs_ymaze_new.mat'],'dirs')
cd(dirs.spikedatadir)
if ~isfolder(dirs.figdir)
    mkdir(dirs.figdir)
end
pfcall = [];
wc_cutoff = .3; 
jd_cutoff = .7; 
coveragecutoff = 0;
d2 = dir('*.mat');
toadd = [];

numcells = NaN(size(d2,1),2);
for id = 1:size(d2,1)    
    load(d2(id).name,'other_cells','hp_cells','hpinterneurons')
    numcells(id,1) = sum(~ismember(hp_cells,hpinterneurons));
    numcells(id,2) = length(other_cells);
end

%%%%% get data out
for id = 1:size(d2,1)
    load(d2(id).name,'other_cells','spikedata','pos')
    spks = spikedata(ismember(spikedata(:,2),other_cells),:);
    EstBin = .04;
    bins = min(pos(:,1)):EstBin:max(pos(:,1));
    FR = NaN(length(bins),length(other_cells));
    for icell = 1:length(other_cells)
        FR(:,icell) = histc(spks(spks(:,2)==other_cells(icell),1),bins);
    end
    tm = histc(pos(:,1),bins);
    FR(tm==0,:) = [];
    bins(tm==0) = [];
    FRz = zscore(FR);

    load(d2(id).name,'CandSeq','CandCorr','CandDis','OutFR','-mat')
    
    CandPassCrit=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff; 

    window = [-1 1];
    windowbin = window./EstBin;
    CandSeq = CandSeq(CandPassCrit,:);
    
    pfc = NaN(size(CandSeq,1),range(windowbin)+1,length(other_cells));
    excl = false(size(CandSeq,1),1);
    for iseq = 1:size(CandSeq,1)-3
        f = find(bins-CandSeq(iseq,1)>0,1,'first')-1;
        if (f+windowbin(1))<=0 || (f+windowbin(2))>size(FRz,1)
            excl(iseq) = 1;
        else
            pfc(iseq,:,:) = FRz(f+windowbin(1):f+windowbin(2),:);
        end
    end
    CandSeq(excl,:) = []; 
    pfc(excl,:,:) = []; 
      
    pfc11 = pfc(diff(CandSeq,[],2)<median(diff(CandSeq,[],2)),:,:);
    pfc22 = pfc(diff(CandSeq,[],2)>median(diff(CandSeq,[],2)),:,:);
    pfc1 = cat(3,squeeze(nanmean(pfc11,1)),squeeze(nanmean(pfc22,1)));
    
    pfcall = cat(2,pfcall,pfc1);
end
%%

%%%%% plot group data
if plot_paper_figures
    figure; hold on

    a = plot(window(1):EstBin:window(2),nanmean(pfcall(:,:,1),2),'b.-','MarkerSize',15);
    mindat = nanmean(pfcall(:,:,1),2)-nanstd(pfcall(:,:,1),[],2)./sqrt(size(pfcall,2));
    maxdat = nanmean(pfcall(:,:,1),2)+nanstd(pfcall(:,:,1),[],2)./sqrt(size(pfcall,2));
    patch([window(1):EstBin:window(2) window(2):-EstBin:window(1)],[mindat' maxdat(end:-1:1)'],'blue','FaceAlpha',.3,'EdgeAlpha',0)

    maxdat = nanmean(pfcall(:,:,2),2)+nanstd(pfcall(:,:,2),[],2)./sqrt(size(pfcall,2));
    mindat = nanmean(pfcall(:,:,2),2)-nanstd(pfcall(:,:,2),[],2)./sqrt(size(pfcall,2));
    b = plot(window(1):EstBin:window(2),nanmean(pfcall(:,:,2),2),'r.-','MarkerSize',15);
    patch([window(1):EstBin:window(2) window(2):-EstBin:window(1)],[mindat' maxdat(end:-1:1)'],'red','FaceAlpha',.3,'EdgeAlpha',0)
    xlabel('Time since start of replay (sec)')
    ylabel('Z-scored FR of PFC cells')
    legend([a b],{'Short replays';'Long replays'})
    axis tight
    set(gca,'FontSize',26)
    set(gcf,'Position',[ 1960         207        1000         598])
    set(gcf, 'Renderer', 'painters');
    helper_savefig([basedir '\Figures\Final\PFC_tolongreplays_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_armcov' num2str(coveragecutoff) toadd])
    helper_saveandclosefig([dirs.figdir '\PFC_tolongreplays_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_armcov' num2str(coveragecutoff) toadd])


    ps = NaN(6,1); Zs = ps;
    for ii = 27:31
        [ps(ii-26,1),~,stats] = signrank(squeeze(pfcall(ii,:,1)),squeeze(pfcall(ii,:,2)),'tail','left');
        Zs(ii-26,1) = stats.zval;
    end
    [ps(6,1),~,stats] = signrank(squeeze(mean(pfcall(29:30,:,1))),squeeze(mean(pfcall(29:30,:,2))),'tail','left');
    Zs(6,1) = stats.zval;

    d1 = pfcall(27:28,:,1);
    d2 = pfcall(27:28,:,2);
    d3 = pfcall(29:30,:,1);
    d4 = pfcall(29:30,:,2);
    dat = [d1(:); d2(:); d3(:); d4(:)];
    shortlong = [ones(1,size(d1(:),1)) 2*ones(1,size(d1(:),1)) ones(1,size(d1(:),1)) 2*ones(1,size(d1(:),1))];
    earlylate = [ones(1,size(d1(:),1)*2) 2*ones(1,size(d1(:),1)*2)];
    [p,tbl] = anovan(dat,[shortlong' earlylate'],'model','interaction','display','off');

    figure; hold on
    text(.1,.9,['80-120ms, Z = ' num2str(Zs(3)) ', P = ' num2str(ps(3)) ', N = 158,158'])
    text(.1,.8,['120-160ms, Z = ' num2str(Zs(4)) ', P = ' num2str(ps(4)) ', N = 158,158'])
    text(.1,.7,['80-160ms, Z = ' num2str(Zs(6)) ', P = ' num2str(ps(6)) ', N = 158,158'])
    text(.1,.6,['Two-way anova between first 80ms and next 80ms, F(1,' num2str(length(dat)) ')'])
    text(.3,.5,['Main effect shortlong, Z = ' num2str(tbl{2,6}) ' p = ' num2str(p(1))])
    text(.3,.4,['Main effect earlylate, Z = ' num2str(tbl{3,6}) ' p = ' num2str(p(2))])
    text(.4,.3,['Interaction, Z = ' num2str(tbl{4,6}) ' p = ' num2str(p(3))])
    helper_saveandclosefig([basedir 'Figures/Final/PFC_tolongreplays_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_armcov' num2str(coveragecutoff) '_stats' toadd])
    helper_saveandclosefig([dirs.figdir 'PFC_tolongreplays_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_armcov' num2str(coveragecutoff) '_stats' toadd])
end
