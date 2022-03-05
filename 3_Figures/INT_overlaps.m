function INT_overlaps(basedir,plot_paper_figures,plot_other_figures)

%%%%% Looking at how the modulation of interneurons chance across laps

%%
%%%%%
%%%%%
%%%%% set up parameters and directories
%%%%%
%%%%%

load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')
cd(dirs.spikedatadir)
if ~isfolder([dirs.figdir 'LookAtINT'])
    mkdir([dirs.figdir 'LookAtINT'])
end
numlaps = 26;
wc_cutoff = .3;
jd_cutoff = .7;
coveragecutoff = 0;
bsline = .3;
negmod = .1;
posmod = .1;
addlab = ['1110_negmod' num2str(negmod) '_posmod' num2str(posmod) '_bsline' num2str(bsline) '_Pearson_1.1abvmean'];
List=dir('*.mat');
excl_ind = true(size(List,1),1);

EstBin = .02;
allbs = []; allbs2 = [];
allr = [];
numcells = NaN(size(List,1),10); % median(numcells(:,9)) is 30
firstcutoff = 8;

%%

%%%%%
%%%%%
%%%%% get data out and plot example neurons
%%%%%
%%%%%

for ListNo= 1:size(List,1)
    if excl_ind(ListNo)==0
        continue
    end
    load(List(ListNo).name,'hp_cells','hpinterneurons')
    if exist('hp_cells','var')
        numcells(ListNo,1) = length(hpinterneurons);        
    end
    load(List(ListNo).name,'CandCorr','CandDis','OutFR','CandStepSize','params','MidTime','-mat')
    CandPassCrit=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; %new
    numcells(ListNo,2)=sum(CandPassCrit); 
    numcells(ListNo,3) = params.Novel;
    numcells(ListNo,9) = length(MidTime);
    numcells(ListNo,10) = sum(~ismember(hp_cells,hpinterneurons));
    
    if numcells(ListNo,1)>0 && params.Novel==1 
       load(List(ListNo).name,'CandSeq','Spike','pos')
        CandSeqsave = CandSeq;
        
        CandSeq = [CandSeq(:,1)-.4 CandSeq(:,1)+.4];
        
        Cand=CandSeq;
        N=round(diff(CandSeq')/EstBin)';
        t=find(mod(N,2));
        Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*floor(N(t)/2)-EstBin/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*floor(N(t)/2)+EstBin/2];
        t=find(~mod(N,2));
        Cand(t,:)=[CandSeq(t,1)/2+CandSeq(t,2)/2-EstBin*N(t)/2 CandSeq(t,1)/2+CandSeq(t,2)/2+EstBin*N(t)/2];
        clear N t

        [~,I]=histc(Spike(:,2),sortrows(Cand(:)));

        B=cumsum([0 diff(Cand')]);
        Index=round(B/(EstBin));


        Orig_Spike=Spike;
    
        for i=1:2:max(I)
            Spike(I==i,2)=Spike(I==i,2)+B(ceil(i/2))-Cand(ceil(i/2),1); 
        end
        Spike(~mod(I,2),:)=[];
        clear Cand I i

        Start_Time=0;
        End_Time=B(end);
        TimeBins=round((End_Time-Start_Time)/EstBin);

        % Bayesian Decoding - 5ms moving step, 20ms window estimate
        binspike=zeros(length(hpinterneurons),TimeBins);
        for CellID=1:length(hpinterneurons)                
            c=histc(Spike(Spike(:,1)==hpinterneurons(CellID),2),Start_Time:EstBin:End_Time);
            binspike(CellID,:)=c(1:TimeBins);            
        end
        binspikes = reshape(binspike,[length(hpinterneurons) size(binspike,2)./size(CandSeq,1) size(CandSeq,1)]);
        [~,I]=histc(CandSeqsave(:,1),[MidTime(1:end);max(pos(:,1))]); 
        firstlaps = I>0 & I<=firstcutoff;
        lastlaps = I>firstcutoff & I<=26;
        bsind = size(binspikes,2)/2-(negmod/EstBin):size(binspikes,2)/2+(posmod/EstBin);
        bsind1 = 1:(bsline/EstBin);
        bs = reshape(mean(binspikes(:,bsind,:),2),[size(binspikes,1) size(binspikes,3)]);
        m = reshape(mean(binspikes(:,size(binspikes,2)/2:size(binspikes,2)/2+20,:),2),[size(binspikes,1) size(binspikes,3)]);
        bs2 = reshape(mean(binspikes(:,bsind1,:),2),[size(binspikes,1) size(binspikes,3)]);
        include = true(size(mean(m,2))); % all cells
        bs = bs(include,:); bs2 = bs2(include,:);
        hpinterneurons = hpinterneurons(include); binspikes = binspikes(include,:,:);
        if sum(I<=numlaps)==0 || sum(include)==0
            continue
        end
        jnk = repmat(I',[size(bs,1) 1]);
        allbs = cat(1,allbs,cat(2,bs(:),jnk(:))); %modtime
        allbs2 = cat(1,allbs2,cat(2,bs2(:),jnk(:))); %baseline
        ind = jnk(1,:)<=numlaps & jnk(1,:)>0;
        [r,p] = corr([(bs(:,ind)-bs2(:,ind))./(bs(:,ind)+bs2(:,ind))]',jnk(:,ind)','rows','complete');
        allr = cat(1,allr,[r(:,1) p(:,1)]);
        
        if plot_other_figures || (plot_paper_figures && ListNo==31) 
            for icell = 1:size(binspikes,1)            
                if plot_other_figures || (plot_paper_figures && ListNo==31 && icell==1) 
                    figure; hold on;                        
                    dat1 = squeeze(binspikes(icell,:,lastlaps));
                    fwd = [mean(dat1,2)-std(dat1,[],2)./sqrt(size(dat1,2))]';
                    rev = [mean(dat1,2)+std(dat1,[],2)./sqrt(size(dat1,2))]';
                    p1 = patch([1:size(binspikes,2) size(binspikes,2):-1:1],[fwd rev(end:-1:1)],'black');
                    p1.FaceAlpha=.1; p1.EdgeAlpha=0;            
                    plot(mean(dat1,2),'k')
                    dat2 = squeeze(binspikes(icell,:,firstlaps));
                    fwd = [mean(dat2,2)-std(dat2,[],2)./sqrt(size(dat2,2))]';
                    rev = [mean(dat2,2)+std(dat2,[],2)./sqrt(size(dat2,2))]';
                    p2 = patch([1:size(binspikes,2) size(binspikes,2):-1:1],[fwd rev(end:-1:1)],'red');
                    p2.FaceAlpha=.1; p2.EdgeAlpha=0;            
                    plot(mean(dat2,2),'r')                                                
                    yl = get(gca,'ylim');

                    plot([size(binspikes,2)/2 size(binspikes,2)/2],[0 yl(2)],'b','LineWidth',2)
                    b1 = plot(mean(squeeze(binspikes(icell,:,firstlaps)),2),'r');            
                    b2 = plot(mean(squeeze(binspikes(icell,:,lastlaps)),2),'k');
                    plot(bsind,zeros(size(bsind)),'b','LineWidth',3)
                    plot(bsind1,zeros(size(bsind1)),'k','LineWidth',3)
                    set(gca,'xtick',[1 size(binspikes,2)/2 size(binspikes,2)],'xticklabel',[-.4 0 .4])
                    xlabel('Seconds from replay start')
                    ylabel('Interneuron FR')
                    legend([b1 b2],{['First ' num2str(firstcutoff) ' laps'];['Laps ' num2str(firstcutoff+1) '-' num2str(numlaps) ]},'Location','northwest')

                    set(gcf,'Position',[2042        -351         802         490])
                    yl = get(gca,'ylim');
                    ylim([0 yl(2)])
                    xlim([1 size(binspikes,2)])
                    set(gca,'FontSize',26)
                    suptitle([List(ListNo).name(1:end-4) ', Cell ' num2str(hpinterneurons(icell)) ' r = ' num2str(round(r(icell,1),2,'significant')) ' p = ' num2str(round(p(icell,1),2,'significant'))])
                    set(gcf,'renderer','Painters')
                    if (plot_paper_figures && ListNo==31 && icell==1) 
                        helper_savefig([dirs.figdir 'Final\OverLaps_' List(ListNo).name(1:end-4) '_Cell' num2str(hpinterneurons(icell)) '_firstcutoff' num2str(firstcutoff)])
                    end
                    helper_saveandclosefig([dirs.figdir 'LookAtINT\OverLaps_' List(ListNo).name(1:end-4) '_Cell' num2str(hpinterneurons(icell)) '_firstcutoff' num2str(firstcutoff)])
                end
            end
        end
    end
end

%%

%%%%%
%%%%%
%%%%% group data histogram of correlation between modulation and lap
%%%%%
%%%%%

if plot_paper_figures
    [r_mod,p_mod] = corr(allbs(allbs(:,2)<=numlaps,1),allbs(allbs(:,2)<=numlaps,2),'rows','complete');
    [r_bsl,p_bsl] = corr(allbs2(allbs(:,2)<=numlaps,1),allbs(allbs(:,2)<=numlaps,2),'rows','complete');

    ['numlaps ' num2str(numlaps) ', mod ' num2str(r_mod) ' ' num2str(p_mod) ', bsl ' num2str(r_bsl) ' ' num2str(p_bsl)];

    [p,~,stats] = signrank(allr(:,1),0,'tail','right');
    figure; hold on; 
    histogram(allr(:,1),-.4:.1:.4,'FaceColor','w', 'LineWidth',3)
    histogram(allr(allr(:,2)<.05,1),-.4:.1:.4, 'LineWidth',3)
    set(gca,'FontSize',26)
    xlabel('Correlation between modulation and lap')
    ylabel('Interneurons')    
    title({['N = ' num2str(size(allr,1)) ', signedrank = ' num2str(stats.signedrank)];['Z = ' num2str(stats.zval) ', P = ' num2str(p)]})
    set(gcf,'Position',[2044          91         751         488])
    set(gcf,'renderer','Painters')
    helper_savefig([dirs.figdir 'LookAtINT\OverLaps_all_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_' num2str(numlaps) 'laps_' addlab])
    helper_saveandclosefig([dirs.figdir 'Final\OverLaps_all_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_' num2str(numlaps) 'laps_' addlab])

    disp(['p = ' num2str(p)])
end