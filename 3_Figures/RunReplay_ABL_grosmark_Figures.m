function RunReplay_ABL_grosmark_Figures(basedir,plot_paper_figures,plot_other_figures)

%%%%% plot interneurons' modulation across laps in grosmark data


%%
%%%%%
%%%%% dirs and params
%%%%%
load([basedir 'dirs_linear_grosmark'],'dirs')
cd(dirs.spikedatadir)
if ~isfolder(dirs.figdir)
    mkdir(dirs.figdir)
end
if ~isfolder([dirs.figdir 'INT/'])
    mkdir([dirs.figdir 'INT/'])
end


%%

%%%%%
%%%%% replicate INT effect (from INT_overlaps.m)
%%%%%

load([basedir 'dirs_linear_grosmark'],'dirs')
cd(dirs.spikedatadir)
List=dir('*.mat');
jd_cutoff = .7; wc_cutoff = .3; coveragecutoff = 0;
EstBin = .02;
allbs = []; allbs2 = [];
allr = [];
numlaps = 26;
numcells = NaN(size(List,1),10);
% min(numcells(:,9)) is 57
% median(numcells(:,9)) is 83
bsline = .3;
negmod = .1;
posmod = .1;
firstcutoff = 8;
addlab = ['_negmod' num2str(negmod) '_posmod' num2str(posmod)...
    '_bsline' num2str(bsline) '_Pearson_1.1abvmean_firstcutoff_'...
    num2str(firstcutoff)];

for ListNo= 1:size(List,1)
    
    load(List(ListNo).name,'hp_cells','hpinterneurons')
    if exist('hp_cells','var')
        numcells(ListNo,1) = length(hpinterneurons);        
    end
    load(List(ListNo).name,'CandCorr','CandDis','OutFR','CandStepSize','params','MidTime','-mat')
    CandPassCrit=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; %new
    numcells(ListNo,2)=sum(CandPassCrit);
    numcells(ListNo,3) = params.Novel;
    numcells(ListNo,9) = length(MidTime)-1;
    numcells(ListNo,10) = sum(~ismember(hp_cells,hpinterneurons));
    
    if numcells(ListNo,1)>0 && params.Novel==1 && sum(CandPassCrit)>1
       load(List(ListNo).name,'CandSeq','Spike','pos')
       CandSeq = CandSeq(CandPassCrit);
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
        lastlaps = I>firstcutoff & I<=numlaps;
        bsind = size(binspikes,2)/2-(negmod/EstBin):size(binspikes,2)/2+(posmod/EstBin);
        bsind1 = 1:(bsline/EstBin);
        bs = reshape(mean(binspikes(:,bsind,:),2),[size(binspikes,1) size(binspikes,3)]);
        m = reshape(mean(binspikes(:,size(binspikes,2)/2:size(binspikes,2)/2+20,:),2),[size(binspikes,1) size(binspikes,3)]);
        bs2 = reshape(mean(binspikes(:,bsind1,:),2),[size(binspikes,1) size(binspikes,3)]);
        include = true(size(m,1),1);
        bs = bs(include,:); bs2 = bs2(include,:);
        hpinterneurons = hpinterneurons(include); binspikes = binspikes(include,:,:);
        if sum(include)==0 || sum(I<=numlaps)==0
            continue
        end
        jnk = repmat(I',[size(bs,1) 1]);
        allbs = cat(1,allbs,cat(2,bs(:),jnk(:)));
        allbs2 = cat(1,allbs2,cat(2,bs2(:),jnk(:)));
        ind = jnk(1,:)<=numlaps & jnk(1,:)>0;
        
        celljnk = [(bs(:,ind)-bs2(:,ind))./(bs(:,ind)+bs2(:,ind))]'; %baseline
        r = NaN(length(hpinterneurons),1); p = r;
        for icell = 1:length(hpinterneurons)
            [r(icell),p(icell)] = corr(celljnk(:,icell),jnk(icell,ind)','rows','complete');
        end
        allr = cat(1,allr,[r(:,1) p(:,1)]);
        
        if plot_other_figures || (plot_paper_figures && ListNo==1) 
            for icell = 1:size(binspikes,1)   
                if plot_other_figures || (plot_paper_figures && icell==8 && ListNo==1) 
                     figure; hold on;                        
                    dat1 = squeeze(binspikes(icell,:,lastlaps));
                    if size(dat1,1)==1; dat1 = dat1'; end
                    fwd = [mean(dat1,2)-std(dat1,[],2)./sqrt(size(dat1,2))]';
                    rev = [mean(dat1,2)+std(dat1,[],2)./sqrt(size(dat1,2))]';
                    p1 = patch([1:size(binspikes,2) size(binspikes,2):-1:1],[fwd rev(end:-1:1)],'black');
                    p1.FaceAlpha=.1; p1.EdgeAlpha=0;            
                    plot(mean(dat1,2),'k')
                    dat2 = squeeze(binspikes(icell,:,firstlaps));
                    if size(dat2,1)==1; dat2 = dat2'; end
                    if ~(size(dat2,2)<2 || size(dat2,1)<2)

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

                        set(gcf,'Position',[ 2042        -351         802         490])
                        yl = get(gca,'ylim');
                        ylim([0 yl(2)])
                        xlim([1 size(binspikes,2)])
                        set(gca,'FontSize',26)
                        suptitle([List(ListNo).name(1:end-4) ', Cell ' num2str(hpinterneurons(icell)) ' r = ' num2str(round(r(icell,1),2,'significant')) ' p = ' num2str(round(p(icell,1),2,'significant'))])
                        set(gcf,'renderer','Painters')
                        if (plot_paper_figures && icell==8 && ListNo==1) 
                            helper_savefig([basedir '\Figures\Final\OverLaps_' List(ListNo).name(1:end-4) '_Cell' num2str(hpinterneurons(icell)) '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) addlab])
                        end                    
                        helper_saveandclosefig([dirs.figdir '\INT\OverLaps_' List(ListNo).name(1:end-4) '_Cell' num2str(hpinterneurons(icell)) '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) addlab])
                    else
                        close all
                    end
                end
            end
        end
    end
end

%%

%%%%%
%%%%% plot
%%%%%

[p,~,stats] = signrank(allr(:,1),0,'tail','right');

if plot_paper_figures
    figure; hold on; 
    mm = max([abs(min(allr(:,1))) abs(max(allr(:,1)))]);
    histogram(allr(:,1),-mm:.1:mm,'FaceColor','w', 'LineWidth',3)
    histogram(allr(allr(:,2)<.05,1),-mm:.1:mm, 'LineWidth',3)
    set(gca,'FontSize',26)
    xlabel('Correlation between modulation and lap')
    ylabel('Interneurons')
    yl = get(gca,'ylim');
    t = text(-.35,round(yl(2)*.75),['p = ' num2str(round(p,2,'significant'))],'FontSize',20);
    % if p<.05
    %     t.Color = 'r';
    % end
    title({['N = ' num2str(size(allr,1)) ', signedrank = ' num2str(stats.signedrank)];['Z = ' num2str(stats.zval) ', P = ' num2str(p)]})
    set(gcf,'Position',[2044          91         751         488])
    xlim([-mm mm])
    set(gcf,'renderer','Painters')
    helper_savefig([dirs.figdir '\OverLaps_all' addlab '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_' num2str(numlaps) 'laps_bsline' num2str(bsline) '_negmod' num2str(negmod) '_posmod' num2str(posmod)])
    helper_saveandclosefig([basedir '\Figures\Final\OverLaps_all_grosmark' addlab '_jd' num2str(jd_cutoff) '_wc' num2str(wc_cutoff) '_cov' num2str(coveragecutoff) '_' num2str(numlaps) 'laps_bsline' num2str(bsline) '_negmod' num2str(negmod) '_posmod' num2str(posmod)])
end