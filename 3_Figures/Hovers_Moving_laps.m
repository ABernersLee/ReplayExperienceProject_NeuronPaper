function Hovers_Moving_laps(basedir,plot_paper_figures,plot_other_figures)

%%%%% Uses shuffles to determine whether hover distributions are more 
%%%%% stable than you would expect by chance across experience/slope
%%%%% changes


%%%%%
%%%%%
%%%%% set up parameters and dirs
%%%%%
%%%%%

close all
wc_cutoff = 0; %.6;
jd_cutoff = 1; %.4;
coveragecutoff = 0;
controlhovernumber = true;

load([basedir 'dirs_linear_lapbylap_addedIN_wholesession.mat'],'dirs')
cd(dirs.spikedatadir)
if ~controlhovernumber
    addlab = ['_wholefields_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_coverage' num2str(coveragecutoff)];
elseif controlhovernumber
    addlab = ['_wholefields_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_coverage' num2str(coveragecutoff) '_ControlHoverNumber'];
end

if ~isfolder([dirs.figdir '/HoversMoving/'])
    mkdir([dirs.figdir '/HoversMoving/'])
end

pvalues = []; RR = []; SS = [];
List=dir('*.mat');
numses = size(List,1);
numperm = 1000;  
numshuff = 1000;

FirstLapP = NaN(numses,1);
FirstLapNum = NaN(numses,1);
FirstLapShuff = NaN(numses,numshuff);
Bin_Size = 2.5;
vellab = {'FirstLap';'FirstTwoLaps';'AllLaps'};
allvelreal = NaN(numses,3);
alloccreal = NaN(numses,3);
allvelp = NaN(numses,3);
allvelshuff = NaN(numses,3,numshuff);
firstlapilaps = NaN(numses,1);

%%%%%
%%%%%
%%%%% data and processing
%%%%%
%%%%%

%%%%% saves out only one example session if plot_other_figures is 0

for ListNo= 1:numses 
    disp(['Session ' num2str(ListNo)])
    load(List(ListNo).name,'InFR',...
        'hover*','move*','Index','InMatrix','OutMatrix',...
        'CandStepSize','Spike',...
        'CandCorr','CandDis','OutFR','CandSeq',...
        'MidTime','params',...
        'hover_sm*','move_sm*','hover_fwd*','move_fwd*',...
        'hovertime2','movetime2','vel','linposcat','pos','-mat') 

    if params.Novel~=1
        continue
    end
    
    tm = pos(:,1);
    
    stepdur = NaN(size(CandSeq,1),1);
    for ii = 1:size(CandSeq,1)
          stepdur(ii,1) = mean(diff([movetime2(movetime2(:,3)==ii,1)]))*5;
    end
    
    maxmove = CandDis/size(OutFR,2);
    CandPassCrit = abs(CandCorr)>wc_cutoff & maxmove<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; 

    [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);
    if sum(I==0)>0
        I(I==0) = NaN;
    end
    I(I>length(MidTime)) = length(MidTime);
        
    Mat=OutMatrix'+InMatrix';
    D=sum([1:size(Mat,1)]'*ones(1,size(Mat,2)).*Mat);
 
    t = find(CandPassCrit); 
    hovpos = []; posall = []; velall = [];
    for it = 1:length(t)
       pos = D(Index(t(it))*4+1:Index(t(it)+1)*4-3);
       posall = [posall;pos'];
        posdiff = diff(pos');        
       velall = [velall;[abs(posdiff) ones(size(abs(diff(pos'))))*I(t(it))]];
       hovertime3 = hovertime2(hovertime2(:,3)==t(it),:);
       ht = hovertime3(:,1:2);
       ind = zeros(size(pos));
       ind(ht(:,1)) = 1;
       ind(ht(:,2)+1) = -1;
       hovpos = [hovpos; pos(cumsum(ind)==1)' ones(sum(cumsum(ind)==1),1)*I(t(it)) ones(sum(cumsum(ind)==1),1)*t(it)];
    end
    if isempty(hovpos) || length(unique(hovpos(:,2)))<=1
        continue
    end
        

    Slopes = 1:length(I);
    Slopes(I>14) = NaN;
    num = 4;
    m = nanmedian(Slopes);
    um = nanmedian(Slopes(Slopes>m));
    lm = nanmedian(Slopes(Slopes<m));
    i = NaN(size(Slopes));
    i(Slopes<=lm) = 1; i(Slopes>lm & Slopes<=m) = 2;
    i(Slopes>m & Slopes<=um) = 3; i(Slopes>=um) = 4;

    
    hind = floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1)));
    hvals = NaN(length(hind)-1,num);
    

    PP = 1.5; PD = .1;
    tolap =0; ilap = 1;
    while tolap ==0 && ilap<15
        if sum(hovpos(:,2)==ilap)>0
            h = histcounts(hovpos(hovpos(:,2)==ilap,1),hind);
            [~,locs1] = findpeaks([0 h 0],'MinPeakProminence',PP,'MinPeakDistance',size(h,2)*PD);
            locs1 = locs1-1;
            tolap = 1;
        else
            ilap= ilap+1;
        end
    end
    
    if ilap<10 
        Lap1replays = unique(hovpos(hovpos(:,2)==ilap,3));
        numhovers_lap1 = length(hovpos(hovpos(:,2)==ilap,3));
        otherreplays = unique(hovpos(hovpos(:,2)>ilap,3));
        
        % changed to only take replays that cumulatively have as many hovers
        % as on lap ilap - Bill's suggestion from lab meeting 4/10/2020

        if length(otherreplays)>length(Lap1replays)
            if nchoosek(length(otherreplays),length(Lap1replays))>=numshuff
                locs_shuff = NaN(numshuff,1);
                for iS = 1:numshuff
                    shuffrep = otherreplays(randperm(length(otherreplays)));
                    newRep = shuffrep(1:length(Lap1replays));
                    takeind = find(ismember(hovpos(:,3),newRep));
                    if controlhovernumber
                        track = 0;
                        while length(takeind)<numhovers_lap1
                            track = track+1;
                            shuffrep = otherreplays(randperm(length(otherreplays)));
                            newRep = shuffrep(1:length(Lap1replays));
                            takeind = find(ismember(hovpos(:,3),newRep));
                            if mod(track,100)==0
                                disp('In while loop at 141 still')
                            end
                        end
                        clear track
                        if length(takeind)>numhovers_lap1
                            inds = randperm(length(takeind));
                            takeind = takeind(inds(1:numhovers_lap1));
                        end
                    end
                    h = histcounts(hovpos(takeind,1),hind);
                    [~,ll] = findpeaks([0 h 0],'MinPeakProminence',PP,'MinPeakDistance',size(h,2)*PD);
                    ll = ll-1;
                    locs_shuff(iS) = length(ll);
                end
                pFirstLap = (1+(sum(locs_shuff<=length(locs1))))./(1+numshuff);
                if plot_other_figures
                    figure; hold on;
                    histogram(locs_shuff,'FaceColor','w','LineWidth',2);
                    yl = get(gca,'ylim');
                    plot([length(locs1); length(locs1)],[yl]','r-','LineWidth',3)                       
                    title(['p = ' num2str(pFirstLap)])
                    set(gcf,'Position',[2143         152         560         420])
                    set(gcf, 'Renderer', 'painters');
                    helper_saveandclosefig([dirs.figdir '\HoversMoving\' List(ListNo).name(1:end-4) '_Hist_FirstLap_' addlab])
                end

                FirstLapP(ListNo) = pFirstLap;
                FirstLapNum(ListNo) = length(locs1);
                FirstLapShuff(ListNo,:) =locs_shuff;

                p2 = 1-binocdf(nansum(FirstLapP<.05)-1,sum(~isnan(FirstLapP)),.05);    
                p3 = (1+sum(nanmean(FirstLapShuff)<=nanmean(FirstLapNum)))./(1+size(FirstLapShuff,2));
                disp(['First Lap: p1 = ' num2str(p2) ', p2 = ' num2str(p3)])
                firstlapilaps(ListNo) = ilap;
            else
                firstlapilaps(ListNo) = 1000;
            end
        else
            firstlapilaps(ListNo) = 100;
        end
    else
        firstlapilaps(ListNo) = 10;
    end
    firstlap = ilap;    
    numlaps = 10;
    
    %is velocity related to where hovers are?
    lapind = getfirstlap(tm,linposcat,MidTime);
    clear pos
    posbinned1 = ceil(((linposcat-min(linposcat))+.01)./Bin_Size);    
    [~,posbinned] = histc(posbinned1, hind);
    posbinned(posbinned1<min(hind) & posbinned1>max(hind)) = NaN;

    posvel2 = NaN(max(posbinned),3);
    
    for ipos = 1:max(posbinned)
        posvel2(ipos,1) = mean(vel(posbinned==ipos & lapind(:,1)));
        posvel2(ipos,2) = mean(vel(posbinned==ipos & (lapind(:,1) | lapind(:,2))));
        posvel2(ipos,3) = mean(vel(posbinned==ipos));
    end
    h1 = histc(hovpos(:,1),hind);
    hOcc1 = histc(posbinned,hind);
    
    range1 = [ceil(length(hind)*.3):floor(length(hind)*.7)];
    posvel = posvel2(range1,:);
    h = h1(range1);
    h(h==0) = NaN;
    hOcc = hOcc1(range1);
    
    velreal = NaN(3,1); 
    occreal = NaN(3,1);
    pvel = NaN(3,1); 
    velshuff = NaN(3,numshuff);
    for ivel = 1:3
        velreal(ivel) = corr(posvel(:,ivel),h,'rows','complete');
        occreal(ivel) = corr(diff(posvel(:,ivel)),diff(h),'rows','complete');
        for ishuff = 1:numshuff
           velshuff(ivel,ishuff) = corr(posvel(:,ivel),h(randperm(length(h))),'rows','complete');
        end
        if velreal(ivel)>0
            pvel(ivel) = (1+sum(velshuff(ivel,:)>=velreal(ivel)))./(1+numshuff);
        elseif velreal(ivel)<0
            pvel(ivel) = -(1+sum(velshuff(ivel,:)<=velreal(ivel)))./(1+numshuff);
        end      
    end
    alloccreal(ListNo,:) = occreal;
    allvelreal(ListNo,:) = velreal;
    allvelp(ListNo,:) = pvel;
    allvelshuff(ListNo,:,:) = velshuff;
    
    
    if plot_other_figures || ListNo==35
        g = figure; hold on
        subplot(numlaps,1,numlaps); hold on    
        h = histogram(hovpos(:,1),hind,'FaceColor','w','LineWidth',2);
        ylabel(['Whole Session'])
        [~,locs,~,~] = findpeaks((h.Values./sum(h.Values))*100,'MinPeakProminence',.4,'MinPeakDistance',size(h.Values,2)*.10);
        locs2 = h.BinEdges(locs)+(h.BinWidth)./2;
        yl = get(gca,'ylim');
        plot([locs2; locs2],[ones(size(locs,2),1)*yl]','r--')

        for ilap = 1:numlaps
            subplot(numlaps+1,1,ilap), hold on
            h = histogram(hovpos(hovpos(:,2)==ilap,1),hind,'FaceColor','w','LineWidth',2);
            yl = get(gca,'ylim');
            plot([locs2; locs2],[ones(size(locs,2),1)*yl]','r--')
            if ilap==firstlap
                locs11 = h.BinEdges(locs1)+(h.BinWidth)./2;
                plot([locs11; locs11],[ones(size(locs1,2),1)*yl]','g--')
            end
            ylabel(['Lap ' num2str(ilap)])
        end    
        set(gcf,'Position',[ 2078        -803         455        1692])       
        a = replace(addlab,'_',' ');
        suptitle(a)
        set(gcf, 'Renderer', 'painters');
        if ListNo==35            
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_Hist_EachLap_' addlab])        
        end
        helper_saveandclosefig([dirs.figdir '\HoversMoving\' List(ListNo).name(1:end-4) '_Hist_EachLap_' addlab])
        close all
    end
    
    f = figure; hold on    
    for islope = 1:num
        if plot_other_figures || ListNo==35
            subplot(num+3,1,islope); hold on
        end
        h = histogram(hovpos(ismember(hovpos(:,3),find(i==islope)),1),hind,'FaceColor','w','LineWidth',2);
        hvals(:,islope) = h.Values;
        if plot_other_figures || ListNo==35
            ylabel('Count')
            xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
            set(gca,'FontSize',26)
        end
    end
    if ~(plot_other_figures || ListNo==35)
        close(f)
    end
    if any(sum(hvals)==0) || max(velall(:,2))<14
        close all
        continue
    end
    
    combs = nchoosek(1:4,2);           
    Rreal2 = corrcoef(hvals);
    Rreal = Rreal2([2:4 7:8 12])';
    
    indforever = repmat(1:size(hvals,1)',[1 1000]);
    pind = randi(1000-size(hvals,1),1,numperm,4);
    indf = NaN(size(hvals,1),4,numperm);
    for id = 1:4
        for it = 1:size(hvals,1)
            %ABL changed from this 5/30/20
%             indf(it,id,:) = indforever(squeeze(pind(1,:,id))+it-1);
            indf(it,id,:) = indforever(squeeze(pind(1,:,id))+it-1)+((id-1)*size(hvals,1)); 
        end
    end
    
    Rshift = NaN(size(combs,1),numperm);
    for id = 1:numperm
        %ABL changed 5/30/20
%         Rshift2 = corrcoef(indf(:,:,id));
        Rshift2 = corrcoef(hvals(indf(:,:,id))); 
        Rshift(:,id) = Rshift2([2:4 7:8 12]);
    end
    if any(isnan(Rreal))
        disp('*')
    end
    thisp = (sum(nanmean(Rshift)>=nanmean(Rreal))+1)./(numperm+1);
    pvalues = [pvalues;thisp];
    RR = cat(1,RR,mean(Rreal));
    SS = cat(1,SS,mean(Rshift));
    
    if plot_other_figures || ListNo==35
        subplot(num+3,1,islope+1); hold on
        histogram(hovpos(:,1),floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k','LineWidth',2)
        ylabel(['Hovers'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',26)

        subplot(num+3,1,islope+2); hold on
        histogram(posall,floor(min(hovpos(:,1))):1:ceil(max(hovpos(:,1))),'FaceColor','k')
        ylabel(['Posterior (all)'])        
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',26)

        subplot(num+3,1,islope+3); hold on
        bar(mean(InFR+OutFR),'FaceColor','r')
        ylabel(['Fields (mean FR)'])        
        xlabel('Decoded Position Bin')
        xlim([floor(min(hovpos(:,1))) ceil(max(hovpos(:,1)))])
        set(gca,'FontSize',26)

        suptitle([List(ListNo).name(1:end-4) ', p = ' num2str(thisp)])
        set(gcf,'Position',[ 1921        -799        1080        1803])
        set(gcf, 'Renderer', 'painters');
        if ListNo==35        
            helper_savefig([dirs.figdir '\Final\' List(ListNo).name(1:end-4) '_Hist_Laps_' addlab])
        end
        helper_saveandclosefig([dirs.figdir '\HoversMoving\' List(ListNo).name(1:end-4) '_Hist_Laps_' addlab])
    end
    
    
    disp(['Hovers Moving: p1 = ' num2str(1-binocdf(sum(pvalues<.05)-1,length(pvalues),.05)) ',p2 = ' num2str((sum(nanmean(SS)>=nanmean(RR))+1)./(sum(~isnan(nanmean(SS)))+1))])
    disp([num2str(sum(~isnan(FirstLapP))) ' Sessions used for FirstLap'])
    
end

%%%%%
%%%%%
%%%%% plot group data
%%%%%
%%%%%
if plot_paper_figures
    figure; hold on
    histogram(nanmean(FirstLapShuff),'FaceColor','w','LineWidth',2);
    yl = get(gca,'ylim');
    plot([nanmean(FirstLapNum) nanmean(FirstLapNum)],yl,'r-','LineWidth',3)
    p3 = (1+sum(nanmean(FirstLapShuff)<=nanmean(FirstLapNum)))./(1+size(FirstLapShuff,2));
    title(['p = ' num2str(p3)])
    set(gcf, 'Renderer', 'painters');    
    helper_savefig([dirs.figdir '\Final\FirstLap_Hist_' addlab])
    helper_saveandclosefig([dirs.figdir '\HoversMoving\FirstLap_Hist_' addlab])

    p = 1-binocdf(sum(pvalues<.05)-1,length(pvalues),.05);
    p2 = (sum(nanmean(SS)>=nanmean(RR))+1)./(sum(~isnan(nanmean(SS)))+1);
    disp(num2str([p p2 sum(pvalues<.05) length(pvalues) sum(pvalues<.05)./length(pvalues)]))
    figure; hold on; 
    histogram(nanmean(SS),'FaceColor','w','LineWidth',2); 
    yl = get(gca,'ylim'); 
    plot([nanmean(RR) nanmean(RR)],yl,'k','LineWidth',3)
    title(['p = ' num2str(round(p,2,'significant')) ', p = ' num2str(round(p2,2,'significant')) ', ' num2str(sum(pvalues<.05)) ' of ' num2str(length(pvalues)) ', '  num2str(round((sum(pvalues<.05)./length(pvalues))*100)) '%'])
    xlabel('Correlation of locations across slopes (r-value)')
    ylabel('Shuffle count')
    set(gca,'FontSize',26)
    set(gcf,'Position',[  2044          58         874         595])
    set(gcf, 'Renderer', 'painters');
    helper_savefig([dirs.figdir '\HoversMoving\Shuff_Hist_' addlab])
    helper_saveandclosefig([dirs.figdir '\Final\Shuff_Hist' addlab])
end
