function SigTestOverLaps_NonOverLapping(basedir)

%%%%%
%%%%% Figures of the significance tests of replays across paramters using
%%%%% replays with non-overlapping bins
%%%%%


%%

%%%%%
%%%%% directiors and get replays numbers out
%%%%%

load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')

if ~isfolder([dirs.figdir 'NonOverlappingLapByLap\'])
   mkdir([dirs.figdir 'NonOverlappingLapByLap\'])
end

cd(dirs.spikedatadir)
d2 = dir('*.mat');
JD = []; JDshuff = []; WC = []; WCshuff = [];lapI = []; idd = [];
for id = 1:size(d2,1)
    load(d2(id).name,'OutMatrix','MidTime','CandSeq','Spike','params','pos')
    if params.Novel==1        
        load(d2(id).name,'CandCorrNOL','TmShuffleCorrNOL','CandDisNOL','TmShuffleDisNOL') %was Pos not Tm
        PosShuffleDisNOL = TmShuffleDisNOL./size(OutMatrix,2);
        CandDisNOL = CandDisNOL./size(OutMatrix,2);
        JDshuff = cat(1,JDshuff,PosShuffleDisNOL);
        JD = cat(1,JD,CandDisNOL);
        WCshuff = cat(1,WCshuff,abs(TmShuffleCorrNOL));
        WC = cat(1,WC,abs(CandCorrNOL));
        [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[min(pos(:,1));MidTime(1:end);max(pos(:,1))]);
        lapI = cat(1,lapI,I);
        idd = cat(1,idd,ones(size(I))*id);
    end
end
%%

%%%%%
%%%%% plot for each combination of paramters wc and jd
%%%%%

n_shuffles = size(WCshuff,2);
lapnum = 16;
wcd = [.9:-.1:0];
jdd = .1:.1:1;
dat = NaN(length(wcd),length(jdd),lapnum); dat2 = dat; dat3 = dat;
quadreal = zeros(lapnum,1); quadshuff = zeros(n_shuffles,lapnum);
for ilap = 1:lapnum
    for wc = 1:length(wcd)
        for jd = 1:length(jdd)
            numinshuff = sum(JDshuff(lapI==ilap,:)<jdd(jd) & WCshuff(lapI==ilap,:)>wcd(wc));
            numreal = sum(JD<jdd(jd) & WC>wcd(wc) & lapI == ilap); 
            dat(wc,jd,ilap) = (sum(numinshuff>=numreal)+1)/(n_shuffles+1);        
            dat2(wc,jd,ilap) = numreal;
            dat3(wc,jd,ilap) = max(numinshuff);
            if wcd(wc)>=.6 && jdd(jd)<=.4
                quadreal(ilap) = quadreal(ilap)+numreal;
                quadshuff(:,ilap) = quadshuff(:,ilap)+numinshuff';
            end
        end
    end
    disp(ilap)
end

N = 64;
n = N/2;
cm = NaN(N,3);

cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)']; 
cm(:,3) = [linspace(0,1,n)';ones(N-n,1)]; 



sz(1) = round(sqrt(lapnum));
sz(2) = ceil(lapnum./sz(1));
figure; hold on; 
for ilap = 1:lapnum
    subplot(sz(1),sz(2),ilap); hold on
    imagesc(log10(dat(:,:,ilap))); axis xy; 

    set(gca,'ytick',1:3:length(wcd),'yticklabel',wcd(1:3:end))
    set(gca,'xtick',2:3:length(jdd),'xticklabel',jdd(2:3:end))
    xlabel('Abs(Max Jump Distance)')
    ylabel('Abs(Weighted Correlation)')
    set(gca,'clim',[log10(.05)*2 0])
    set(gcf,'colormap',cm)
    colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5])
    x = [.5 4.5 4.5 .5];
    y = [.5 .5 4.5 4.5] ; 
    patch(x,y,'w','FaceColor','none','EdgeColor','green','LineWidth',2)
    plot(3,3,'*g','MarkerSize',5)
    plot(4,3,'*g','MarkerSize',5)
    plot(4,4,'*g','MarkerSize',5)
    axis tight
    title(['Lap ' num2str(ilap-1)])
end
set(gcf,'Position',[ 1747        -137        1548         963])
helper_saveandclosefig([dirs.figdir 'NonOverlappingLapByLap\SigTestOverShortLaps2_Tm'])

ind = 0:.05:1;
hJDshuff = histc(JDshuff(:),ind);
hJD = histc(JD(:),ind);
hWCshuff = histc(WCshuff(:),ind);
hWC = histc(WC(:),ind);
figure; hold on;
subplot(1,2,1)
b = bar(ind,[hJDshuff./sum(hJDshuff) hJD./sum(hJD)]);
b(1).FaceColor = 'b'; b(2).FaceColor = 'r';
xlabel('abs(Max Jump Distance)')
ylabel('Histogram')
set(gca,'FontSize',18)
subplot(1,2,2)
c = bar(ind,[hWCshuff./sum(hWCshuff) hWC./sum(hWC)]);
c(1).FaceColor = 'b'; c(2).FaceColor = 'r';
xlabel('abs(Weighted Correlation)')
ylabel('Histogram')
legend(c,{'Shuffled','Real'})
set(gcf,'Position',[ 1966         480         814         323])
set(gca,'FontSize',18)
helper_saveandclosefig([dirs.figdir 'NonOverlappingLapByLap\SigTest_WCandJDdist_Tm'])

% figure; hold on;
quadp = (sum(quadshuff'>=quadreal,2)+1)./(size(quadshuff,1)+1);
figure; hold on
for iq = 1:length(quadp)
   text(iq/length(quadp),iq/length(quadp),['Lap ' num2str(iq) ', p = ' num2str(quadp(iq))])
end
helper_saveandclosefig([dirs.figdir 'NonOverlappingLapByLap\SigTest_QuadPvalues_Tm'])