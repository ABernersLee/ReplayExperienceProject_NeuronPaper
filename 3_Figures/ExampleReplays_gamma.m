function ExampleReplays_gamma(basedir)

%%%%% found good replays and plot example 

wc_cutoff = .6; jd_cutoff = .4; coveragecutoff = .5;

load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')
fieldlab = 'LapbyLap';

if ~isfolder([dirs.figdir 'ExampleReplays\'])
    mkdir([dirs.figdir 'ExampleReplays\'])
end

cd(dirs.spikedatadir)
d2 = dir('*.mat');
EstBin = .02;


for ListNo = 38 %[35 36 38 41 43] %1:size(d2,1)
    cd(dirs.spikedatadir)
    load(d2(ListNo).name,'rawspikedata','hp_cells','hpinterneurons',...
        'Index','StepGamma','Radjusted','PosShuffleCorr',...
        'PosShuffleDis','CandStepSize','InMatrix','OutMatrix','Spike',...
        'CandCorr','CandDis','OutFR','CandSeq','CandRS','MidTime',...
        'params','hovertime2','movetime2','move2','hover2','-mat')        

    cellstouse = setdiff(hp_cells,hpinterneurons);
    Tetrode = unique(rawspikedata(ismember(rawspikedata(:,2),cellstouse),3));
    No = 1;
    cd(dirs.cscdatadir)
    
     if contains(d2(ListNo).name,'cpp') || contains(d2(ListNo).name,'sal')
         if Tetrode(No)<10
             EEGname=[d2(ListNo).name(1:6) 'tt0' num2str(Tetrode(No))...
                 d2(ListNo).name(6:end)];
         else
             EEGname=[d2(ListNo).name(1:6) 'tt' num2str(Tetrode(No))...
                 d2(ListNo).name(6:end)];
         end 
         load([EEGname],'csctimes','cscsamples','Fs','-mat')         

         if exist('Times','var')
             Fs = 3255;
             csctimes = Times; clear Times
             cscsamples = LFP_Samples; clear LFP_Samples
         end
     elseif contains(d2(ListNo).name,'Janni_2010408_Run') &&...
             str2double(d2(ListNo).name(end-4))>1
        EEGname=['Janni_tt' num2str(Tetrode(No)) d2(ListNo).name(6:end-4)]; 
         load([EEGname '.mat'],'csctimes','cscsamples','Fs')
     else
        EEGname=[d2(ListNo).name(1:end-4) '_csc_tt' num2str(Tetrode(No))];    
        if ~isfile([EEGname '.mat'])
            dEEG = dir([d2(ListNo).name(1:end-4)...
                '*_csc_tt' num2str(Tetrode(No)) '.mat']);
            if size(dEEG,1)>1
                csctimes_temp = []; cscsamples_temp = [];
                for de = 1:size(dEEG,1)
                   load(dEEG(de).name,'csctimes','cscsamples','Fs','-mat')   
                   csctimes_temp = cat(1,csctimes_temp,csctimes);
                   cscsamples_temp = cat(1,cscsamples_temp,cscsamples);
                end
                csctimes = csctimes_temp; cscsamples = cscsamples_temp;
            else
                 load(dEEG(1).name,'csctimes','cscsamples','Fs','-mat')  
            end
        else
            load([EEGname '.mat'],'csctimes','cscsamples','Fs')
        end
     end
%     end
    
    cd(dirs.spikedatadir)
    LFP_Frequency=Fs;
    LFP=cscsamples;
    LFPTime=csctimes;        

    [b,a]=butter(3,[25/LFP_Frequency*2 50/LFP_Frequency*2],'bandpass');
    LFP=filtfilt(b,a,LFP);
    
    D=OutMatrix'+InMatrix';

    [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:2:end);max(Spike(:,2))]);
    maxmove = CandDis/size(OutFR,2);
    pc = ((sum(abs(PosShuffleCorr)>=abs(CandCorr),2)+1)./(size(PosShuffleCorr,2)+1))<.05;
    pd = ((sum(PosShuffleDis<=CandDis,2)+1)./(size(PosShuffleDis,2)+1))<.05;
    CandPassCrit = abs(CandCorr)>wc_cutoff & maxmove<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff; % & pc & pd;
    I2 = I(CandPassCrit);                
    C = find(CandPassCrit);
    [~,II]=histc(LFPTime,[0:EstBin/4:Index(end)*EstBin]+CandSeq(1,1)); 
    
    for c=37 %1:length(C)
        DiffLoc = StepGamma(StepGamma(:,4)==C(c),3);
        G = StepGamma(StepGamma(:,4)==C(c),1);
        idat = II(ismember(II,4*Index(C(c))+1:4*Index(C(c)+1)-3));
        %microvolts to millivolts (mV)
        gdat = LFP(ismember(II,4*Index(C(c))+1:4*Index(C(c)+1)-3))/1000; 
        [~,locs1]=findpeaks(-gdat);
        locs = ((locs1-1)./length(idat))*range(idat)+1;

        f = figure; hold on
        subplot(2,1,1); hold on                        
        imagesc(D(:,4*Index(C(c))+1:4*Index(C(c)+1)-3)); colormap hot; axis xy; axis tight; 
        set(gca,'clim',[0 .1])
        xlabel('Timebin')
        ylabel('Decoded Position')
%         title(['Lap ' num2str(I(C(c))) ', ' num2str(C(c)) ', Slope: ' num2str(CandRS(C(c),2))...
%             ', StepD: ' num2str(round(CandStepSize(C(c),3,2),2,'significant')) ', NLin: ' ...
%             num2str(round(Radjusted(C(c),3)/Radjusted(C(c),1),2,'significant'))])
%          title(['Lap: ' num2str(I(C(c))) ', Slpe: ' num2str(CandRS(C(c),2))...
%                 ', StpD: ' num2str(round(CandStepSize(C(c),3,2),2,'significant')) ', Hvr: ' ...
%                 num2str(sum(hover2(hover2(:,4)==C(c),1))) ', Mve: ' num2str(sum(move2(move2(:,4)==C(c),1))) ...
%                 ', Stps: ' num2str(CandStepSize(C(c),2,2))]) % ', ' num2str(tm(C(c))) ' ms'])
%         title(['Lap: ' num2str(I(C(c))) ', Slpe: ' num2str(CandRS(C(c),2))...
%                         ', StpD: ' num2str(round(CandStepSize(C(c),3,2),2,'significant')) ', C: ' num2str(C(c)) ', c: ' num2str(c) ...
%                         ', Stps: ' num2str(CandStepSize(C(c),2,2))])
        yl = get(gca,'ylim');
        set(gca,'ytick',[1 yl(2)],'yticklabel',[1 size(D,1)*2.5])
        ylabel('Position (cm)')        
        set(gca,'xtick',[1 length(DiffLoc)+1],'xticklabel',[1 round((EstBin/4)*(length(DiffLoc)+1)*1000)])
        xlabel('Time (ms)')
        set(gca,'FontSize',20)

        subplot(2,1,2); hold on;
        
        line([locs locs]',[ones(size(locs))*[0 (max(DiffLoc)+1)]]','Color','b','LineStyle','--')        
        yyaxis left
        plot(DiffLoc,'k.-','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',20)       
        yl = get(gca,'ylim');
        ylim([yl(1) max(DiffLoc)+1])
        set(gca,'ytick',[1 yl(2)],'yticklabel',[1 round(yl(2)*2.5)])
        ylabel('Step size (cm)')
        set(gca,'xtick',[1 length(DiffLoc)],'xticklabel',[1 round((EstBin/4)*(length(DiffLoc)+1)*1000)])
        xlabel('Time (ms)')
        set(gca,'FontSize',20)
        axis tight
        
        yyaxis right
        plot((([1:length(idat)]-1)./length(idat))*range(idat)+1,gdat,'b')        
        yl = get(gca,'ylim');
        ylabel('Slow Gamma (mV)')
        set(gca,'xtick',[1 length(DiffLoc)],'xticklabel',[1 round((EstBin/4)*(length(DiffLoc)+1)*1000)])
        xlabel('Time (ms)')
        set(gca,'FontSize',20)
        axis tight
        
%         subplot(4,1,4); hold on
%         plot(G,'k.')
%         yl = get(gca,'ylim');
%         line([locs locs]',[ones(size(locs))*yl]','Color','k','LineStyle','--')
%         axis tight

        set(gcf,'Position',[ 2194        -131         697         866])
        helper_saveandclosefig([dirs.figdir 'ExampleReplays\' d2(ListNo).name(1:end-4) '_' num2str(C(c)) '_' fieldlab '_Gamma'])

    end    
end