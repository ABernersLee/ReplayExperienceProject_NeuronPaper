function brad_effect(basedir,plot_paper_figures,plot_other_figures)

%%%% params and dirs
load([basedir 'dirs_linear_lapbylap_addedIN.mat'],'dirs')
cd(dirs.spikedatadir)
d2 = dir('*.mat');
if ~isfolder([dirs.figdir '/Gamma/'])
    mkdir([dirs.figdir '/Gamma/'])
end
    
jd_cutoff = .4; 
wc_cutoff = .6;
numlaps = 30;
allmrv = []; allmrv2 = []; 
alldat = [];
numbins = 50; smoothf= 30;


%%%% get data out
for id = 1:size(d2,1)
    iident = d2(id).name;
    load(iident,'StepGamma','CandSeq','MidTime','Spike','params','CandCorr','CandDis','OutFR')
    if params.Novel==1
        CandPassCrit=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff;
        
        daydat = [];
        [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:end);max(Spike(:,2))]);
        GammaLaps = I(StepGamma(:,4));
        StepGamma(CandPassCrit(StepGamma(:,4))==0,1) = NaN;
        mrv = NaN(1,numlaps); mrv2 = mrv;
        for ilap = 1:numlaps
            dat = StepGamma(~isnan(StepGamma(:,1)) & GammaLaps==ilap,[1 3]);
            if ~isempty(dat)
                 [~,c,i] = histcounts(dat(:,1),-pi:(2*pi/numbins):pi);
                 mm = NaN(length(c)-1,1);
                for ideg = 1:length(c)-1; mm(ideg,1) = nanmean(dat(i==ideg,2)); end
                mrv(1,ilap) = circ_r(c(1:end-1)',mm);
                
                   mm = NaN(length(c)-1,1);
                for ideg = 1:length(c)-1; mm(ideg,1) = nanmedian(dat(i==ideg,2)); end                                         
                mrv2(1,ilap) = circ_r(c(1:end-1)',mm);               
                daydat = cat(1,daydat,dat); 
            end
        end
        allmrv = [allmrv;mrv];
        allmrv2 = [allmrv2;mrv2];
        
        
        if ~isempty(daydat) && plot_other_figures

            [mm2,c] = histcounts(daydat(:,1),-pi:(2*pi/numbins):pi);
            m2 = circ_mean(c(1:end-1)',mm2');
            daydat(:,1) = mod(daydat(:,1)-m2+2*pi,2*pi)-pi;

            [mm2,c,i] = histcounts(daydat(:,1),-pi:(2*pi/numbins):pi);
            mm = NaN(length(c)-1,1);
            for ideg = 1:length(c)-1; mm(ideg,1) = nanmedian(daydat(i==ideg,2)); end
            output2 = smoothts([mm2 mm2 mm2],'b',smoothf); 
            output = smoothts([mm;mm;mm]','b',smoothf);
            h = output(length(mm)+1:(2*length(mm)));
            h2 = output2(length(mm2)+1:(2*length(mm2)));
            hh = ([h h(1)]-min(h))./(range(h));
            hh2 = ([h2 h2(1)]-min(h2))./(range(h2));
            m = circ_mean(c(1:end-1)',mm);
            m2 = circ_mean(c(1:end-1)',mm2');
            figure; hold on;
            subplot(1,3,1)
            polarplot(c,hh,'b-',[m m],[0 1],'b-',c,hh2,'r-',[m2 m2],[0 1],'r-')
            subplot(1,3,2)
            bar(c,[h2 h2(1)])
            subplot(1,3,3)
            bar(c,[h h(1)])

            [pval2,z2] = circ_rtest(c(1:end-1)',mm2');

    %             [pval,z] = circ_rtest(c(1:end-1)',mm);
            [~,pval] = circ_corrcl(daydat(:,1),daydat(:,2));
            suptitle([iident(1:end-4) ', SpikeZ = ' num2str(round(z2,2,'significant')) ', p = ' num2str(round(pval2,2,'significant'))...
                ', StepRhoPval = ' num2str(round(pval,2,'significant'))])
    %             suptitle([ 'p = ' num2str(pval) ', Z = ' num2str(z) ', ' iident(1:end-4)])
            set(gcf,'Position',[680   506   899   472])
            helper_saveandclosefig([dirs.figdir '/Gamma/Brad_Effect_' iident(1:end-4)])
        end

        alldat = cat(1,alldat,daydat);
    end    
end

load(iident,'replayparams')



%%%%% plot brad effect in the same way
if plot_paper_figures
    numbins = 25; smoothf= 10;
    [mm2,c,i] = histcounts(alldat(:,1),-pi:(2*pi/numbins):pi); %spike prob
    [pval3,z3] = circ_rtest(mm2,1:numbins,(2*pi/numbins)); %spike prob
    mm = NaN(length(c)-1,1); %step size
    for ideg = 1:length(c)-1; mm(ideg,1) = ...
            nanmedian(alldat(i==ideg,2)*replayparams.Bin_Size); end
    [pval1,z1] = circ_rtest(mm,1:numbins,(2*pi/numbins)); %step size

    output2 = smoothts([mm2 mm2 mm2],'b',smoothf); 
    output = smoothts([mm;mm;mm]','b',smoothf);
    h = output(length(mm)+1:(2*length(mm)));
    h2 = output2(length(mm2)+1:(2*length(mm2)));
    hh = ([h h(1)]-min(h))./(range(h));
    hh2 = ([h2 h2(1)]-min(h2))./(range(h2));
    m = circ_mean(c(1:end-1)',mm);
    m2 = circ_mean(c(1:end-1)',mm2');
    close all
    pp = polarplot(c,hh,'b-',[m m],[0 1],'b-',c,hh2,'r-',[m2 m2],[0 1],...
        'r-','LineWidth',3);
    set(gca,'ThetaZeroLocation','top','ThetaDir','clockwise',...
        'ThetaTick',0:30:360,'RGrid','off','RTick',[],'FontSize',26);

    legend([pp(1) pp(3)],{['Step Size, N = ' num2str(size(alldat,1))...
        ', Z = ' num2str(z1) ', P = ' num2str(pval1)];['Spike Prob, N = '...
        num2str(size(alldat,1)) ', Z = ' num2str(z3) ', P = ' ...
        num2str(pval3)]},'Location','southoutside')
    set(gcf,'Position',[ 1907        -186        1055         774])
    helper_savefig([dirs.figdir 'Final\Gamma_Brad_Effect_Polar'])
    helper_saveandclosefig([dirs.figdir 'Gamma\Brad_Effect_Polar'])
end
