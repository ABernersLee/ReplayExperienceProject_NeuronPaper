function ExampleReplays_Position(basedir)

%%%% found good replays and plot example sessions
%%%% figures that were used in the paper get saved out in Final

dirs.figdir = [basedir '/Figures/'];
dirs.spikedata = [basedir '/spikedata/'];
if ~isfolder([dirs.figdir '/Final/'])
    mkdir([dirs.figdir '/Final/'])
end

if ~isfolder([dirs.figdir '/ExampleReplays/'])
    mkdir([dirs.figdir '/ExampleReplays/'])
end
   
    
for itype = 1 
    if itype == 1
        savetag = 'lapbylap';
        cd([dirs.spikedata 'final_linear_lapbylap_addedIN'])
    elseif itype == 2
        savetag = 'wholesession';
        cd([dirs.spikedata 'final_linear_lapbylap_addedIN_wholesession\'])
    elseif itype >2
        savetag = ['downsampled_' num2str(itype)];
        cd([dirs.spikedata 'final_linear_lapbylap_addedIN_downsampled_sameevents\' num2str(itype) '\2\'])
    end



    d2 = dir('*.mat');
    
    for ilist = 1:5
        if ilist == 1
            ListNo = 38; C = [96 284 292 383 557 574]; %for ListNo=38;
        elseif ilist==2
            ListNo=43; C = [8 19 52 104 155 187 257]; %for ListNo=43;
        elseif ilist ==3
        ListNo=41; C = [28 453 471 533 569]; %for ListNo=41;
        elseif ilist == 4
            ListNo=35; C = [60 109 150 285]; %for ListNo=35;
        elseif ilist ==5
            ListNo=36; C = [7 62 108 159 213 242 268]; %for ListNo=36;
        end
        
        clear CandRS hover2
        load(d2(ListNo).name,'PosShuffleDis','PosShuffleCorr','hp_cells','hpinterneurons','pos','linposcat','Index',...
            'Radjusted','CandStepSize','hover2','move2','InMatrix','OutMatrix','Spike','CandCorr','CandDis','OutFR','InFR','CandSeq','CandRS','MidTime','params','hovertime2','movetime2','-mat')   
        if ~exist('hover2','var')
            continue
        end
        if ~exist('CandRS','var')
            CandRS = NaN(size(CandCorr,1),4);
        end
        [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:2:end);max(Spike(:,2))]);
        [~,Ipos]=histc(pos(:,1),[MidTime(1:2:end);max(Spike(:,2))]);
        maxmove = CandDis/size(OutFR,2);
        n_shuffles = size(PosShuffleDis,2);
        sig2 = (sum(PosShuffleDis<=CandDis,2)+1)./(n_shuffles+1);
        sig3 = (sum(PosShuffleCorr>=CandCorr,2)+1)./(n_shuffles+1);
        
        % I2 = I(CandPassCrit);        

        D=OutMatrix'+InMatrix';
        % C = find(CandPassCrit);/=
        % C = C([1 3 5:6 9:10 14:15 25 27 43 45]); %for fig
        % C = [8 19 47 52 60 77 100 104 155 187 310 322]; %new2 for ListNo=43;
        % C = [8 19 52 60 104 155 187 188 257 310]; %for ListNo=43;
        % C = [8 19 52 60 104 155 187 257 310]; %for ListNo=43;
        % C = [266 270 646]; %for ListNo=59;
        % C = [28 453 471 533 569]; %for ListNo=41;
        % C = [60 109 150 285]; %for ListNo=35;
        % C = [7 62 108 159 213 242 268]; %for ListNo=36;

        reptime = (CandSeq(C,1)+CandSeq(C,2))/2;    
        [~,b] = histc(reptime,pos(:,1));

        numover=4;
        sz(1) = ceil(length(C)/numover);
        if sz(1)==0
            sz(1)=length(C);
        end
        toadd = 1;
        tm = round((CandSeq(:,2)-CandSeq(:,1))*1000);
        
        f = figure; hold on
        
        for c=1:length(C)
            if numover+toadd>(numover*2)
                set(gcf,'Position',[ 396         505        1869         277])
                helper_saveandclosefig([dirs.figdir '\ExampleReplays\Position_' d2(ListNo).name(1:end-4) '_' num2str(c) '_' savetag])
                f = figure; hold on;
                toadd = 1;
            end
            subplot(1,numover,toadd); hold on
            toadd = toadd+1;
            imagesc(D(:,4*Index(C(c))+1:4*Index(C(c)+1)-3)); colormap hot; axis xy; axis tight; 
            set(gca,'clim',[0 .1])
            xlabel('Timebin')
            ylabel('Decoded Position')
            title(['Lap: ' num2str(I(C(c))) ', Slpe: ' num2str(CandRS(C(c),2))...
                ', StpD: ' num2str(round(CandStepSize(C(c),3,2),2,'significant')) ', Hvr: ' ...
                num2str(sum(hover2(hover2(:,4)==C(c),1))) ', Mve: ' num2str(sum(move2(move2(:,4)==C(c),1))) ...
                ', Stps: ' num2str(CandStepSize(C(c),2,2))]) % ', ' num2str(tm(C(c))) ' ms'])
            xlim([1 45])
            hdat = hovertime2(hovertime2(:,3)==C(c),1:2); 
            hl = line(hdat',zeros(size(hdat))','color','b','LineWidth',4);
            mdat = movetime2(movetime2(:,3)==C(c),1:2); 
            line(mdat',zeros(size(mdat))','color','r','LineWidth',4);
            axis off

        end
        
        if f.isvalid
            set(gcf,'Position',[ 396         505        1869         277])
            helper_saveandclosefig([dirs.figdir '\ExampleReplays\Position_' d2(ListNo).name(1:end-4) '_' num2str(c+1) '_' savetag])
        end
        
        if itype<3
            bbad1 = find(abs(CandCorr)<.2 & maxmove<.3 & CandStepSize(:,5,2)>.5);
            bbad2 = find(abs(CandCorr)>.7 & maxmove>.6 & CandStepSize(:,5,2)>.5);
            C2 = [C(end) bbad1(1) bbad2(1)];
            f = figure; hold on        
            for c=1:length(C2)        
                subplot(1,3,c); hold on            
                imagesc(D(:,4*Index(C2(c))+1:4*Index(C2(c)+1)-3)); colormap hot; axis xy; axis tight; 
                set(gca,'clim',[0 .1])
                xlabel('Timebin')
                ylabel('Decoded Position')
                set(gca,'FontSize',20)
                title({['Correlation: ' num2str(round(CandCorr(C2(c)),2,'significant'))];['Max Jump Dist: ' num2str(round(maxmove(C2(c)),2,'significant'))]},'FontWeight','normal') 
            end
            set(gcf,'Position',[1846         382        1202         378])
            helper_savefig([dirs.figdir '\Final\' savetag '_' d2(ListNo).name(1:end-4) '_good2badreplay'])            
            helper_saveandclosefig([dirs.figdir '\ExampleReplays\' savetag '_' d2(ListNo).name(1:end-4) '_good2badreplay'])
        end

        if itype==1

            figure; hold on; plot(pos(:,1),linposcat,'k','LineWidth',2)
            
            yl = get(gca,'ylim');
            mt = MidTime(1:end)+(diff([MidTime(1:end);max(pos(:,1))]))/8;
            for ilap = 3:2:length(mt)
                text(mt(ilap),yl(2)*1.01,num2str(ilap),'FontSize',12);
            end
            text(MidTime(1)+diff(MidTime(1:2))/8,yl(2)*1.03,'Lap','FontSize',12)
            plot(reptime,linposcat(b),'ro','MarkerSize',10,'LineWidth',3)
            
            yl(2) = yl(2)*1.05;
            plot([MidTime(1:end) MidTime(1:end)]',[ones(size(MidTime(1:end),1),1)*yl]','k--')
            
            for irep = 1:length(reptime)
                text(reptime(irep),linposcat(b(irep))+40,num2str(irep),'Color','r','FontSize',12)
            end
            set(gcf,'Position',[ 151         444        1621         296])
            xlim([min(pos(:,1)) max(pos(:,1))])
            xt = get(gca,'xtick');
            set(gca,'xticklabel',round((xt-min(pos(:,1)))/60,2,'significant'))
            set(gca,'ytick',[])
            set(gca,'FontSize',26)
            xlabel('Minutes')
            ylabel('Position')
            if itype<3
                helper_savefig([dirs.figdir '\Final\' savetag '_' d2(ListNo).name(1:end-4) '_Pos'])
            end
            helper_saveandclosefig([dirs.figdir '\ExampleReplays\' savetag '_' d2(ListNo).name(1:end-4) '_Pos'])


            FR = OutFR+InFR;
            [~,m] = max(FR,[],2);
            [~,ii] = sort(m);
            [~,iii] = sort(ii);
            hp_cells2 = NaN(size(hp_cells));
            hp_cells2(~ismember(hp_cells,hpinterneurons)) = iii;
            hpcells = hp_cells(~ismember(hp_cells,hpinterneurons));
            CellFR = NaN(length(hp_cells),length(C));
            dur = NaN(length(C),1);

            s = Spike(Spike(:,2)>CandSeq(C(3),1) & Spike(:,2)<CandSeq(C(3),2) & ismember(Spike(:,1),hpcells),1:2);
            cc = hp_cells2(s(:,1));
            numc = 3;
            rc = cc(2:round(length(cc)/numc):end);
            ccol = {'r';'g';'b';'m'};
            figure; hold on
            for c=1:length(C)
                s = Spike(Spike(:,2)>CandSeq(C(c),1) & Spike(:,2)<CandSeq(C(c),2) & ismember(Spike(:,1),hpcells),1:2);

                subplot(2,length(C),c); hold on    
                imagesc(D(:,4*Index(C(c))+1:4*Index(C(c)+1)-3)); colormap hot; axis xy; 
%                 xlim([1 45])
                xlim([-1 45]) %changed 6/8/2021 to line up (10ms on either side)
                set(gca,'clim',[0 .06])
                axis off
                title([num2str(tm(C(c))) ' ms'],'FontSize',26,'FontWeight','normal')


                subplot(2,length(C),c+length(C)); hold on
                xdat = [[s(:,2)-CandSeq(C(c),1)]*ones(1,2)]'; ydat = [hp_cells2(s(:,1))-.5 hp_cells2(s(:,1))+1.5]';
                l = line(xdat,ydat,'Color','k','LineWidth',3);
                for icell = 1:numc
                    l = line(xdat(:,ydat(1,:)+.5==rc(icell)),ydat(:,ydat(1,:)+.5==rc(icell)),'Color',ccol{icell},'LineWidth',3);
                end
                ylim([0 length(hpcells)+1]) 
                xlim([0 .225])
                set(gca,'ytick',[])

                CellFR(:,c) = histc(s(:,1),hp_cells);
                dur(c) = CandSeq(C(c),2)-CandSeq(C(c),1);

            end
            set(gcf,'Position',[ 641         180        2319         497])
            if ListNo == 41 || ListNo==35
             set(gcf,'Position',[ 1943         246        1013         559])
            end
            if itype<3
                helper_savefig([dirs.figdir '\Final\' savetag '_' d2(ListNo).name(1:end-4) '_Cells'])
            end
            helper_saveandclosefig([dirs.figdir '\ExampleReplays\' savetag '_' d2(ListNo).name(1:end-4) '_Cells'])
            


            figure; hold on
            for c=1:length(C)
                s = Spike(Spike(:,2)>CandSeq(C(c),1) & Spike(:,2)<CandSeq(C(c),2) & ismember(Spike(:,1),hpcells),1:2);

                subplot(3,length(C),c); hold on    
                imagesc(D(:,4*Index(C(c))+1:4*Index(C(c)+1)-3)); colormap hot; axis xy; 
                xlim([1 45])
                set(gca,'clim',[0 .06])
                axis off
                title(['Nlin: ' num2str(round(Radjusted(C(c),3)-Radjusted(C(c),1),2,'significant'))],'FontSize',12)

                Z = D(:,4*Index(C(c))+1:4*Index(C(c)+1)-3);
                I2=sum([1:size(Z,1)]'*ones(1,size(Z,2)).*Z);
                x = find(~isnan(I2));
                y = I2;

                pp = 1;
                subplot(3,length(C),c+length(C)); hold on
                [yy,S1] = polyfit(x,y,pp);
                [yyy] = polyval(yy,x);
                plot(x,y,'k.')
                plot(x,yyy,'r')        
                R2_1 = 1 - (S1.normr/norm(y - mean(y)))^2;
                Ra = 1-( (1- R2_1)*(length(x)-1)./(length(x)-pp-1));
                xlim([1 45])
                title([num2str(round(Ra,2,'significant'))],'FontSize',12)
                axis off

                pp = 3;
                subplot(3,length(C),c+length(C)*2); hold on
                [yy,S1] = polyfit(x,y,pp);
                [yyy] = polyval(yy,x);
                plot(x,y,'k.')
                plot(x,yyy,'r')    
                R2_1 = 1 - (S1.normr/norm(y - mean(y)))^2;
                Ra = 1-( (1- R2_1)*(length(x)-1)./(length(x)-pp-1));
                plot(x,y,'k.')
                plot(x,y,'k.')
                xlim([1 45])
                title([num2str(round(Ra,2,'significant'))],'FontSize',12)

                axis off

                CellFR(:,c) = histc(s(:,1),hp_cells);
                dur(c) = CandSeq(C(c),2)-CandSeq(C(c),1);

            end
            set(gcf,'Position',[ 641         180        2319         497])
            if ListNo == 41 || ListNo==35
             set(gcf,'Position',[ 1943         246        1013         559])
            end
            helper_saveandclosefig([dirs.figdir '\ExampleReplays\' savetag '_' d2(ListNo).name(1:end-4) '_CellsShape'])



            CellFR(ismember(hp_cells,hpinterneurons),:) = [];
            CellFR = CellFR';


            CellFR2 = CellFR;
            CellFR2(CellFR==0) = NaN;
            CellFR3 = CellFR>0;
            jit = rand(size(CellFR,2),1);
            jitt = repmat([(jit-.5)*.005]',[size(dur,1) 1]);
            jit = rand(size(CellFR,2),1);
            jitt2 = repmat([(jit-.8)*.5]',[size(dur,1) 1]);
            dur2 = repmat(dur,[1 size(CellFR,2)])+jitt;
            xl = [min(dur)-.01 max(dur)+.01];

            rc2 = rc;

            cnv = hp_cells(~ismember(hp_cells,hpinterneurons));
            for icell = 1:numc
                jnk = find(rc(icell)==hp_cells2,1,'first');
                rc2(icell) = find(cnv==jnk,1,'first');
            end

            figure; hold on;
            
            subplot(5,1,1); hold on
            plot(dur2,CellFR+jitt2,'.k')
            errorbar(dur+.005,mean(CellFR,2),std(CellFR,[],2)./sqrt(size(CellFR,1)),'k.','LineWidth',2)
            for icell = 1:numc
                plot(dur2(:,rc2(icell)),CellFR(:,rc2(icell))+jitt2(:,rc2(icell)),'.-','color',ccol{icell},'MarkerSize',10)
            end
            [r,p] = corr(dur,nanmean(CellFR,2),'type','Pearson');
            ylabel(['Spikes per event'])
            xlabel('Duration of event')
            set(gca,'FontSize',16)
            title(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])
            set(gca,'xlim',xl)


            subplot(5,1,2); hold on
            errorbar(dur,mean(CellFR3,2),std(CellFR3,[],2)./sqrt(size(CellFR3,1)),'k.','LineWidth',2)
            [r,p] = corr(dur,nanmean(CellFR3,2),'type','Pearson');
            ylabel(['Cells active'])
            xlabel('Duration of event')
            set(gca,'xlim',xl)
            set(gca,'FontSize',16)
            title(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])

            subplot(5,1,3); hold on
            plot(dur2,(CellFR+jitt2)./dur,'.k')
            for icell = 1:numc
                plot(dur2(:,rc2(icell)),(CellFR(:,rc2(icell))+jitt2(:,rc2(icell)))./dur,'.-','color',ccol{icell},'MarkerSize',10)
            end
            errorbar(dur+.005,mean(CellFR./dur,2),std(CellFR./dur,[],2)./sqrt(size(CellFR,1)),'k.','LineWidth',2)
            [r,p] = corr(dur,nanmean(CellFR./dur,2),'type','Pearson');
            ylabel(['FR per event'])
            xlabel('Duration of event')
            set(gca,'xlim',xl)
            set(gca,'FontSize',16)
            title(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])


            subplot(5,1,4); hold on
            plot(dur2,(CellFR2+jitt2)./dur,'.k')
            for icell = 1:numc
                nonandat = (CellFR2(:,rc2(icell))+jitt2(:,rc2(icell)))./dur;
                plot(dur2(~isnan(nonandat),rc2(icell)),nonandat(~isnan(nonandat)),'.-','color',ccol{icell},'MarkerSize',10)
            end
            errorbar(dur+.005,nanmean(CellFR2./dur,2),nanstd(CellFR2./dur,[],2)./sqrt(size(CellFR2,1)),'.k','LineWidth',2)
            [r,p] = corr(dur,nanmean(CellFR2./dur,2),'type','Pearson');
            ylabel(['FR of active cells'])
            xlabel('Duration of event')
            set(gca,'xlim',xl)
            set(gca,'FontSize',16)
            title(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])


            subplot(5,1,5); hold on
            errorbar(dur,mean(CellFR3,2)./dur,std(CellFR3./dur,[],2)./sqrt(size(CellFR3,1)),'k.','LineWidth',2)
            [r,p] = corr(dur,nanmean(CellFR3./dur,2),'type','Pearson');
            ylabel(['Cells active per second'])
            xlabel('Duration of event')
            set(gca,'xlim',xl)
            set(gca,'FontSize',16)
            title(['r = ' num2str(round(r,2,'significant')) ', p = ' num2str(round(p,2,'significant'))])
            set(gcf,'Position',[2172        -542         328        1278])
            helper_saveandclosefig([dirs.figdir '\ExampleReplays\' savetag '_' d2(ListNo).name(1:end-4) '_Cells_dur'])

        end
    end
end