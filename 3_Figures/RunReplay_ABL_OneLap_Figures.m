function RunReplay_ABL_OneLap_Figures(basedir,plot_paper_figures,plot_other_figures)

%%%%%
%%%%%
%%%%% A lot of different plots of the replays after one lap of running.
%%%%% Including example sessions, group statistics of changes over laps and
%%%%% number of replays in each epoch and significanve
%%%%%
%%%%%

%%

%%%%%
%%%%% dirs
%%%%%
load([basedir 'dirs_linear_OneLap'],'dirs')

if ~isfolder(dirs.figdir)
    mkdir(dirs.figdir)
end
if ~isfolder([dirs.figdir '\Paper\'])
    mkdir([dirs.figdir '\Paper\'])
end
%% 

%%%%%
%%%%% plots fields and behavior and replay example sessions
%%%%%

allCandDis0 = []; allCandDis1 = [];
cd(dirs.spikedatadir)
d2 = dir('*.mat');
seslab = {'Sleep0';'Run1';'Sleep1';'Run2';'Sleep2'};
jd_cutoff = .4; wc_cutoff = .6; coveragecutoff = 0;
for id = 1:size(d2,1)

    load(d2(id).name,'InFR','OutFR','params','pos','linposcat','dirdat',...
        'vel','FREpochOne','PosShuffleCorrEpoch1','PosShuffleDisEpoch1',...
        'CandCorrEpoch1','CandDisEpoch1','CandSeq','MidTime','Spike',...
        'replayparams','params','MatrixEpoch1','Index','CandStepSizeE1')
    
    OutFR = FREpochOne; clear InFR
    t= CandStepSizeE1(:,5,2)>.4; 
    tt=CandDisEpoch1/size(OutFR,2)<jd_cutoff & abs(CandCorrEpoch1)>wc_cutoff;
    ts0 = CandSeq(:,1)>params.Sleep_Times(1,1) & CandSeq(:,1)<params.Sleep_Times(1,2);
    t1 = CandSeq(:,1)>params.Run_Times(1,1) & CandSeq(:,1)<params.Run_Times(1,2);
    ts1 = CandSeq(:,1)>params.Sleep_Times(2,1) & CandSeq(:,1)<params.Sleep_Times(2,2);
    t2 = CandSeq(:,1)>params.Run_Times(2,1) & CandSeq(:,1)<params.Run_Times(2,2);
    ts2 = CandSeq(:,1)>params.Sleep_Times(3,1) & CandSeq(:,1)<params.Sleep_Times(3,2);
    ts = cat(2,ts0 &t,t1 &t,ts1 &t,t2 &t,ts2 &t);

    if plot_paper_figures
       figure; hold on
       irun = 1;    
       plot(pos(pos(:,1)>params.Run_Times(irun,1) & pos(:,1)<params.Run_Times(irun,2),1) - params.Run_Times(irun,1),linposcat(pos(:,1)>params.Run_Times(irun,1) & pos(:,1)<params.Run_Times(irun,2)),'.k')   
       irun = 2;
       plot(pos(pos(:,1)>params.Run_Times(irun,1) & pos(:,1)<params.Run_Times(irun,2),1) + 3*range(params.Run_Times(1,1:2)) - params.Run_Times(irun,1),linposcat(pos(:,1)>params.Run_Times(irun,1) & pos(:,1)<params.Run_Times(irun,2)),'.k')       
       set(gcf,'position',[ 680   553   997   425])   
       helper_saveandclosefig([dirs.figdir '\Paper\behavior2_' params.ident])
    end


   allCandDis0 = cat(1,allCandDis0,[CandDisEpoch1(ts0)./size(OutFR,2) CandCorrEpoch1(ts0) CandStepSizeE1(ts0,5,2)]);
   allCandDis1 = cat(1,allCandDis1,[CandDisEpoch1(ts1)./size(OutFR,2) CandCorrEpoch1(ts1) CandStepSizeE1(ts1,5,2)]);


  for isession = 1:size(ts,2)      
      rep = find(ts(:,isession));
      [~,o] = sort(abs(CandCorrEpoch1(rep)) .* (1-((CandDisEpoch1(rep)/size(OutFR,2)))),'descend');
      rep2 = rep(o);
      r2 = abs(CandCorrEpoch1(rep));
      r = r2(o);
      dd = ((CandDisEpoch1(rep)/size(OutFR,2)));
      d = dd(o);


      if plot_paper_figures
          figure; hold on
          for i = 1:25       
              if length(rep2)>=i
                  Mat = MatrixEpoch1(Index(rep2(i))+1:Index(rep2(i)+1),:);
                  subplot(5,5,i)
                  imagesc(Mat'); colormap hot; axis xy; axis tight
                  set(gca,'clim',[0 .06])              
                  xlabel(['r = ' num2str(round(r(i),2,'significant')) ', d = ' num2str(round(d(i),2,'significant'))])
    %               set(gca,'ytick',[],'xtick',[])
              end
          end
          suptitle(seslab{isession})
          set(gcf,'renderer','Painters')
          set(gcf,'Position',[  680   349   615   629])
          helper_saveandclosefig([dirs.figdir 'Paper\ReplayExamples_sorted_' seslab{isession} '_' params.ident '_ax'])     
      end
      
  end

  if id == 1 && plot_paper_figures
      figure; hold on
      imagesc(Mat'); colormap hot; axis xy; axis tight
      set(gca,'clim',[0 .06])            
      c = colorbar('location','eastoutside','Ticks',[.01 .03 .05],'TickLength',.06,'TickDirection','both');
      set(gcf,'position',[2183         371         533         139])
       helper_saveandclosefig([dirs.figdir 'Paper\ReplayExamples_Colorbar_' params.ident])
  end
    
end


%%

%%%%%
%%%%% plots duration and slope across presleep, first experience, 
%%%%% sleep, second experience, post-sleep,  using run1 or 2 fields
%%%%%

if plot_paper_figures
    fieldlab = {'Run1Fields';'Run2Fields';'BothFields'};
    shufftoo = false;
    shufftoolab = {'_overlap';'_Shufftoo'};
    ifield=3;    
    jd_cutoff = .4; wc_cutoff = .6; coveragecutoff = 0;
    cd(dirs.spikedatadir)
    d2 = dir('*.mat');
    for mlaptimeskip = [1 4] %1 for the stats, 4 for the figure
        clearvars -except d2 jd_cutoff wc_cutoff coveragecutoff plot_paper_figures plot_other_figures dirs fieldlab shufftoo shufftoolab ifield mlaptimeskip
        cd(dirs.spikedatadir)
        
        mlaptimeskip2 =mlaptimeskip;
        sleep2replays = []; sleep0replays = []; sleep1replays = []; run2replays = []; run1replays = []; session_rat = [];
        numlaps_epoch = NaN(size(d2,1),5);
        for idd = 1:size(d2,1)
    % 
            if ifield == 2
                load(d2(idd).name,'TmShuffleCorr','TmShuffleDis','CandCorr_overlap','OutFR','CandDis_overlap','CandRS_overlap', ...
                    'CandSeq_overlap','MidTime','Spike','CandStepSize','replayparams_overlap','params')
                PosShuffleCorr = TmShuffleCorr;
                PosShuffleDis = TmShuffleDis;
                CandCorr = CandCorr_overlap;
                CandDis = CandDis_overlap;
                CandRS = CandRS_overlap;
                CandSeq = CandSeq_overlap;
                CandStepSize = ones(size(CandDis,1),5,2);

            elseif ifield == 3 || ifield==1
                load(d2(idd).name,'TmShuffleCorrEpoch1','TmShuffleDisEpoch1','CandCorrEpoch1_overlap','FREpochOne_overlap',...
                    'CandDisEpoch1_overlap','CandRSEpoch1_overlap','CandSeq_overlap','MidTime','Spike','CandStepSizeE1','replayparams_overlap','params')
               PosShuffleCorr = TmShuffleCorrEpoch1;
                PosShuffleDis = TmShuffleDisEpoch1;
                CandCorr = CandCorrEpoch1_overlap;
                OutFR = FREpochOne_overlap;
                CandDis = CandDisEpoch1_overlap;
                CandRS = CandRSEpoch1_overlap;
                CandStepSize = CandStepSizeE1;
                CandSeq = CandSeq_overlap;
            end
            if shufftoo
                sigst1 = (sum(((PosShuffleDis./size(OutFR,2)))<=CandDis/size(OutFR,2),2)+1)./(size(PosShuffleCorr,2)+1);
                sigst2 = (sum(abs(PosShuffleCorr)>=abs(CandCorr),2)+1)./(size(PosShuffleCorr,2)+1);
            else
                sigst1 = zeros(size(CandDis));
                sigst2 = sigst1;
            end


            mlaptime = median(diff(MidTime));

                 %1/26/2022 change* need to change unites of slope!! 5*old#     *replayparams_overlap.Bin_Size/100*1000/5
                 
            %pre-epoch sleep (sleep0)
            MidSleep = params.Sleep_Times(1,1):mlaptime:params.Sleep_Times(2,2);
            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidSleep(1:mlaptimeskip:end)]);
            numlaps_epoch(idd,1) =max(I);                     
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff & CandSeq(:,1)>params.Sleep_Times(1,1) & CandSeq(:,1)<params.Sleep_Times(1,2) & sigst1<.05 & sigst2<.05;
           sleepreplays00 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams_overlap.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams_overlap.Bin_Size CandStepSize(t,4,2)*replayparams_overlap.Bin_Size CandSeq(t,1)-params.Sleep_Times(1,1)];

            %get any run1 replays
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff & CandSeq(:,1)>params.Run_Times(1,1) & CandSeq(:,1)<params.Run_Times(1,2) & sigst1<.05 & sigst2<.05;
           run1replays1 = [ones(sum(t),1) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams_overlap.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams_overlap.Bin_Size CandStepSize(t,4,2)*replayparams_overlap.Bin_Size CandSeq(t,1)-params.Run_Times(1,1)];
            numlaps_epoch(idd,2) =NaN;  

            %get sleep 1 (btwn run1 and run2)  lap-equivilants
            endsleep1time = params.Sleep_Times(2,2); 
            MidSleep = params.Sleep_Times(2,1):mlaptime:endsleep1time;
            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidSleep(1:mlaptimeskip:end)]);
            numlaps_epoch(idd,3) =max(I);
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff & CandSeq(:,1)>params.Sleep_Times(2,1) & CandSeq(:,1)<params.Sleep_Times(2,2) & sigst1<.05 & sigst2<.05;
            sleepreplays11 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams_overlap.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams_overlap.Bin_Size CandStepSize(t,4,2)*replayparams_overlap.Bin_Size CandSeq(t,1)-params.Sleep_Times(2,1)];


            if ifield == 2 || ifield==3
                load(d2(idd).name,'TmShuffleCorr','TmShuffleDis','CandCorr_overlap','OutFR','CandDis_overlap','CandRS_overlap', ...
                    'CandSeq_overlap','MidTime','Spike','CandStepSize','replayparams_overlap','params')
                PosShuffleCorr = TmShuffleCorr;
                PosShuffleDis = TmShuffleDis;
                CandCorr = CandCorr_overlap;
                CandDis = CandDis_overlap;
                CandRS = CandRS_overlap;
                CandSeq = CandSeq_overlap;
                CandStepSize = ones(size(CandDis,1),5,2);

            elseif ifield==1
                load(d2(idd).name,'TmShuffleCorrEpoch1','TmShuffleDisEpoch1','CandCorrEpoch1_overlap','FREpochOne_overlap',...
                    'CandDisEpoch1_overlap','CandRSEpoch1_overlap','CandSeq_overlap','MidTime','Spike','CandStepSizeE1','replayparams_overlap','params')
               PosShuffleCorr = TmShuffleCorrEpoch1;
                PosShuffleDis = TmShuffleDisEpoch1;
                CandCorr = CandCorrEpoch1_overlap;
                OutFR = FREpochOne_overlap;
                CandDis = CandDisEpoch1_overlap;
                CandRS = CandRSEpoch1_overlap; 
                CandStepSize = CandStepSizeE1;
                CandSeq = CandSeq_overlap;
            end

            if shufftoo
                sigst1 = (sum(((PosShuffleDis./size(OutFR,2)))<=CandDis/size(OutFR,2),2)+1)./(size(PosShuffleCorr,2)+1);
                sigst2 = (sum(abs(PosShuffleCorr)>=abs(CandCorr),2)+1)./(size(PosShuffleCorr,2)+1);
            else
                sigst1 = zeros(size(CandDis));
                sigst2 = sigst1;
            end

            %get run2 replays
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff & CandSeq(:,1)>params.Run_Times(2,1) & CandSeq(:,1)<params.Run_Times(2,2) & sigst1<.05 & sigst2<.05;

            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime(1:mlaptimeskip2:end);max(Spike(:,2))]);
            numlaps_epoch(idd,4) =max(I);  

                 %1/26/2022 change* need to change unites of slope!! 5*old#     *replayparams_overlap.Bin_Size/100*1000/5
            I(I>length(MidTime)) = length(MidTime);
            lapreplays1 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams_overlap.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams_overlap.Bin_Size CandStepSize(t,4,2)*replayparams_overlap.Bin_Size CandSeq(t,1)-params.Run_Times(2,1)];

            %run post epoch2 (sleep2)
            MidSleep = params.Sleep_Times(3,1):mlaptime:params.Sleep_Times(3,2);
            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidSleep(1:mlaptimeskip:end)]);
            numlaps_epoch(idd,5) =max(I);    
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,2)>coveragecutoff & CandSeq(:,1)>params.Sleep_Times(3,1) & CandSeq(:,1)<params.Sleep_Times(3,2) & sigst1<.05 & sigst2<.05;
            sleepreplays22 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams_overlap.Bin_Size/100*1000/5 CandStepSize(t,3,2)*replayparams_overlap.Bin_Size CandStepSize(t,4,2)*replayparams_overlap.Bin_Size CandSeq(t,1)-params.Sleep_Times(3,1)];

            sleep0replays = cat(1,sleep0replays,sleepreplays00);
            run1replays = cat(1,run1replays,run1replays1);
            sleep1replays = cat(1,sleep1replays,sleepreplays11);    
            run2replays = cat(1,run2replays,lapreplays1);
            sleep2replays = cat(1,sleep2replays,sleepreplays22);
            session_rat = cat(1,session_rat,idd*ones(size(sleepreplays11,1),1));
        end

      
        %cut long sessions
        cutofflaps = median(numlaps_epoch);
        sleep0replays(sleep0replays(:,1)>cutofflaps(1),:) = [];        
        session_rat(sleep1replays(:,1)>cutofflaps(3),:) = [];
        sleep1replays(sleep1replays(:,1)>cutofflaps(3),:) = [];
        run2replays(run2replays(:,1)>cutofflaps(4),:) = [];
        sleep2replays(sleep2replays(:,1)>cutofflaps(5),:) = [];
    

       
        % Duration and Slope new plot
        tlab = {'Duration';'Slope';'rangeDtime'; 'mean(abs(moves))'};
        for itype = 2:3 
            xt = [1]; xtl = [1];
            
            figure; hold on;
            numnow = 1;

            %run1
            errorbar(numnow,mean(run1replays(:,itype)),std(run1replays(:,itype))./sqrt(size(run1replays,1)),'k','LineWidth',3)
            numnow = numnow+8;

            %sleep1
            A = NaN(max(sleep1replays(:,1)),2);
            for ilap = 1:max(sleep1replays(:,1))
                A(ilap,1) = nanmean(sleep1replays(sleep1replays(:,1)==ilap,itype));
                A(ilap,2) = nanstd(sleep1replays(sleep1replays(:,1)==ilap,itype))./sqrt(sum(sleep1replays(:,1)==ilap));
            end
            A(isnan(A(:,1)),:) = [];

              Filter=fspecial('gaussian',[10 1],4);
                pd = 10;
                A0 = A;
                A2 = filter2(Filter,[A0(1:pd,1);A0(:,1);A0(end-(pd-1):end,1)]);
                A3 = filter2(Filter,[A0(1:pd,2);A0(:,2);A0(end-(pd-1):end,2)]);
                A2 = A2(pd+1:end-pd,1);
                A3 = A3(pd+1:end-pd,1);
                fwd = [A2-A3]'; rev = [A2+A3]';
                p2 = patch([numnow:size(A0,1)+numnow-1 size(A0,1)+numnow-1:-1:numnow],[fwd rev(end:-1:1)],'black');
                p2.FaceAlpha=.1; p2.EdgeAlpha=0;
                plot(numnow:size(A0,1)+numnow-1,A2,'k','LineWidth',3)   
                xtn = get(gca,'xtick');
            xt = [xt numnow:median(diff(xtn)):size(A0,1)+numnow-1];
            xtl = [xtl [1:median(diff(xtn)):size(A0,1)+numnow-1]*mlaptimeskip];
                
             
            numnow = size(A0,1)+numnow+8;
                
            %run2
            B = NaN(max(run2replays(:,1)),2);
            for ilap = 1:max(run2replays(:,1))
                B(ilap,1) = nanmean(run2replays(run2replays(:,1)==ilap,itype));
                B(ilap,2) = nanstd(run2replays(run2replays(:,1)==ilap,itype))./sqrt(sum(run2replays(:,1)==ilap));
            end
            B(isnan(B(:,1)),:) = [];
            errorbar(numnow:size(B,1)+numnow-1,B(:,1),B(:,2),'k','LineWidth',3)
            xt = [xt numnow:2:size(B,1)+numnow-1];
            xtl = [xtl [1:2:size(B,1)]*mlaptimeskip2];
            numnow = numnow+size(B,1)+8;
            ylabel(tlab{itype-1})
            set(gca,'FontSize',18)

            if itype==2
                [s2r,s2p] = corr(sleep2replays(:,1),sleep2replays(:,itype),'rows','complete','tail','right');
                [s1r,s1p] = corr(sleep1replays(:,itype),sleep1replays(:,1),'rows','complete','tail','right');
                [r2r,r2p] = corr(run2replays(:,itype),run2replays(:,1),'rows','complete','tail','right');
            elseif itype==3                
                [s2r,s2p] = corr(sleep2replays(:,1),sleep2replays(:,itype),'rows','complete','tail','left');
                [s1r,s1p] = corr(sleep1replays(:,itype),sleep1replays(:,1),'rows','complete','tail','left');
                [r2r,r2p] = corr(run2replays(:,itype),run2replays(:,1),'rows','complete','tail','left');
            end
            [p3,tbl] = anovan([sleep1replays(:,itype); run2replays(:,itype)],[[ones(size(sleep1replays,1),1)*2; ones(size(run2replays,1),1)*3] [sleep1replays(:,1); run2replays(:,1)]],'model','interaction','continuous',2,'display','off');
            p6 = anovan([sleep2replays(:,itype); run2replays(:,itype)],[[ones(size(sleep2replays,1),1)*2; ones(size(run2replays,1),1)*3] [sleep2replays(:,1); run2replays(:,1)]],'model','interaction','continuous',2,'display','off');
            
            session_rat(:,2) = NaN(size(session_rat(:,1)));
            session_rat(session_rat(:,1)==1,2) = 1;
            session_rat(session_rat(:,1)>1 & session_rat(:,1)<4 ,2) = 2;
            session_rat(session_rat(:,1)>3 ,2) = 3;
            [rSR,pSR] = anovan([sleep1replays(:,itype)],[sleep1replays(:,1) session_rat(:,1:2)],'model','interaction','continuous',1,'nested',[[0 1 0];[0 0 1];[0 0 0]]);
            
            set(gca,'xlim',[-4 numnow+4])
            set(gcf,'Position',[  1934         298         983         326])
            set(gcf, 'Renderer', 'painters');
             title(['n1 = ' num2str(sum(~isnan(sleep1replays(:,1,1)))) ' r='  num2str(round(s1r,2,'significant')) ' p='  num2str(round(s1p,2,'significant')) ...
            ', n2=' num2str(sum(~isnan(run2replays(:,1,1)))) ' r='  num2str(round(r2r,2,'significant')) ' p='  num2str(round(r2p,2,'significant'))])
            set(gca,'xtick',xt,'xticklabel',xtl)
            helper_saveandclosefig([dirs.figdir '\Paper\NonOverlapAcrossEpochs_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_' fieldlab{ifield} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1}])


            figure; hold on;
            plab = {'Group';'Lap';'Interaction'};
            for ip = 1:3
                text(0,ip/4,[plab{ip} ': F = ' num2str(tbl{ip+1,6}) ' df = ' num2str(tbl{ip+1,3}) ',' num2str(sum(~isnan([sleep1replays(:,itype); run2replays(:,itype)]))-tbl{ip+1,3}) ' , p = ' num2str(p3(ip))])
            end
            set(gcf,'Position',[1919         133         861         403])
            axis off
            helper_saveandclosefig([dirs.figdir '\Paper\NonOverlapAcrossEpochs_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_ANOVAstats_' tlab{itype-1} '_' fieldlab{ifield} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1}])

            disp('Correlations (r,p):')        
            disp(round([s1r s1p; r2r r2p; s2r s2p] ,2,'significant'))
            disp('Interaction:')
            disp(round([p3(3) p6(3)] ,2,'significant'))

            totest1 = sleep1replays(~isnan(sleep1replays(:,1)),:);
            totest2 = run2replays(~isnan(run2replays(:,1)),:);
            totest = [totest1; totest2];
            

            %change point analysis
            ipt = findchangepts(totest(:,itype),'MaxNumChanges',1);
            if isempty(ipt)
                ipt = findchangepts(totest(:,itype),'MaxNumChanges',2);
            end
            
            if plot_other_figures && ~isempty(ipt)
                figure; hold on
                for ip = 1:length(ipt)
                    if ipt(ip)<=size(totest1,1)
                        thelap(ip,1) = totest1(ipt(ip),1);
                        oflaps(ip,1) = max(totest1(:,1));
                        ws(ip,1) = 1;
                        text(.1,.1*ip,['Change point in Sleep 1 on "lap" ' num2str(thelap(ip,1)) ' of ' num2str(oflaps(ip,1))])
                    else
                        thelap(ip,1) = totest2(ipt(ip)-size(totest1,1),1);
                        oflaps(ip,1) = max(totest2(:,1));
                        ws(ip,1) = 2;
                        text(.1,.1*ip,['Change point in Run 2 on lap ' num2str(thelap(ip,1)) ' of ' num2str(oflaps(ip,1))])
                    end

                end
                helper_saveandclosefig([dirs.figdir 'ChangePointAnalysis_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_' fieldlab{ifield} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1}])
            end
            
            %Chow test
                %         Chow tests rely on:
                % 
                % Independent, Gaussian-distributed innovations
                % 
                % Constancy of the innovations variance within subsamples
                % 
                % Constancy of the innovations across any structural breaks

%             Mdl = fitlm(totest(:,1),totest(:,itype));
%             res = Mdl.Residuals.Raw;
%             %should be false
%             [hARCH,pValueARCH] = archtest(res); 
%             % should be false, seems to be not a normal distribution so do a permutation test too
%             [hKS,pValueKS] = kstest(res/std(res)); 
            
            [~,p1,teststat] = chowtest(totest(:,1),totest(:,itype),size(totest1,1)+1,'Coeffs',[0 1]); %,'Display','summary');
            numperm = 1000; permstat = NaN(numperm,1);
            for ip = 1:numperm            
                [~,~,permstat(ip)] = chowtest(totest(:,1),totest(randperm(size(totest,1)),itype),size(totest1,1)+1,'Coeffs',[0 1]);
            end
            p = (1+(sum(permstat>=teststat)))./(1+numperm);
            
            figure; hold on; 
            histogram(permstat,'FaceColor','none','LineWidth',2); 
            yl = get(gca,'ylim'); plot([teststat teststat],yl,'-r','LineWidth',3)
            xl = get(gca,'xlim');
            text(xl(2)*.5,yl(2)*.8,['F = ' num2str(teststat) ', N = ' num2str(size(totest,1)) ', bp = ' num2str(size(totest1,1)+1)])
            ylabel('Permutations')
            xlabel('Chow Test''s F')
            set(gca,'FontSize',18)
            title({['Permutation test on chow test statistic'],['Lap and ' tlab{itype-1} ', OGp = ' num2str(round(p1,2,'significant')) ', PERMp = ' num2str(round(p,2,'significant'))]})
            set(gcf, 'Renderer', 'painters');
            set(gcf,'Position',[2156        -676         666         544])
            helper_saveandclosefig([dirs.figdir '\Paper\ChowTest_Lap_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_' fieldlab{ifield} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1}])
            disp(p)
            
            %BY TIME
%            Mdl = fitlm(totest(:,end),totest(:,itype));
%             res = Mdl.Residuals.Raw;
%             %should be false
%             [hARCH,pValueARCH] = archtest(res); 
%             % should be false, seems to be not a normal distribution so do a permutation test too
%             [hKS,pValueKS] = kstest(res/std(res)); 
            
            [~,p1,teststat] = chowtest(totest(:,end),totest(:,itype),size(totest1,1)+1,'Coeffs',[0 1]); %,'Display','summary');
            numperm = 1000; permstat = NaN(numperm,1);
            for ip = 1:numperm            
                [~,~,permstat(ip)] = chowtest(totest(:,end),totest(randperm(size(totest,1)),itype),size(totest1,1)+1,'Coeffs',[0 1]);
            end
            p = (1+(sum(permstat>=teststat)))./(1+numperm);

            figure; hold on; histogram(permstat,'FaceColor','none','LineWidth',2); 
            yl = get(gca,'ylim'); plot([teststat teststat],yl,'-r','LineWidth',3)
            xl = get(gca,'xlim');
            text(xl(2)*.5,yl(2)*.8,['F = ' num2str(teststat) ', N = ' num2str(size(totest,1)) ', bp = ' num2str(size(totest1,1)+1)])
            ylabel('Permutations')
            xlabel('Chow Test''s F')
            set(gca,'FontSize',18)        
            title({['Permutation test on chow test statistic'],['Time and ' tlab{itype-1} ', OGp = ' num2str(round(p1,2,'significant')) ', PERMp = ' num2str(round(p,2,'significant'))]})
            set(gcf, 'Renderer', 'painters');
            set(gcf,'Position',[2156        -676         666         544])
            helper_saveandclosefig([dirs.figdir '\Paper\ChowTest_Time_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_' fieldlab{ifield} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1}])
            disp(p)

            %change between sleep and start of run 2
            figure; hold on
            text(.1,.9,tlab{itype-1})
            if sum(totest2(:,1)==1)>0
                p = ranksum(totest1(:,itype),totest2(totest2(:,1)==1,itype));
                text(.1,.8,['Ranksum between Sleep1 and Lap 1 of run 2, p = ' num2str(p)])
                text(.1,.7,['Sleep 1 N = ' num2str(sum(~isnan(totest1(:,itype)))) ', Run 2 N = ' num2str(sum(~isnan(totest2(totest2(:,1)==1,itype))))])
            end
            p = ranksum(totest1(:,itype),totest2(totest2(:,1)<5,itype));
            text(.1,.6,['Ranksum between Sleep1 and Laps 1-4 of run 2, p = ' num2str(p)])
            text(.1,.5,['Sleep 1 N = ' num2str(sum(~isnan(totest1(:,itype)))) ', Run 2 N = ' num2str(sum(~isnan(totest2(totest2(:,1)<5,itype))))])
            helper_saveandclosefig([dirs.figdir '\Paper\Change_SleepRun_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_' fieldlab{ifield} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1}])

        end    
    end    
    
end
%%

%%%%%
%%%%% checking to see if there is a significant number of events in 
%%%%% sleep using run1 or 2 fields USING TIMEBIN SHUFFLE
%%%%%

if plot_other_figures || plot_paper_figures
  
    fieldlab = {'Run1 Fields';'Run2 Fields'};
    controllab = {'';'NumControl';'TimeControl'}; %For Neuron reviewers
    wcd = [.9:-.1:0];
    jdd = .1:.1:1;
    cd(dirs.spikedatadir)
    
     d2 = dir('*.mat');
     star = [];
     quadreal = zeros(6,2,3);
     quadshuff = zeros(6,5000,2,3);
     real_all = zeros(length(wcd),length(jdd),6,2);
     sh_all = zeros(length(wcd),length(jdd),6,2,5000);
     
     real_all_controls = zeros(length(wcd),length(jdd),6,2,2);
     sh_all_controls = zeros(length(wcd),length(jdd),6,2,5000,2);
     
    for idtest = 1:size(d2,1)
        for icontrols = 1:3
            dat = NaN(length(wcd),length(jdd),6,2);
            for wc = 1:length(wcd)
                for jd = 1:length(jdd)            
                    for ifield = 1
                        clearvars -except ifield fieldlab wc jd dat wcd jdd idtest star d2 quadreal quadshuff ...
                            plot_paper_figures plot_other_figures dirs real_all sh_all icontrols real_all_controls ...
                            sh_all_controls controllab

                        cd(dirs.spikedatadir)            
                        jd_cutoff = jdd(jd); wc_cutoff = wcd(wc);            
                        allreal = []; allshuff = []; alltotals = []; allsig= [];

                        for idd = idtest 

                            if ifield==1
                                load(d2(idd).name,'Index','TmShuffleCorrEpoch1','TmShuffleDisEpoch1','CandCorrEpoch1','FREpochOne','CandDisEpoch1','CandSeq','MidTime','params')
                                len = diff(Index); candind = len>=0;
                                TmShuffleCorr = TmShuffleCorrEpoch1(candind,:);
                                TmShuffleDis = TmShuffleDisEpoch1(candind,:);
                                CandCorr = CandCorrEpoch1(candind);
                                OutFR = FREpochOne;
                                CandDis = CandDisEpoch1(candind);
                                CandSeq = CandSeq(candind,:);
                            elseif ifield == 2
                                load(d2(idd).name,'TmShuffleCorr','TmShuffleDis','CandCorr','OutFR','CandDis','CandSeq','MidTime','params')

                            end
                            
                            tms1 = CandSeq(:,1)>params.Sleep_Times(1,1) & CandSeq(:,2)<params.Sleep_Times(1,2) & CandSeq(:,2)<params.Run_Times(1,1);
                            tms2 = CandSeq(:,1)>params.Sleep_Times(2,1) & CandSeq(:,1)<params.Sleep_Times(2,2);
                         
                            
                            %get sleep0 (presleep)
                            if icontrols==1
                                tms = tms1;
                            elseif icontrols==2
                                [minnumtms,o] = min([sum(tms1),sum(tms2)]);
                                tms = tms1;
                                if o==2
                                   toscr = find(tms1); toscr2 = toscr; 
                                   tms(toscr2(1:(sum(tms1)-minnumtms))) = false;
                                end
                            elseif icontrols ==3
                                [mintime,o2] = min([diff(params.Sleep_Times(1,:)),diff(params.Sleep_Times(2,:))]);
                                if o2==2
                                    tms = CandSeq(:,1)>params.Sleep_Times(1,1) & CandSeq(:,2)<params.Sleep_Times(1,1)+mintime & CandSeq(:,2)<params.Run_Times(1,1);
                                else
                                    tms = tms1;
                                end
                            end
                            tts0 = sum(tms);
                            ts0=sum(abs(CandCorr)>wc_cutoff & (CandDis/size(OutFR,2))<jd_cutoff & tms); 
                            ts0_shuff = sum((TmShuffleDis(tms,:)./size(OutFR,2))<jd_cutoff & abs(TmShuffleCorr(tms,:))>wc_cutoff,1);       
                            sts0=sum(abs(CandCorr(tms))>wc_cutoff & CandDis(tms)/size(OutFR,2)<jd_cutoff); 

                            %get any run1 replays
                            tms = CandSeq(:,1)>params.Run_Times(1,1) & CandSeq(:,1)<params.Run_Times(1,2);    
                            str1=sum(abs(CandCorr(tms))>wc_cutoff & CandDis(tms)/size(OutFR,2)<jd_cutoff ); 

                            %before and after one pass across the track
                            ss = find(MidTime>params.Run_Times(1,1) & MidTime<params.Run_Times(1,2));
                            tms = CandSeq(:,1)<MidTime(ss(1)) & CandSeq(:,1)>params.Run_Times(1,1) & CandSeq(:,1)<params.Run_Times(1,2);    
                            ttr11 = sum(tms);
                            tr11=sum(abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & tms); 
                            tr11_shuff = sum((TmShuffleDis(tms,:)./size(OutFR,2))<jd_cutoff & abs(TmShuffleCorr(tms,:))>wc_cutoff,1);                
                            tms = CandSeq(:,1)>MidTime(ss(1)) & CandSeq(:,1)>params.Run_Times(1,1) & CandSeq(:,1)<params.Run_Times(1,2);    
                            ttr12 = sum(tms);
                            tr12=sum(abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & tms); 
                            tr12_shuff = sum((TmShuffleDis(tms,:)./size(OutFR,2))<jd_cutoff & abs(TmShuffleCorr(tms,:))>wc_cutoff,1);

                            %get sleep1 (btwn run1 and run2)
                            if icontrols==1
                                tms = tms2;
                            elseif icontrols==2
                                tms = tms2;                                    
                                if o==1
                                   toscr = find(tms2); toscr2 = toscr; 
                                   tms(toscr2(1:(sum(tms2)-minnumtms))) = false;
                                end
                                clear minnumtms o
                            elseif icontrols==3
                                if o2==1
                                    tms = CandSeq(:,1)>params.Sleep_Times(2,1) & CandSeq(:,2)<params.Sleep_Times(2,1)+mintime;
                                else
                                    tms = tms2;
                                end
                                clear mintime o2
                            end
                            tts1 = sum(tms);
                            ts1=sum(abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & tms); 
                            ts1_shuff = sum((TmShuffleDis(tms,:)./size(OutFR,2))<jd_cutoff & abs(TmShuffleCorr(tms,:))>wc_cutoff,1);        
                            sts1=sum(abs(CandCorr(tms))>wc_cutoff & CandDis(tms)/size(OutFR,2)<jd_cutoff); 

                            %get run2 replays
                            tms = CandSeq(:,1)>params.Run_Times(2,1) & CandSeq(:,1)<params.Run_Times(2,2);    
                            ttr2 = sum(tms);
                            tr2=sum(abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & tms); % & CandStepSize(:,5,2)>coveragecutoff
                            tr2_shuff = sum((TmShuffleDis(tms,:)./size(OutFR,2))<jd_cutoff & abs(TmShuffleCorr(tms,:))>wc_cutoff,1);     
                            str2=sum(abs(CandCorr(tms))>wc_cutoff & CandDis(tms)/size(OutFR,2)<jd_cutoff ); 

                            %get sleep2 (post run2)
                            tms = CandSeq(:,1)>params.Sleep_Times(3,1) & CandSeq(:,1)<params.Sleep_Times(3,2);        
                            tts2 = sum(tms);
                            ts2=sum(abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & tms); 
                            ts2_shuff = sum((TmShuffleDis(tms,:)./size(OutFR,2))<jd_cutoff & abs(TmShuffleCorr(tms,:))>wc_cutoff,1);       
                            sts2=sum(abs(CandCorr(tms))>wc_cutoff & CandDis(tms)/size(OutFR,2)<jd_cutoff ); 

                           
                            allreal = cat(1,allreal,[ts0,tr11,tr12,ts1,tr2,ts2]);
                            alltotals = cat(1,alltotals,[tts0,ttr11,ttr12,tts1,ttr2,tts2]);
                            allsig = cat(1,allsig,[sts0,str1,sts1,str2,sts2]);
                            allshuff = cat(1,allshuff,cat(3,ts0_shuff,tr11_shuff,tr12_shuff,ts1_shuff,tr2_shuff,ts2_shuff));                                            
                           
                        end

                        rl = sum(allreal,1);
                        isz = sum(alltotals,1)==0;
                        sh = squeeze(sum(allshuff,1));
                        clear p allp p2
                        for ip = 1:size(allreal,2)
                            p(ip,1) = (1+sum(sh(:,ip)>=rl(ip)))./(1+size(sh,1));            
                        end                
                        p(isz) = NaN;

                        if icontrols ==1
                            dat(wc,jd,:,ifield) = p;
                            real_all(wc,jd,:,ifield) = squeeze(real_all(wc,jd,:,ifield))+rl';                    
                            sh_all(wc,jd,:,ifield,:) = squeeze(sh_all(wc,jd,:,ifield,:))+sh';
                        else
                            real_all_controls(wc,jd,:,ifield,icontrols-1) = squeeze(real_all(wc,jd,:,ifield))+rl';                    
                            sh_all_controls(wc,jd,:,ifield,:,icontrols-1) = squeeze(sh_all(wc,jd,:,ifield,:))+sh';
                        end
                        
                        if wcd(wc)>=.6 && jdd(jd)<=.4
                            quadreal(:,ifield,icontrols) = quadreal(:,ifield,icontrols)+rl';
                            quadshuff(:,:,ifield,icontrols) = quadshuff(:,:,ifield,icontrols)+sh';
                        end
                        disp(['added ' num2str(wc) ' ' num2str(jd) ' ' num2str(ifield) ' ' controllab{icontrols} ' tms ' num2str(tts0) ' ' num2str(tts1)])
                    end
                end
            end
            star = cat(3,star,[dat(4,4,1,1) dat(4,3,1,1) dat(3,3,1,1); dat(4,4,4,1) dat(4,3,4,1) dat(3,3,4,1)]);


            numshuff = size(allshuff,2);
            quadshuff = quadshuff(:,1:numshuff,:,:);             
        end
    end
    
    if plot_paper_figures
        epochlab = {'Sleep0';'Run1 before run';'Run1 after run';'Sleep1';'Run2';'Sleep2'};
        
        for icontrol = 1:3
            if icontrol==1
                dat_ALL = (1+sum(sh_all>=real_all,5))./(1+size(sh_all,5));
            else                
                dat_ALL = (1+sum(sh_all_controls(:,:,:,:,:,icontrol-1)>=real_all_controls(:,:,:,:,icontrol-1),5))./(1+size(sh_all,5));
            end
            
            N = 64;
            n = N/2;
            cm = NaN(N,3);

            cm(:,1) = [ones(n,1);linspace(1,0,N-n)';];
            cm(:,2) = [linspace(0,1,n)';linspace(1,0,N-n)']; 
            cm(:,3) = [linspace(0,1,n)';ones(N-n,1)]; 

            ifield = 1;
            figure; hold on
            set(gcf,'Position',[  66         518        2105         289])    
            for iepoch = 1:size(epochlab,1)
                subplot(1,size(epochlab,1),iepoch); hold on
                imagesc(log10(dat_ALL(:,:,iepoch,ifield))); axis xy;         

                set(gca,'ytick',1:3:length(wcd),'yticklabel',wcd(1:3:end))
                set(gca,'xtick',2:3:length(jdd),'xticklabel',jdd(2:3:end))
                xlabel('Abs(Max Jump Distance)')
                ylabel('Abs(Weighted Correlation)')
                set(gca,'clim',[log10(.05)*2 0])
                set(gcf,'colormap',cm)
                colorbar('Direction','reverse','Ticks',[log10(.005),log10(.05),log10(.5)],'TickLabels',[.005,.05,.5])
                x = [.5 4.5 4.5 .5];
                y = [.5 .5 4.5 4.5] ; 
                [xx,yy] = find(isnan(dat_ALL(:,:,iepoch,ifield)));        
                xxr = [xx-.5 xx+.5 xx+.5 xx-.5];
                yyr = [yy-.5 yy-.5 yy+.5 yy+.5];        
                patch(xxr',yyr','k','FaceColor','k','EdgeColor','k','LineWidth',1)
                patch(x,y,'w','FaceColor','none','EdgeColor','green','LineWidth',2)
                plot(3,3,'*g','MarkerSize',5)
                plot(4,3,'*g','MarkerSize',5)
                plot(4,4,'*g','MarkerSize',5)
                axis tight
                title([epochlab{iepoch}])
                set(gca,'FontSize',12)
            end
             helper_saveandclosefig([dirs.figdir 'Paper\SigAcrossEpochs_Run1in2parts_' fieldlab{ifield} '_' num2str(numshuff) '_TimeBinShuffle_' controllab{icontrol}])
         
        end         
         
    end
            
    if plot_paper_figures
        for icontrol = 1:3
            quadp = (sum(quadshuff(:,:,1,icontrol)>=quadreal(:,1,icontrol),2)+1)./(size(quadshuff,2)+1);
            figure; hold on
            for iq = 1:length(quadp)
               text(iq/10,iq/10,[epochlab{iq} ', p = ' num2str(quadp(iq))])
            end
            helper_saveandclosefig([dirs.figdir 'SigAcrossEpochs_Run1in2parts_' fieldlab{ifield} '_' num2str(numshuff) '_TimeBinShuffle_QuadPvalues_' controllab{icontrol}])
        end
    end 
    
    
    
end

%% 

%%%%%
%%%%% Figures of the significance tests of replays across paramters using
%%%%% replays with non-overlapping bins
%%%%%
SigTestOverLaps_NonOverLapping(basedir)

%%

%%%%%
%%%%% Get significance values off of figures if I already made them in 
%%%%% SigTestOverLaps_NonOverLapping.m
%%%%%

if plot_paper_figures
    Ps = NaN(6,3);
    uiopen([dirs.figdir 'Paper\SigAcrossEpochs_Run1in2parts_Run1 Fields_5000_TimeBinShuffle.fig'],1)

    f = gcf;
    Sleep0 = 10.^(f(1).Children(12).Children(5).CData);
    Ps(1,:) = round([Sleep0(4,4) Sleep0(3,4) Sleep0(3,3)],2,'significant');
    Sleep1 = 10.^(f(1).Children(6).Children(5).CData); 
    Ps(2,:) = round([Sleep1(4,4) Sleep1(3,4) Sleep1(3,3)],2,'significant');
    close all
    
    %and for across novel laps
    if ~isfile([dirs.homedir '\Figures\NonOverlappingLapByLap\SigTestOverShortLaps2_Tm.fig'])
        SigTestOverLaps_NonOverLapping(dirs.homedir)
    end
    uiopen([dirs.homedir '\Figures\NonOverlappingLapByLap\SigTestOverShortLaps2_Tm.fig'],1)

    f = gcf;
    for ilap = 0:3
        Lap = 10.^(f.Children(32-(ilap*2)).Children(5).CData); % Lap 0-3
        Ps(ilap+3,:) = round([Lap(4,4) Lap(3,4) Lap(3,3)],2,'significant');
    end
    close all

    epochlab = {'Sleep 0';'Sleep 1';'Lap 0';'Lap 1';'Lap 2';'Lap 3'};
    figure; hold on
    for iep = 1:length(epochlab)
        text(iep/10,iep/10,[epochlab{iep} ': ' num2str(Ps(iep,:))])
    end
     helper_saveandclosefig([dirs.figdir 'Paper\SigAcrossEpochs_Pvalues'])
end
