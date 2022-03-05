function RunReplay_ABL_RemoteSleep_Figures(basedir,plot_paper_figures,plot_other_figures)

%%%%%
%%%%% This looks at the replays after many laps within the data where there
%%%%% is a rest after one novel track before the next one
%%%%%

%%

%%%%%
%%%%% directories
%%%%%

load([basedir 'dirs_remotesleep'],'dirs')

if ~isfolder([dirs.figdir '\Paper\'])
    mkdir([dirs.figdir '\Paper\'])
end
%% 

%%%%%
%%%%% look at fields and behavior including example figure for paper
%%%%%


cd(dirs.spikedatadir)
d2 = dir('*.mat');

for id = 1:size(d2,1)
    daydir = d2(id).name;
    clear InFR OutFR
    load(daydir,'InFR','OutFR','pos','linposcat','params','rundat','MidTime')   
    
    for irun = 1:2
        if plot_other_figures
        [m,mx] = max(InFR(:,:,irun),[],2);
        [~,o] = sort(mx);
        figure; hold on
        subplot(1,2,1); hold on
        imagesc(InFR(o,:,irun)./m(o)); axis ij; colormap hot; axis tight
        title('In'); xlabel('Position Bin'); ylabel('Neuron')
        [m,mx] = max(OutFR(:,:,irun),[],2);
        [~,o] = sort(mx);        
        subplot(1,2,2); hold on
        imagesc(OutFR(o,:,irun)./m(o)); axis ij; colormap hot; axis tight
        title('Out'); xlabel('Position Bin'); ylabel('Neuron')
        suptitle([daydir(1:end-4) ' Run ' num2str(irun) ' fields'])
        helper_saveandclosefig([dirs.figdir '\fields_' daydir(1:end-4) '_run' num2str(irun) 'fields'])
        end
        
        if plot_paper_figures
        figure; hold on; 
        plot((pos(rundat(:,irun),1) - min(pos(rundat(:,irun),1)))/60,linposcat(rundat(:,irun))-min(linposcat(rundat(:,irun))),'k'); 
    %     hold on; plot(MidTime(:),min(linposcat)+range(linposcat)/2,'r*')
        set(gcf,'Position',[ 1065         285        1921         366])
        set(gcf,'renderer','Painters')
        helper_saveandclosefig([dirs.figdir '\Paper\laps_new_' daydir '_run' num2str(irun)])    
        end
    end
    
end

%%  

%%%%% looking at duration and slope across presleep, first experience, 
%%%%% sleep, second experience, post-sleep,  using run1 or 2 fields

if plot_paper_figures

    shufftoo = 0;
    shufftoolab = {'';'_Shufftoo'};
    mlaptimeskip = 1; 
    mlaptimeskip2 = 1;

    cd(dirs.spikedatadir)
    d2 = dir('*.mat');

    numcells = NaN(size(d2,1),1); lapsall = [];
    for idd = 1:size(d2,1)
         load([d2(idd).name],'hp_cells','hpinterneurons','MidTime','params')
            if exist('hp_cells','var')
                numcells(idd,1) = sum(~ismember(hp_cells,hpinterneurons));
                clear hp_cells hpinterneurons
            end

            MidTime1 = MidTime(MidTime>params.Run_Times(1,1) & MidTime<params.Run_Times(1,2));
            MidTime2 = MidTime(MidTime>params.Run_Times(2,1) & MidTime<params.Run_Times(2,2));
        lapsall = cat(1,lapsall,[length(MidTime1)-1 length(MidTime2)-1]);
    end

    cd(dirs.spikedatadir)
    jd_cutoff = .4; wc_cutoff = .6; coveragecutoff = 0; 
    sleep2replays = []; 
    run2replays = []; run1replays = [];
    numlaps_epoch = NaN(size(d2,1),3,2);

    for idd = 1:size(d2,1)
        run1replays2 = []; run2replays2 = []; sleep2replays2 = [];
        for ifield = 1:2
            load(d2(idd).name,'CandCorr','OutFR','CandDis','CandRS','CandSeq','MidTime','CandStepSize','params','pos')
            %didn't save replayparams out but as seen in
            %decode_spikedensity_events_remote_sleep.m, the Bin_Size is 2.5, so
            %making variable here (get to in RunReplay_ABL_RemoteSleep.m)
            replayparams.Bin_Size = 2.5;
            CandCorr = CandCorr(:,ifield);
            CandStepSize = CandStepSize(:,:,ifield);
            CandDis = CandDis(:,ifield);
            CandRS = CandRS(:,:,ifield);
            OutFR = OutFR(:,:,ifield);

            MidTime1 = MidTime(MidTime>params.Run_Times(1,1) & MidTime<params.Run_Times(1,2));
            MidTime2 = MidTime(MidTime>params.Run_Times(2,1) & MidTime<params.Run_Times(2,2));
            mlaptime = median([diff(MidTime1); diff(MidTime2)]);


                     %1/26/2022 change* need to change unites of slope!! 5*old#     *replayparams.Bin_Size/100*1000/5

            %sleep 2
            MidSleep = params.Sleep_Times(1,1):mlaptime:params.Sleep_Times(1,2);
            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidSleep(1:mlaptimeskip:end)]);        
            numlaps_epoch(idd,2,ifield) =max(I);       
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,1)>=coveragecutoff & CandSeq(:,1)>params.Sleep_Times(1,1) & CandSeq(:,1)<params.Sleep_Times(1,2); % & sigst1<.05 & sigst2<.05;
            sleep2replays1 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 CandStepSize(t,3,1)*replayparams.Bin_Size CandStepSize(t,4,1)*replayparams.Bin_Size CandSeq(t,1)-params.Sleep_Times(1,1); NaN(sum(~t),6)];


            %get any run1 replays
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,1)>=coveragecutoff & CandSeq(:,1)>params.Run_Times(1,1) & CandSeq(:,1)<params.Run_Times(1,2); % & sigst1<.05 & sigst2<.05;
            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime1(1:mlaptimeskip2:end);params.Run_Times(1,2)]);
            I(I>length(MidTime)) = length(MidTime);
            numlaps_epoch(idd,1,ifield) =max(I); 
            run1replays1 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 CandStepSize(t,3,1)*replayparams.Bin_Size CandStepSize(t,4,1)*replayparams.Bin_Size CandSeq(t,1)-params.Run_Times(1,1); NaN(sum(~t),6)];

            %get run2 replays
            t=abs(CandCorr)>wc_cutoff & CandDis/size(OutFR,2)<jd_cutoff & CandStepSize(:,5,1)>=coveragecutoff & CandSeq(:,1)>params.Run_Times(2,1) & CandSeq(:,1)<params.Run_Times(2,2); % & sigst1<.05 & sigst2<.05;
            [~,I]=histc(CandSeq(:,1)/2+CandSeq(:,2)/2,[MidTime2(1:mlaptimeskip2:end);params.Run_Times(2,2)]);
            I(I>length(MidTime)) = length(MidTime);
            numlaps_epoch(idd,3,ifield) =max(I);         
            run2replays1 = [I(t) CandSeq(t,2)-CandSeq(t,1) abs(CandRS(t,2))*replayparams.Bin_Size/100*1000/5 CandStepSize(t,3,1)*replayparams.Bin_Size CandStepSize(t,4,1)*replayparams.Bin_Size CandSeq(t,1)-params.Run_Times(2,1); NaN(sum(~t),6)];

            run1replays2 = cat(3,run1replays2,run1replays1);
            run2replays2 = cat(3,run2replays2,run2replays1);
            sleep2replays2 = cat(3,sleep2replays2,sleep2replays1);
        end

        run1replays = cat(1,run1replays,run1replays2);        
        run2replays = cat(1,run2replays,run2replays2);
        sleep2replays = cat(1,sleep2replays,sleep2replays2);
    end

    cutofflaps = squeeze(median(numlaps_epoch,1));
    cutofflaps([1 3],:) = round(14./mlaptimeskip2);

    run1replays(run1replays(:,1,1)>cutofflaps(1,1),:,1) = NaN;
    sleep2replays(sleep2replays(:,1,1)>cutofflaps(2,1),:,1) = NaN;
    run2replays(run2replays(:,1,1)>cutofflaps(3,1),:,1) = NaN;
    run1replays(run1replays(:,1,2)>cutofflaps(1,2),:,2) = NaN;
    sleep2replays(sleep2replays(:,1,2)>cutofflaps(2,2),:,2) = NaN;
    run2replays(run2replays(:,1,2)>cutofflaps(3,2),:,2) = NaN;
    run1replays(sum(sum(isnan(run1replays),2),3)==12,:,:) = [];
    run2replays(sum(sum(isnan(run2replays),2),3)==12,:,:) = [];
    sleep2replays(sum(sum(isnan(sleep2replays),2),3)==12,:,:) = [];

    % Duration and Slope
    tlab = {'Duration';'Slope';'rangeDtime'; 'mean(abs(moves))'};
    for itype = 2:3 %:5
        xt = []; xtl = [];
        figure; hold on;
        numnow = 1;
        xlabpos = NaN(3,1);
        for ifield = 1 %:2

                %run1 replays
                A = NaN(max(run1replays(:,1,ifield)),2);
                for ilap = 1:max(run1replays(:,1,ifield))
                    A(ilap,1) = nanmean(run1replays(run1replays(:,1,ifield)==ilap,itype,ifield));
                    A(ilap,2) = nanstd(run1replays(run1replays(:,1,ifield)==ilap,itype,ifield))./sqrt(sum(run1replays(:,1,ifield)==ilap));
                end
                A(isnan(A(:,1)),:) = [];

                if ifield == 1
                    errorbar(1:size(A,1),A(:,1),A(:,2),'k','LineWidth',3)


                    xtn = get(gca,'xtick');
                    xt = [xt numnow+1:median(diff(xtn)):size(A,1)+numnow];
                    xtl = [xtl [2:median(diff(xtn)):size(A,1)+numnow]*mlaptimeskip2];

                    numnow = size(A,1)+8;
                    xlabpos(1) = round(size(A,1)./2);
                else
    %                 errorbar(1:size(A,1),A(:,1),A(:,2),'color',[.5 .5 .5],'LineWidth',3)
                end        



                %sleep2
                A = NaN(max(sleep2replays(:,1,ifield)),2);
                for ilap = 1:max(sleep2replays(:,1,ifield))
                    A(ilap,1) = nanmean(sleep2replays(sleep2replays(:,1,ifield)==ilap,itype,ifield));
                    A(ilap,2) = nanstd(sleep2replays(sleep2replays(:,1,ifield)==ilap,itype,ifield))./sqrt(sum(sleep2replays(:,1,ifield)==ilap));
                end
                A(isnan(A(:,1)),:) = [];


                %toplot
    %             A0 = NaN(max(sleep2replays(:,1,ifield)),2);
    %             for ilap = 1:2:max(sleep2replays(:,1,ifield))
    %                 A0(ilap,1) = nanmean(sleep2replays(sleep2replays(:,1,ifield)==ilap | sleep2replays(:,1,ifield)==ilap,itype,ifield));
    %                 A0(ilap,2) = nanstd(sleep2replays(sleep2replays(:,1,ifield)==ilap | sleep2replays(:,1,ifield)==ilap,itype,ifield))./sqrt(sum(sleep2replays(:,1,ifield)==ilap | sleep2replays(:,1,ifield)==ilap));
    %             end
    %             A0(isnan(A0(:,1)),:) = [];


                if ifield ==1
                    Filter=fspecial('gaussian',[10 1],4);
                    pd = 10/mlaptimeskip;
                    A0 = A;
                    A2 = filter2(Filter,[A0(1:pd,1);A0(:,1);A0(end-(pd-1):end,1)]);
                    A3 = filter2(Filter,[A0(1:pd,2);A0(:,2);A0(end-(pd-1):end,2)]);
                    A2 = A2(pd+1:end-pd,1);
                    A3 = A3(pd+1:end-pd,1);
    %                 A2 = A(:,1); A3 = A(:,2);
                    fwd = [A2-A3]'; rev = [A2+A3]';
                    p2 = patch([numnow:size(A0,1)+numnow-1 size(A0,1)+numnow-1:-1:numnow],[fwd rev(end:-1:1)],'black');
                    p2.FaceAlpha=.1; p2.EdgeAlpha=0;
                    plot(numnow:size(A0,1)+numnow-1,A2,'k','LineWidth',3)                
                    xlabpos(2) = numnow+round((size(A0,1))./2);


                    xtn = get(gca,'xtick');
                    xt = [xt numnow:median(diff(xtn)):size(A0,1)+numnow-1];
                    xtl = [xtl [1:median(diff(xtn)):size(A0,1)+numnow-1]*mlaptimeskip2];

                    numnow = size(A0,1)+numnow;
                else
    %                 errorbar(numnow:size(A,1)+numnow-1,A(:,1),A(:,2),'color',[.5 .5 .5],'LineWidth',3)
                end         


        end

    %       

        ylabel(tlab{itype-1})
        set(gca,'FontSize',18)


         %is the first run changing across runs differently than the sleep and run2 with
         %the run1 fields?     
        p1 = anovan([run1replays(:,itype,1);sleep2replays(:,itype,1); run2replays(:,itype,1)],[[ones(size(run1replays,1),1)*1; ones(size(sleep2replays,1),1)*2; ones(size(run2replays,1),1)*3]...
            [ run1replays(:,1,1); sleep2replays(:,1,1); run2replays(:,1,1)]],'model','interaction','continuous',2,'display','off');
         %specifically to the sleep, is the first run chaning across laps more
         %than during subsequent sleep
        [p2,tbl] = anovan([run1replays(:,itype,1);sleep2replays(:,itype,1)],[[ones(size(run1replays,1),1)*1;ones(size(sleep2replays,1),1)*2;]...
            ,[ run1replays(:,1,1); sleep2replays(:,1,1)]],'model','interaction','continuous',2,'display','off');
         %during run2 are the local replays changing more than the remote
         a = [run2replays(:,itype,1);run2replays(:,itype,2)];
         b = [[ones(size(run2replays,1),1)*1;ones(size(run2replays,1),1)*2;]...
            ,[ run2replays(:,1,1); run2replays(:,1,2)]];
        p3 = anovan(a(sum(isnan(b),2)==0),b(sum(isnan(b),2)==0,:),'model','interaction','continuous',2,'display','off');

         % are the replays a different length between the sleep and subsequent remote experience for run1 fields
        p0 = ranksum(run2replays(:,itype,1),sleep2replays(:,itype,1));

        if itype==2
            [r1,rp1] = corr(run1replays(:,1,1),run1replays(:,itype,1),'rows','complete','tail','right');
            [r2,rp2] = corr(sleep2replays(:,1,1),sleep2replays(:,itype,1),'rows','complete','tail','right');
        elseif itype==3        
            [r1,rp1] = corr(run1replays(:,1,1),run1replays(:,itype,1),'rows','complete','tail','left');
            [r2,rp2] = corr(sleep2replays(:,1,1),sleep2replays(:,itype,1),'rows','complete','tail','left');
        end
        set(gca,'xlim',[0 numnow])
        set(gcf,'Position',[  1934         298         983         326])
        set(gcf, 'Renderer', 'painters');

        title(['n1 = ' num2str(sum(~isnan(run1replays(:,1,1)))) ' r='  num2str(round(r1,2,'significant')) ' p='  num2str(round(rp1,2,'significant')) ...
            ', n2=' num2str(sum(~isnan(sleep2replays(:,1,1)))) ' r='  num2str(round(r2,2,'significant')) ' p='  num2str(round(rp2,2,'significant'))])

        set(gca,'xtick',xt,'xticklabel',xtl)
        helper_saveandclosefig([dirs.figdir '\Paper\AcrossEpochs_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1}  '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1} '_cutoffsleep'])


        figure; hold on;
        plab = {'Group';'Lap';'Interaction'};
        for ip = 1:3
            text(0,ip/4,[plab{ip} ': F = ' num2str(tbl{ip+1,6}) ' df = ' num2str(tbl{ip+1,3}) ',' num2str(sum(~isnan([run1replays(:,itype,1);sleep2replays(:,itype,1)]))-tbl{ip+1,3}) ' , p = ' num2str(p2(ip))])
        end
        set(gcf,'Position',[1919         133         861         403])
        axis off
        helper_saveandclosefig([dirs.figdir '\Paper\AcrossEpochs_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_ANOVA_' tlab{itype-1}  '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1} '_cutoffsleep'])


        figure; hold on;
        text(.1,.8,tlab{itype-1})
        text(.1,.7,['Is the first run changing across runs differently than the sleep and run2 with the run1 fields? Ip='  num2str(round(p1(3),3,'significant'))])
        text(.1,.6,['Specifically to the sleep, is the first run chaning across laps more than during subsequent sleep?, Ip='  num2str(round(p2(3),3,'significant'))])
        text(.1,.5,['During run2, are the local replays changing more than the remote? Ip='  num2str(round(p3(3),3,'significant'))])
        text(.1,.4,['Are the replays a different length between the sleep and subsequent remote experience for run1 fields? p='  num2str(round(p0,3,'significant'))])
        set(gcf,'Position',[1927         229        1024         405])
        helper_saveandclosefig([dirs.figdir '\Paper\AcrossEpochs_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1}  '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1} '_stats' '_cutoffsleep'])

        totest1 = run1replays(~isnan(run1replays(:,1,1)),:,1);

        totest2 = sleep2replays(~isnan(sleep2replays(:,1,1)),:,1);
        totest = [totest1; totest2];

         %Chow test
                %         Chow tests rely on:
                % 
                % Independent, Gaussian-distributed innovations
                % 
                % Constancy of the innovations variance within subsamples
                % 
                % Constancy of the innovations across any structural breaks

            %BY LAP
    %         Mdl = fitlm(totest(:,1),totest(:,itype));
    %         res = Mdl.Residuals.Raw;
    %         [hARCH,pValueARCH] = archtest(res); %should be false
    %         [hKS,pValueKS] = kstest(res/std(res)); % should be false, seems to be not a normal distribution so do a permutation test too

            [~,p1,teststat] = chowtest(totest(:,1),totest(:,itype),size(totest1,1)+1,'Coeffs',[0 1]); %,'Display','summary');
            numperm = 1000; permstat = NaN(numperm,1);
            for ip = 1:numperm            
                [~,~,permstat(ip)] = chowtest(totest(:,1),totest(randperm(size(totest,1)),itype),size(totest1,1)+1,'Coeffs',[0 1]);
            end
            p = (1+(sum(permstat>=teststat)))./(1+numperm);
            figure; hold on; histogram(permstat,'FaceColor','none','LineWidth',2); 
            yl = get(gca,'ylim'); plot([teststat teststat],yl,'-r','LineWidth',3)
            xl = get(gca,'xlim');
            text(xl(2)*.5,yl(2)*.8,['F = ' num2str(teststat) ', N = ' num2str(size(totest,1)) ', bp = ' num2str(size(totest1,1)+1)])
            ylabel('Permutations')
            xlabel('Chow Test''s F')
            set(gca,'FontSize',18)
            title({['Permutation test on chow test statistic'],['Lap and ' tlab{itype-1} ', OGp = ' num2str(round(p1,2,'significant')) ', PERMp = ' num2str(round(p,2,'significant'))]})
            set(gcf, 'Renderer', 'painters');
            set(gcf,'Position',[2156        -676         666         544])
            helper_saveandclosefig([dirs.figdir '\Paper\ChowTest_Lap_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1} '_cutoffsleep'])
            disp(p)
    %         [h,p,teststat] = chowtest(log(totest(:,1)),log(totest(:,itype)),size(totest1,1)+1,'Coeffs',[0 1]); %,'Display','summary');

            %BY TIME
    %         Mdl = fitlm(totest(:,end),totest(:,itype));
    %         res = Mdl.Residuals.Raw;
    %         [hARCH,pValueARCH] = archtest(res); %should be false
    %         [hKS,pValueKS] = kstest(res/std(res)); % should be false, seems to be not a normal distribution so do a permutation test too

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
            helper_saveandclosefig([dirs.figdir '\Paper\ChowTest_Time_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1} '_cutoffsleep'])
            disp(p)

            %change between sleep and start of run 2
            figure; hold on
            text(.1,.9,tlab{itype-1})
            p = ranksum(totest1(totest1(:,1)>(10/mlaptimeskip2),itype),totest2(totest2(:,1)==1,itype));

            text(.1,.8,['Ranksum between Run and Lap 1 of Sleep, p = ' num2str(p)])
            text(.1,.7,['Run 1 N = ' num2str(sum(~isnan(totest1(:,itype)))) ', Sleep N = ' num2str(sum(~isnan(totest2(totest2(:,1)==1,itype))))])
            p = ranksum(totest1(totest1(:,1)>(10/mlaptimeskip2),itype),totest2(totest2(:,1)<(5/mlaptimeskip),itype));

            text(.1,.6,['Ranksum between Run laps > 10 and Laps 1-4 of Sleep, p = ' num2str(p)])
            text(.1,.5,['Run 1 N = ' num2str(sum(~isnan(totest1(:,itype)))) ', Sleep N = ' num2str(sum(~isnan(totest2(totest2(:,1)<5,itype))))])

            set(gcf,'Position',[ 2109         -29         785         326])
            helper_saveandclosefig([dirs.figdir '\Paper\Change_SleepRun_wc' num2str(wc_cutoff) '_jd' num2str(jd_cutoff) '_covcutoff' num2str(coveragecutoff) '_laps_forppt_' tlab{itype-1} '_mlaptimeskip' num2str(mlaptimeskip)  '_mlaptimeskip2_' num2str(mlaptimeskip2) shufftoolab{shufftoo+1} '_cutoffsleep'])            
    end
end
