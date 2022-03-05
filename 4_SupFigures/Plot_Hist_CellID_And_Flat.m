function Plot_Hist_CellID_And_Flat(basedir)

%%%%%
%%%%%
%%%%% This takes the data and figures generated in
%%%%% Hovers_Moving_FlatPosterior.m and Hovers_Moving_Shuffle and
%%%%% reorganizes them into new figures and stats for the supplement.
%%%%%
%%%%%

uiopen([basedir 'Figures\HoversMovingShuffle\corr_60cell_all_wholefields_wc0_jd1_coverage0_CellIDShuffle.fig'],1)
a1 = gca;
a1data1 = a1.Children(2).Data;
a1data = reshape(a1data1,[100 size(a1data1,1)/100])';
uiopen([basedir 'Figures\HoversMovingShuffle\corr_60cell_all_wholefields_wc0_jd1_coverage0_Flat_FamandNov.fig'],1)
a2 = gca;
a2data2 = a2.Children(2).Data;
a2data = reshape(a2data2,[100 size(a2data2,1)/100])';

figure; hold on; 
h1 = histogram(mean(a1data,2),-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor',[.5 .5 .5],'Normalization','probability'); 
h2 = histogram(mean(a2data,2),-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor','k','Normalization','probability'); 
yl = get(gca,'ylim');
plot([median(mean(a1data,2)) median(mean(a1data,2))],yl,'--','Color',[.5 .5 .5],'LineWidth',2)
text(median(mean(a1data,2))+.01,yl(2)*.98,num2str(round(median(mean(a1data,2)),2)),'Color',[.5 .5 .5],'FontSize',12)
plot([median(mean(a2data,2)) median(mean(a2data,2))],yl,'--','Color','k','LineWidth',2)
text(median(mean(a2data,2))+.01,yl(2)*.98,num2str(round(median(mean(a2data,2)),2)),'Color','k','FontSize',12)
legend([h1 h2],'Cell ID shuffle','Flat field distribution','Location','northwest')
xlabel('Correlation')
ylabel('Proportion of sessions')
set(gca,'FontSize',18)
set(gcf,'renderer','Painters')
helper_saveandclosefig([basedir 'Figures\HoversMovingShuffle\Hist_CellID_And_Flat_sessions'])

[p,~,stats] = signrank(mean(a1data,2),mean(a2data,2));
figure; hold on; 
text(.1,.9,['Cell ID Shuffle, N = ' num2str(length(mean(a1data,2))) ' Mean = ' num2str(mean(mean(a1data,2))) ', SEM = ' num2str(std(mean(a1data,2))./sqrt(sum(~isnan(mean(a1data,2)))))])
text(.1,.8,['Flat field distribution, N = ' num2str(length(mean(a2data,2))) ' Mean = ' num2str(mean(mean(a2data,2))) ', SEM = ' num2str(std(mean(a2data,2))./sqrt(sum(~isnan(mean(a2data,2)))))])
text(.1,.7,['signedrank = ' num2str(stats.signedrank) ', Z = ' num2str(stats.zval) ', P = ' num2str(p)])

helper_savefig([basedir 'Figures\Final\Hist_CellID_And_Flat_stats'])
helper_saveandclosefig([basedir 'Figures\HoversMovingShuffle\Hist_CellID_And_Flat_sessions_stats'])

%%
uiopen([basedir 'Figures\HoversMovingShuffle\corr_60cell_mean_wholefields_wc0_jd1_coverage0_CellIDShuffle.fig'],1)
a1 = gca;
a1data = a1.Children(2).Data;
uiopen([basedir 'Figures\HoversMovingShuffle\corr_60cell_mean_wholefields_wc0_jd1_coverage0_Flat_FamandNov.fig'],1)
a2 = gca;
a2data = a2.Children(2).Data;

figure; hold on; 
h1 = histogram(a1data,-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor',[.5 .5 .5],'Normalization','probability'); 
h2 = histogram(a2data,-.5:.05:1,'FaceColor','none','LineWidth',2,'EdgeColor','k','Normalization','probability'); 
yl = get(gca,'ylim');
plot([median(a1data) median(a1data)],yl,'--','Color',[.5 .5 .5],'LineWidth',2)
text(median(a1data)+.01,yl(2)*.98,num2str(round(median(a1data),2)),'Color',[.5 .5 .5],'FontSize',12)
plot([median(a2data) median(a2data)],yl,'--','Color','k','LineWidth',2)
text(median(a2data)+.01,yl(2)*.98,num2str(round(median(a2data),2)),'Color','k','FontSize',12)
legend([h1 h2],'Cell ID shuffle','Flat field distribution','Location','northwest')
xlabel('Correlation')
ylabel('Proportion of shuffles')
set(gca,'FontSize',18)
set(gcf,'renderer','Painters')
helper_savefig([basedir 'Figures\Final\Hist_CellID_And_Flat_mean'])
helper_saveandclosefig([basedir 'Figures\HoversMovingShuffle\Hist_CellID_And_Flat_mean'])

[p,~,stats] = ranksum(a1data,a2data);

figure; hold on; 
text(.1,.9,['Cell ID Shuffle, N = ' num2str(length(a1data)) ' Mean = ' num2str(mean(a1data)) ', SEM = ' num2str(std(a1data)./sqrt(sum(~isnan(a1data))))])
text(.1,.8,['Flat field distribution, N = ' num2str(length(a2data)) ' Mean = ' num2str(mean(a2data)) ', SEM = ' num2str(std(a2data)./sqrt(sum(~isnan(a2data))))])
text(.1,.7,['rank-sum, Z = ' num2str(stats.zval) ', P = ' num2str(p)])
helper_saveandclosefig([basedir 'Figures\HoversMovingShuffle\Hist_CellID_And_Flat_stats_mean'])

close all
