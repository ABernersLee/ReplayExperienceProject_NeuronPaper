function schematic(basedir)
%%%% creatiting example replays with different step sizes and such
%%%% for schematic


%%
mult = 10;
dat1 = [0 pi; pi 2*pi; 2*pi 3*pi; 3*pi 4*pi; 4*pi 5*pi; 5*pi 6*pi; 6*pi 7*pi; 7*pi 8*pi; 8*pi 9*pi; 9*pi 10*pi];
dat2 = [-2 -2; -2 -3; -3 -3; -3 -4; -4 -4; -4 -5; -5 -5; -5 -6; -6 -6; -6 -7];
figure; hold on; plot(0:.1:mult*pi,-cos([0:.1:mult*pi]),'k','LineWidth',3); 
line(dat1(1:2:end,:)',dat2(1:2:end,:)','color','k','LineWidth',3)
line(dat1(2:2:end,:)',dat2(2:2:end,:)','color','b','LineWidth',3)
line([0:2*pi:mult*pi; 0:2*pi:mult*pi],[[-7 2]'*ones(size([0:2*pi:mult*pi]))],'color',[.5 .5 .5],'LineStyle','--')
xlim([-1 mult*pi])
plot([-1 -1],[-2 -7],'k','LineWidth',3)
axis off
helper_saveandclosefig([basedir 'Figures\Final\schematic1'])


%%
mult = 6;
dat1 = [0 pi; pi 2*pi; 2*pi 3*pi; 3*pi 4*pi; 4*pi 5*pi; 5*pi 6*pi; 6*pi 7*pi; 7*pi 8*pi; 8*pi 9*pi; 9*pi 10*pi];
dat2 = [-2 -2; -2 -2-(5/3); -2-(5/3) -2-(5/3); -2-(5/3) -2-(10/3); -2-(10/3) -2-(10/3);-2-(10/3) -2-(15/3) ];
figure; hold on;
plot(0:.1:mult*pi,-cos([0:.1:mult*pi]),'k','LineWidth',3); 
line(dat1(1:2:mult,:)',dat2(1:2:mult,:)','color','k','LineWidth',3)
line(dat1(2:2:mult,:)',dat2(2:2:mult,:)','color','b','LineWidth',3)
line([0:2*pi:mult*pi; 0:2*pi:mult*pi],[[-7 2]'*ones(size([0:2*pi:mult*pi]))],'color',[.5 .5 .5],'LineStyle','--')
xlim([-1 10*pi])
plot([-1 -1],[-2 -7],'k','LineWidth',3)
axis off
helper_saveandclosefig([basedir 'Figures\Final\schematic_shortreplay'])
set(gcf,'Position',[ 2042         329         864         420])
helper_saveandclosefig([basedir 'Figures\Final\schematic_shortreplay_stretch'])
%%

mult = 6;

dat1 = [0 pi/2; pi/2 2*pi; 2*pi 2.5*pi; 2.5*pi 4*pi; 4*pi 4.5*pi; 4.5*pi 6*pi; 6*pi 6.5*pi; 6.5*pi 8*pi; 8*pi 8.5*pi; 8.5*pi 10*pi];
% dat2 = [-2 -2; -2 -3; -3 -3; -3 -4; -4 -4; -4 -5; -5 -5; -5 -6; -6 -6; -6 -7];
dat2 = [-2 -2; -2 -2-(5/3); -2-(5/3) -2-(5/3); -2-(5/3) -2-(10/3); -2-(10/3) -2-(10/3);-2-(10/3) -2-(15/3) ];
figure; hold on; plot(0:.1:mult*pi,-cos([0:.1:mult*pi]),'k','LineWidth',3); 
line(dat1(1:2:mult,:)',dat2(1:2:mult,:)','color','k','LineWidth',3)
line(dat1(2:2:mult,:)',dat2(2:2:mult,:)','color','b','LineWidth',3)
line([0:2*pi:mult*pi; 0:2*pi:mult*pi],[[-7 2]'*ones(size([0:2*pi:mult*pi]))],'color',[.5 .5 .5],'LineStyle','--')
xlim([-1 10*pi])
plot([-1 -1],[-2 -7],'k','LineWidth',3)
axis off
helper_saveandclosefig([basedir 'Figures\Final\schematic_longjump'])
%%

% mult = 10;
dat1 = [0 pi*1.5; pi*1.5 2*pi; 2*pi 3.5*pi; 3.5*pi 4*pi; 4*pi 5.5*pi; 5.5*pi 6*pi; 6*pi 7.5*pi; 7.5*pi 8*pi; 8*pi 9.5*pi; 9.5*pi 10*pi];
dat2 = [-2 -2; -2 -2-(5/3); -2-(5/3) -2-(5/3); -2-(5/3) -2-(10/3); -2-(10/3) -2-(10/3);-2-(10/3) -2-(15/3) ];
figure; hold on; plot(0:.1:mult*pi,-cos([0:.1:mult*pi]),'k','LineWidth',3); 
line(dat1(1:2:mult,:)',dat2(1:2:mult,:)','color','k','LineWidth',3)
line(dat1(2:2:mult,:)',dat2(2:2:mult,:)','color','b','LineWidth',3)
line([0:2*pi:mult*pi; 0:2*pi:mult*pi],[[-7 2]'*ones(size([0:2*pi:mult*pi]))],'color',[.5 .5 .5],'LineStyle','--')
xlim([-1 10*pi])
plot([-1 -1],[-2 -7],'k','LineWidth',3)
axis off
% set(gcf,'Position',[ 1149         275         788         420])
helper_saveandclosefig([basedir 'Figures\Final\schematic_longhover'])
%%
mult = 10;
dat1 = [0 pi; pi 2*pi; 2*pi 3*pi; 3*pi 4*pi; 4*pi 5*pi; 5*pi 6*pi; 6*pi 7*pi; 7*pi 8*pi; 8*pi 9*pi; 9*pi 10*pi];
dat2 = [-2 -2; -2 -3; -3 -3; -3 -4; -4 -4; -4 -5; -5 -5; -5 -6; -6 -6; -6 -7];
figure; hold on;
plot(0:.1:mult*pi,-cos([0:.1:mult*pi]),'k','LineWidth',3); 
xlim([-1 mult*pi])
ylim([-7 2])
axis off
helper_saveandclosefig([basedir 'Figures\Final\schematic_raw'])
%%

mult = 10;
dat1 = [0 pi*1.5; pi*1.5 2*pi; 2*pi 3.5*pi; 3.5*pi 4*pi; 4*pi 5.5*pi; 5.5*pi 6*pi; 6*pi 7.5*pi; 7.5*pi 8*pi; 8*pi 9.5*pi; 9.5*pi 10*pi];
dat2 = [-2 -2; -2 -3; -3 -3; -3 -4; -4 -4; -4 -5; -5 -5; -5 -6; -6 -6; -6 -7];
figure; hold on; plot(0:.1:mult*pi,-cos([0:.1:mult*pi]),'k','LineWidth',3); 
line(dat1(1:2:end,:)',dat2(1:2:end,:)','color','k','LineWidth',3)
line(dat1(2:2:end,:)',dat2(2:2:end,:)','color','b','LineWidth',3)
line([0:2*pi:mult*pi; 0:2*pi:mult*pi],[[-7 2]'*ones(size([0:2*pi:mult*pi]))],'color',[.5 .5 .5],'LineStyle','--')
xlim([-1 mult*pi])
plot([-1 -1],[-2 -7],'k','LineWidth',3)
axis off
helper_saveandclosefig([basedir 'Figures\Final\schematic5'])