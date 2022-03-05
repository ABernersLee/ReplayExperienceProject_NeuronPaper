function [armpos,linposcat,linposnorm,linposcatnan,dirdat,cm_conv,rundat] = GUI_make_linear_position_v2(pos,params)

%%%%%
%%%%%
%%%%% GUI to make position and velocity data and check that there aren't
%%%%% any problems with the position tracking/data. Version used only for
%%%%% RunReplay_ABL_RemoteSleep.m
%%%%%
%%%%%

st = (repmat(pos(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(pos(:,1))]);
nd = (repmat(pos(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(pos(:,1))]);
rundat = (st&nd)';
a1 = []; p1 = []; p21 = [];
% Filter=fspecial('gaussian',[1 20],4);
for irun = 1:size(rundat,2)
    pos1 = pos(rundat(:,irun),:);
    armposdat = params.arms{irun};
    %get the dimention with the largest variablity
%     [~,I] = max(range(pos1(:,2:3)));    
%     p2 = pos1(:,I+1)-min(pos1(:,I+1))+.001;
%     p = pos1(:,1:3); 
    
    
    if ~isfield(params,'armslength')
        params.armslength = 160;
    end
    
    x = pos1(:,2); y = pos1(:,3);     
    len=sqrt((armposdat(2,1)-armposdat(1,1))^2+(armposdat(2,2)-armposdat(1,2))^2);
    xlen=abs((armposdat(2,1)-armposdat(1,1)))/len;
    ylen=abs((armposdat(2,2)-armposdat(1,2)))/len;
    pp=[x y]*[xlen ylen]';
     ss=[armposdat(1,1) armposdat(1,2)]*[xlen ylen]';
     ee=[armposdat(2,1) armposdat(2,2)]*[xlen ylen]';
     frac=min(1,max(0,((pp-ss)./(ee-ss))));
     qx=armposdat(1)+(frac*(armposdat(2,1)-armposdat(1,1)));
     qy=armposdat(2)+(frac*(armposdat(2,2)-armposdat(1,2)));
     distance=sqrt(((x-qx).^2)+((y-qy).^2));
     projpoint=frac*len;
    

    %normalize arm to arms length
    dat = projpoint; 
%     dat = distance;
    p2=((dat-min(dat))./(max(dat)-min(dat)))*params.armslength;
    cm_conv = params.armslength./range(projpoint);
    
    p2(p2==0) = p2(p2==0)+.0001;   
    
    p = pos1(:,1:3);
    a = ones(size(p2,1),1);
    p1 = [p1;p];
    a1 = [a1;a];
%     p2 = filter2(Filter,p2);
    p21 = [p21;p2];
    
    
end
p2 = p21; a=a1;p=p1;
po_new = p2;

% if 0 
%plots
figure; hold on
set(gcf,'Position',[417 179 1366 668])
subplot(1,2,1); hold on;
for i = 1:4
    plot(p(a==i,2),p(a==i,3),'.')
end
subplot(1,2,2)
plot(p2,'.')  


%checks if this is okay
prompt = 'Is this okay? y/n';
title(prompt)
yn = input([prompt ':  '],'s');
if strcmp(yn,'n')
    return
else
    close gcf
end
% end

armpos = a;
linposcat = nansum(p2,2);
linposnorm = nansum(po_new,2);
linposcatnan = [pos(:,1) p2];

close all

    

% skip = floor((length(linposnorm)/(max(linposnorm(:,1))-min(linposnorm(:,1))))/10);
skip = 50;
movement = NaN(size(linposnorm,1),1);
for j = 1:skip-1       
    dat = linposnorm(j:skip:end-j-1);
    movement(j:skip:end-j-1)=diff([dat;dat(end,:)]);                        
end
    
dirdat2 = NaN(size(linposnorm,1),1);
dirdat2(movement>0) = 2; %true; %new
dirdat2(movement<0) = 1; %false; %new

dirdat2(isnan(dirdat2)) = 1.5; %new
[FilterA] =fspecial('gaussian',[20 1],10); %new 30 10
Smoothed_Velocity=filtfilt(FilterA,1,dirdat2); %new    
vel = round(Smoothed_Velocity)==2; %new
dirdat = NaN(size(linposnorm));
dirdat(~isnan(linposnorm)) = vel(~isnan(linposnorm));


% if 0 
figure; hold on

plot(pos(dirdat==1,1),linposnorm(dirdat==1),'r.')
plot(pos(dirdat==0,1),linposnorm(dirdat==0),'b.')
plot(pos(isnan(dirdat),1),linposnorm(isnan(dirdat)),'k.')
    

%checks if this is okay
prompt = 'Is this okay? y/n';
title(prompt)
yn = input([prompt ':  '],'s');
if strcmp(yn,'n')
    return
else
    close gcf
end
% end