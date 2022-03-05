function [armpos,linposcat,linposnorm,linposcatnan,dirdat,cm_conv] = GUI_make_linear_position(pos,params)

%%%%%
%%%%% GUI that linearizes position data and assigns arms and checks to make
%%%%% sure it all looks good
%%%%%

if params.Track_Type==2

    %projects points onto the params.arms
    if params.numarms==3
        p = pos;
        x = p(:,2); y = p(:,3); 

        distance = NaN(size(p,1),size(params.arms,1));
        projpoint = distance;
        for k = 1:size(params.arms,1)
             len(k)=sqrt((params.arms(k,3)-params.arms(k,1))^2+(params.arms(k,4)-params.arms(k,2))^2);
             xlen=(params.arms(k,3)-params.arms(k,1))/len(k);
             ylen=(params.arms(k,4)-params.arms(k,2))/len(k);
             pp=[x y]*[xlen ylen]';
             ss=[params.arms(k,1) params.arms(k,2)]*[xlen ylen]';
             ee=[params.arms(k,3) params.arms(k,4)]*[xlen ylen]';
             frac=min(1,max(0,((pp-ss)./(ee-ss))));
             qx=params.arms(k,1)+(frac*(params.arms(k,3)-params.arms(k,1)));
             qy=params.arms(k,2)+(frac*(params.arms(k,4)-params.arms(k,2)));
             distance(:,k)=sqrt(((x-qx).^2)+((y-qy).^2));
             projpoint(:,k)=frac*len(k);
        end
        cl=cumsum(len);
        offset=[0 cl(1:end-1)];
        [~,a]=min(distance,[],2); 
        D=diag(ones(size(params.arms,1),1));
        po = sum((projpoint.*D(a,:)),2);
        posave = po;
    elseif params.numarms == 4

        p = bsxfun(@minus,pos(:,2:3),center);
        isarm = false(length(p),3);
        isarm(:,1) = abs(p(:,1))>abs(p(:,2));
        isarm(:,2:3) = bsxfun(@ge,p(:,1:2),[0 0]);

        %the y axis is flipped from the intuitive
        a = NaN(length(p),1);
        a(isarm(:,1)&isarm(:,2),1) = 1;
        a(isarm(:,1)&~isarm(:,2),1) = 3;
        a(~isarm(:,1)&~isarm(:,3),1) = 4;
        a(~isarm(:,1)&isarm(:,3),1) = 2;


        po = NaN(length(p),4);
        po(a==1,1) = abs(p(a==1,1));
        po(a==3,3) = abs(p(a==3,1));
        po(a==2,2) = abs(p(a==2,2));
        po(a==4,4) = abs(p(a==4,2));
        toadd = flip(cumsum([0 max(flip(po,2))]),2);
        p2 = bsxfun(@plus,po,toadd(2:end));
        p = [pos(:,1) p2];
    end
        
    al = params.armslength;

    %normalize each arm to 100% of itself
    po_new = po;
    cm_conv1 = NaN(size(al,1),1);
    for mii = 1:size(al,1)
       dat = po_new(a==mii);
       po(a==mii)=((dat-min(dat))./(max(dat)-min(dat)))*al(mii);
       po_new(a==mii)=((dat-min(dat))./(max(dat)-min(dat)))*100;
       cm_conv1(mii,1) = al(mii)./range(dat);
    end    
    po(po==0) = po(po==0)+.0001;   
 
    %normalizes each arm to armslength        
    offset = [0;cumsum(al(1:end-1))];
    p2 = NaN(size(po,1),size(params.arms,1));
    for mii = 1:size(params.arms,1)
        p2(a==mii,mii)=po(a==mii)+offset(mii);
    end
    cm_conv = mean(cm_conv1);

elseif params.Track_Type==1

    %get the dimention with the largest variablity
    [~,I] = max(range(pos(:,2:3)));    
    p2 = pos(:,I+1)-min(pos(:,I+1))+.001;
    p = pos(:,1:3); 
    
    
    if ~isfield(params,'armslength')
        params.armslength = 160;
    end
    x = pos(:,2); y = pos(:,3);     
    len=sqrt((params.arms(2,1)-params.arms(1,1))^2+(params.arms(2,2)-params.arms(1,2))^2);
    xlen=abs((params.arms(2,1)-params.arms(1,1)))/len;
    ylen=abs((params.arms(2,2)-params.arms(1,2)))/len;
    pp=[x y]*[xlen ylen]';
     ss=[params.arms(1,1) params.arms(1,2)]*[xlen ylen]';
     ee=[params.arms(2,1) params.arms(2,2)]*[xlen ylen]';
     frac=min(1,max(0,((pp-ss)./(ee-ss))));
     qx=params.arms(1)+(frac*(params.arms(2,1)-params.arms(1,1)));
     qy=params.arms(2)+(frac*(params.arms(2,2)-params.arms(1,2)));
     distance=sqrt(((x-qx).^2)+((y-qy).^2));
     projpoint=frac*len;
    

    %normalize arm to arms length
    dat = projpoint; 
%     dat = distance;
    p2=((dat-min(dat))./(max(dat)-min(dat)))*params.armslength;
    cm_conv = params.armslength./range(projpoint);
    
    p2(p2==0) = p2(p2==0)+.0001;   
    po_new = p2;
    p = pos(:,1:3);
    a = ones(size(p2,1),1);
end

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