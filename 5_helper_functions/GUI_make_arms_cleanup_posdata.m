function [rawpos,params,vel,deltm] = GUI_make_arms_cleanup_posdata(rawpos,params)

%%%%%
%%%%%
%%%%% GUI to make position and velocity data and check that there aren't
%%%%% any problems with the position tracking/data
%%%%%
%%%%%

rawpos(rawpos(:,2)==0 & rawpos(:,3)==0,2:3) = NaN;
if ~isfield(params,'armslength') %not for grosmark data
rawpos(:,2) = fillmissing(rawpos(:,2),'linear');
rawpos(:,3) = fillmissing(rawpos(:,3),'linear');
end

figure; hold on
plot(rawpos(:,2),rawpos(:,3),'.k')
if ~isfield(params,'armslength') 
    xlim([0 max(rawpos(:,2))+50])
    ylim([0 max(rawpos(:,3))+50])
end

if ~isfield(params,'armslength')
% get the ends of the arms of a plus maze, y maze or linear track
if params.Track_Type == 2

    prompt = 'Type the length of that arm (in cm), press Enter';

    title('Select Middle of Maze, Press Enter')
    [xMiddle,yMiddle] = getpts; % press enter after final pt    
    x = NaN(params.numarms,1); y = x; len = x;
    for iarm = 1:params.numarms        
        title(['Select Arm ' num2str(iarm) ' (XY: 1:long/right 2:left 3:center/bottom) , Press Enter'])
        [x(iarm),y(iarm)] = getpts;            
        title(prompt)
        len(iarm) = input([prompt ':  ']);
    end
    arms = [repmat([xMiddle yMiddle],[params.numarms 1]) x y];
    params.arms = arms;
    params.armslength = len;

elseif params.Track_Type == 1

    prompt = 'Type the length of the linear track (in cm), press Enter';

    title('Select the bottom OR left of the linear track, Press Enter')
    [xbottom,ybottom] = getpts; % press enter after final pt    
    title('Select the top OR right of the linear track, Press Enter')
    [xtop,ytop] = getpts; % press enter after final pt    
    
    title(prompt)
    len = input([prompt ':  ']);
    
    params.arms = [xtop ytop;xbottom ybottom];
    params.armslength = len;

end
end

deltm = false(size(rawpos,1),1);
badpos = 1;
badtime = 1;
close all
while badtime
     figure; hold on     
    plot(rawpos(:,2),rawpos(:,3),'k.')
    figure; hold on
    a = subplot(2,1,1); hold on
    plot(rawpos(:,1),rawpos(:,2),'k.')
    b = subplot(2,1,2); hold on
    plot(rawpos(:,1),rawpos(:,3),'k.')
    prompt = 'Are the start and end times okay?';
    title(prompt)
    yn = input([prompt ':  '],'s');           
    if strcmp(yn,'n')                           
        title({'Click on the real start time,';'points before this point will be deleted'})
        [x,~] = getpts;
        yl = get(gca,'ylim');
        plot([x x],yl,'r--')
        plot(a,rawpos(rawpos(:,1)<x,1),rawpos(rawpos(:,1)<x,2),'b*')
        plot(b,rawpos(rawpos(:,1)<x,1),rawpos(rawpos(:,1)<x,3),'b*')
        prompt = 'Keep this time mask?';
        title(prompt)
        yn = input([prompt ':  '],'s');  
        if strcmp(yn,'y')
            deltm(rawpos(:,1)<x) = true;
            rawpos(rawpos(:,1)<x,:) = [];
        end
        title({'Click on the real end time,';'points after this point will be deleted'})
        [x,~] = getpts;
        yl = get(gca,'ylim');
        plot([x x],yl,'r--')
        plot(a,rawpos(rawpos(:,1)>x,1),rawpos(rawpos(:,1)>x,2),'b*')
        plot(b,rawpos(rawpos(:,1)>x,1),rawpos(rawpos(:,1)>x,3),'b*')
        prompt = 'Keep this time mask?';
        title(prompt)
        yn = input([prompt ':  '],'s');  
        if strcmp(yn,'y')
            deltm(rawpos(:,1)>x) = true;
            rawpos(rawpos(:,1)>x,:) = [];
        end                
    elseif strcmp(yn,'y')
        badtime = 0;
    end
end
        
while badpos == 1

    figure; hold on
    plot(rawpos(:,2),rawpos(:,3),'.k')
    xlim([min(rawpos(:,2))-50 max(rawpos(:,2))+50])
    ylim([min(rawpos(:,3))-50 max(rawpos(:,3))+50])

    prompt = 'Are there any bad position points? y/n?';
    title(prompt)
    yn = input([prompt ':  '],'s');
    if strcmp(yn,'y')
        badpos = 1;
    elseif strcmp(yn,'n')
        badpos = 0;
    end
    
    if badpos == 1
        xypt = NaN(3,2);
        title({'Click on top left point of a box to use as mask,';'points within that box will be deleted'})
        [xypt(1,1),xypt(1,2)] = getpts;
        plot(xypt(1,1),xypt(1,2),'+r','MarkerSize',15)
        title({'Click on top right point of a box to use as mask,'; 'points within that box will be deleted'})
        [xypt(2,1),~] = getpts;
        plot(xypt(2,1),xypt(1,2),'+r','MarkerSize',15)
        title({'Click on bottom left point of a box to use as mask,' ; 'points within that box will be deleted'})
        [~,xypt(3,2)] = getpts;
        plot(xypt(1,1),xypt(3,2),'+r','MarkerSize',15)
%         title('Click on bottom right point of a box to use as mask, points within that box will be deleted')
%         [xypt(4,1),xypt(4,2)] = getpts;
        plot(xypt(2,1),xypt(3,2),'+r','MarkerSize',15)

        xind = rawpos(:,2)<=xypt(2,1) & rawpos(:,2)>=xypt(1,1);
        yind = rawpos(:,3)<=xypt(1,2) & rawpos(:,3)>=xypt(3,2);
        ind = xind & yind;
        plot(rawpos(ind,2),rawpos(ind,3),'*b')
        clear xind yind

        prompt = 'Keep this mask? y/n';
        title(prompt)
        yn = input([prompt ':  '],'s');

        if strcmp(yn,'y')            
            rawpos(ind,2:3) = NaN;
            rawpos(:,2) = fillmissing(rawpos(:,2),'linear');
            rawpos(:,3) = fillmissing(rawpos(:,3),'linear');
        end        
    end
end

close all

vel = get_velocity(rawpos);



