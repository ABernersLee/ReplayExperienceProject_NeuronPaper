function vel = get_velocity(pos)

%velbin is the time bin that I use and then interpolate from 
velbin = 2;

[FilterA] =fspecial('gaussian',[100 1],10); 


pre_thinned_dat = pos;
skip = floor(velbin*(length(pre_thinned_dat)/range(pre_thinned_dat(:,1)))); %1 second * velbin (2 seconds here)
pre_thinned_dat(:,4:6) = NaN(size(pre_thinned_dat,1),3);
for j = 1:skip+1
    endind = size(pre_thinned_dat,1)-j-1;
    dat = pre_thinned_dat(j:skip:endind,:);
    Time_Change1=diff(dat(:,1));
    Time_Change1(end+1)=Time_Change1(end);            
    Time_Change1(Time_Change1<=0)=min(Time_Change1(Time_Change1>0))/10;
    X_Movement=diff(dat(:,2));
    X_Movement(end+1)=X_Movement(end);
    Y_Movement=diff(dat(:,3));
    Y_Movement(end+1)=Y_Movement(end);
    Distance_Moved1=sqrt((X_Movement.^2)+(Y_Movement.^2));
    Raw_Velocity=abs(Distance_Moved1./Time_Change1);        
    pre_thinned_dat(j:skip:endind,4) = Raw_Velocity;
    pre_thinned_dat(j:skip:endind,5) = Time_Change1;
    pre_thinned_dat(j:skip:endind,6) = Distance_Moved1;
    if j == skip+1
        lasts = find(isnan(pre_thinned_dat(:,4)));
%                 for i = 1:length(lasts)
        dat = pre_thinned_dat(lasts(1):end,:);
        Time_Change1=diff(dat(:,1));
        Time_Change1(end+1)=Time_Change1(end);            
        Time_Change1(Time_Change1<=0)=min(Time_Change1(Time_Change1>0))/10;
        X_Movement=diff(dat(:,2));
        X_Movement(end+1)=X_Movement(end);
        Y_Movement=diff(dat(:,3));
        Y_Movement(end+1)=Y_Movement(end);
        Distance_Moved1=sqrt((X_Movement.^2)+(Y_Movement.^2));
        Raw_Velocity=abs(Distance_Moved1./Time_Change1);        
        pre_thinned_dat(lasts(1):end,4) = Raw_Velocity;
        pre_thinned_dat(lasts(1):end,5) = Time_Change1;
        pre_thinned_dat(lasts(1):end,6) = Distance_Moved1;
%                 end
    end

end           
pre_thinned_dat(end,:) = pre_thinned_dat(end-1,:);
vel=filtfilt(FilterA,1,pre_thinned_dat(:,4));
