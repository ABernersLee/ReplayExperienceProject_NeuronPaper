function get_midpoint(iident)

%%%%%
%%%%% gets midpoint of lap in order to determine stopping periods on either
%%%%% end
%%%%%

load(iident,'MidTime')


if ~exist('MidTime','var')

    load(iident,'laps_coverspace','laps_singlepass','params','pos',...
        'linposcat','behave_ind','behavior')


    if params.Track_Type==1
        laps = laps_coverspace;    
        MidTime = NaN(max(laps),2);
        mp = NaN(max(laps),2);
        cut1 = (max(linposcat(behavior==3))-min(linposcat(behavior==3)))*.56+min(linposcat(behavior==3));
        cut2 = (max(linposcat(behavior==3))-min(linposcat(behavior==3)))*.44+min(linposcat(behavior==3));
        for ilap = 1:max(laps)
            ind = find(laps==ilap & behave_ind==2 & linposcat<cut1 & linposcat>cut2);
            mp(ilap,1) = round((max(ind)-min(ind))/2+min(ind));
            ind = find(laps==ilap & behave_ind==6 & linposcat<cut1 & linposcat>cut2);
            mp(ilap,2) = round((max(ind)-min(ind))/2+min(ind));
            MidTime(ilap,1) = pos(mp(ilap,1),1);
            MidTime(ilap,2) = pos(mp(ilap,2),1);
        end
        
  

    elseif params.Track_Type==2

        laps = laps_singlepass;
        MidTime = NaN(max(laps),1);
        inmiddle = nansum(behavior==1,2); 
        mp = NaN(max(laps),1);
        for ilap = 1:max(laps)
            ind = find(laps==ilap & inmiddle,1,'first');
            mp(ilap,1) =  ind;
            MidTime(ilap,1) = pos(mp(ilap,1),1);
        end
    end


    vc = varycolor(nanmax(laps)); 
    vc = vc(randperm(size(vc,1)),:); 
    figure; hold on;
    for ilap = 1:max(laps)
         plot(pos(laps==ilap,1),linposcat(laps==ilap),'.','MarkerEdgeColor',vc(ilap,:))
    end
    plot(MidTime,linposcat(mp),'r*','MarkerSize',20)
    axis tight
    close gcf

    MidTime = sort(MidTime(:));

    save(iident,'MidTime','-append')
end


%%%%% to check     
load(iident,'pos','linposcat')
figure; hold on;
[~,I]=histc(MidTime,pos(:,1));
plot(pos(:,1),linposcat,'.k')   
plot(MidTime,linposcat(I),'r*','MarkerSize',20)
