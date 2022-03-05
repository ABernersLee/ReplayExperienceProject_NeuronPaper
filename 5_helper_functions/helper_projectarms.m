function [Dlinpos,armind] = helper_projectarms(p,params,Bin_Size)

x = p(:,1); y = p(:,2); 

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
[~,a]=min(distance,[],2); 
D=diag(ones(size(params.arms,1),1));
po = sum((projpoint.*D(a,:)),2);

% figure; hold on; plot(x(a==1),y(a==1),'r'); plot(x(a==2),y(a==2),'b'); plot(x(a==3),y(a==3),'g')

al = params.armslength;
for mii = 1:size(al,1)
   dat = po(a==mii);
   po(a==mii)=((dat-min(dat))./(max(dat)-min(dat)))*(al(mii)-2); %so they dont bleed      
end    
   

offset = [0;cumsum(al(1:end-1))];
p2 = NaN(size(po,1),size(params.arms,1));
for mii = 1:size(params.arms,1)
    p2(a==mii,mii)=po(a==mii)+offset(mii);    
end

p2(p2==0) = .001;
Dlinpos = ceil(nansum(p2,2)./Bin_Size);    
[a,b] = find(~isnan(p2)==1);
armind(a) = b;

% armind = false(max(Dlinpos),size(params.arms,1));
% for iarm = 1:size(params.arms,1)
%     DL = ceil(p2(:,iarm)./Bin_Size);   
%     armind(unique(DL(~isnan(DL))),iarm) = true;
% end
% figure; imagesc(armind)
% figure; hold on; plot(p2(:,1),'r'); plot(p2(:,2),'b'); plot(p2(:,3),'g')