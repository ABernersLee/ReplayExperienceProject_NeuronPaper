function I = getfirstlap(pos,Pos,MidTime)

%%%%% Gets first lap index for looking at how hovers appear on the first
%%%%% lap in Hovers_Moving_laps.m for sup fig 6.

numlaps = 2; % length(MidTime)-1;
[~,I4]=histc(pos(:,1),[MidTime(1:end)]);
%same as in decode_spikedensity_events_lapbylap

cutoff1 = range(Pos)*(9/10)+min(Pos); %top
cutoff2 = (range(Pos)*(1/10))+min(Pos); % bottom
LapTime = min(pos(:,1));
Ilap = false(size(pos,1),numlaps);
I = Ilap;
Iup = NaN(size(pos,1),numlaps);
for ilap = 1:2      
    if sum(Pos(I4==ilap)<=cutoff2)>0
        %rat just ran down to the bottom on this lap (down)
        endlap = mean(pos(I4==ilap & Pos<=cutoff2));
        Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
        Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 0;
    elseif sum(Pos(I4==ilap)>=cutoff1)>0
        %rat just ran up to the top on this lap (up)
        endlap = mean(pos(I4==ilap & Pos>=cutoff1));
        Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
        Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 1;
    elseif sum(Pos(I4==ilap)>=cutoff1)==0 && sum(Pos(I4==ilap)<=cutoff2)==0
        cutoff11 = range(Pos)*(8/10)+min(Pos); %top
        cutoff22 = (range(Pos)*(2/10))+min(Pos); % bottom
         if sum(Pos(I4==ilap)<=cutoff22)>0
            %rat just ran down to the bottom on this lap (down)
            endlap = mean(pos(I4==ilap & Pos<=cutoff22));
            Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
            Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 0;
         elseif sum(Pos(I4==ilap)>=cutoff11)>0
            %rat just ran up to the top on this lap (up)
            endlap = mean(pos(I4==ilap & Pos>=cutoff11));
            Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
            Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 1;
         elseif sum(Pos(I4==ilap)<=cutoff22)==0 && sum(Pos(I4==ilap)>=cutoff11)==0
             badlap = 1; n1 = 8; n2 = 2;
             while badlap == 1
                 n1 = n1-1; n2 = n2+1;
                 cutoff11 = range(Pos)*(n1/10)+min(Pos); %top
                 cutoff22 = (range(Pos)*(n2/10))+min(Pos); % bottom
                  if sum(Pos(I4==ilap)<=cutoff22)>0
                    %rat just ran down to the bottom on this lap (down)
                    endlap = mean(pos(I4==ilap & Pos<=cutoff22));
                    Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
                    Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 0;
                    badlap = false;
                 elseif sum(Pos(I4==ilap)>=cutoff11)>0
                    %rat just ran up to the top on this lap (up)
                    endlap = mean(pos(I4==ilap & Pos>=cutoff11));
                    Ilap(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = true; 
                    Iup(pos(:,1)>=LapTime(end) & pos(:,1)<endlap,ilap) = 1;
                    badlap = false;
                  end
                  disp(['Bad Lap ' num2str(ilap)])
             end
             
         end        
    end
    LapTime = cat(1,LapTime,endlap);
    
    I(:,ilap) = Ilap(:,ilap);
end