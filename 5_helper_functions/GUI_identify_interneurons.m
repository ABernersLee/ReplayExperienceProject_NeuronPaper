function [spikedata,hp_cells,other_cells,hpinterneurons] = GUI_identify_interneurons(rawspikedata,pos,params)

%%%%%
%%%%%
%%%%% GUI to seperate interneurons and check to make sure it looks right
%%%%%
%%%%%

%rawspikedata is time, cell #, tetrode, and spike width
hp_cells = unique(rawspikedata(ismember(rawspikedata(:,3),params.HPtt),2));
other_cells = unique(rawspikedata(ismember(rawspikedata(:,3),params.Othertt),2));

figure; 

rawspikedata(ismember(rawspikedata(:,2),hp_cells),4);

h = histcounts(rawspikedata(:,2),[hp_cells-.5;hp_cells(end)]);
hh = NaN(length(hp_cells),1);
for icell = 1:length(hp_cells)
    hh(icell,1) = mean(rawspikedata(rawspikedata(:,2)==hp_cells(icell),4));
end

% define interneruons by spike width and number of spikes 
done = 0;
while done == 0
    figure; hold on
    plot(hh,h,'*k')
    yl = get(gca,'ylim');
    xlabel('Spike Width')
    ylabel('# of Spikes')
    plot([11 11],yl,'r--')    
    text(11.1,yl(2)*.9,'Ting''s cutoff','Color','r')

    prompt = 'Are there any interneurons? y/n';
    title(prompt)
    yn = input([prompt ':  '],'s');

    if strcmp(yn,'n')            
        hpinterneurons = [];
        close gcf
        done = 1;
    end

    if done == 0
        xypt = NaN(3,2);
        title({'Click on top left point of a box to use as mask,';'points within that box will be interneurons'})
        [xypt(1,1),xypt(1,2)] = getpts;
        plot(xypt(1,1),xypt(1,2),'+r','MarkerSize',15)
        title({'Click on top right point of a box to use as mask,'; 'points within that box will be interneurons'})
        [xypt(2,1),~] = getpts;
        plot(xypt(2,1),xypt(1,2),'+r','MarkerSize',15)
        title({'Click on bottom left point of a box to use as mask,' ; 'points within that box will be interneurons'})
        [~,xypt(3,2)] = getpts;
        plot(xypt(1,1),xypt(3,2),'+r','MarkerSize',15)
        plot(xypt(2,1),xypt(3,2),'+r','MarkerSize',15)

        xind = hh<=xypt(2,1) & hh>=xypt(1,1);
        yind = h<=xypt(1,2) & h>=xypt(3,2);
        ind = xind & yind';
        plot(hh(ind),h(ind),'*b')

        prompt = 'Is this right? y/n';
        title(prompt)
        yn = input([prompt ':  '],'s');

        if strcmp(yn,'y')            
            hpinterneurons = hp_cells(ind);
            close gcf
            done = 1;
        end
    end
end

% spiketoposindex 

bks = [0; find(diff(pos(:,1))>60); size(pos(:,1),1)];
spikedata = [rawspikedata(:,1:2) NaN(size(rawspikedata,1),1)];
for ibk = 1:size(bks,1)-1
    rind = rawspikedata(:,1)>=pos(bks(ibk)+1,1) & rawspikedata(:,1)<=pos(bks(ibk+1),1);
    rs = rawspikedata(rind,1);
    Time = pos(pos(:,1)>pos(bks(ibk)+1,1) & pos(:,1)<=pos(bks(ibk+1),1),1);
    [~,i] = histc(rs(:,1),[Time(1)-.001; Time(1:end-1)+diff(Time)/2 ; Time(end)+0.001]);
    spikedata(rind,3) = i+ bks(ibk);
%     figure; histogram(i)
end
tothrow = isnan(spikedata(:,3)) | spikedata(:,3)==0;
notthrow = false(size(tothrow));
for isleep = 1:size(params.Sleep_Times,1)
    notthrow = spikedata(:,1)>=params.Sleep_Times(isleep,1) & spikedata(:,1)<=params.Sleep_Times(isleep,2) | notthrow;
end
spikedata(tothrow & (~notthrow),:) = [];

close gcf