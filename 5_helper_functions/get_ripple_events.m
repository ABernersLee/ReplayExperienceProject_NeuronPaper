function get_ripple_events(iident,dirs)

%%%%% Ripple - using average ripple power (zscored) across 3 tt with the most cells
%%%%% some of this code is edited from Brad's old ripple code

load(iident,'hp_cells','hpinterneurons','rawspikedata','params')
cellstouse = setdiff(hp_cells,hpinterneurons);
Tetrode = unique(rawspikedata(ismember(rawspikedata(:,2),cellstouse),3));

numcells = NaN(length(Tetrode),1);
for itt = 1:length(Tetrode)
    numcells(itt) = length(unique(rawspikedata(rawspikedata(:,3)==Tetrode(itt),2)));
end
[~,b] = sort(numcells,'descend');

% The following values are used to build the bandpass filter for ripple detection
Ripple_Stop_Low=130;
Ripple_Pass_Low=150;             % Most papers use between 150 and 250 Hz to identify ripples
Ripple_Pass_High=250;
Ripple_Stop_High=275;
Stop_Band_Attenuation_One=60;    % This was the default, I think.  
Pass_Band=1;                     % This was the default, I think.
Stop_Band_Attenuation_Two=80;    % This was the default, I think.
 
if contains(iident,'cpp') || contains(iident,'sal')     
     EEGname=dir([dirs.cscdatadir '\' iident(1:6) '*tt*' iident(6:end)]);
else
    EEGname=dir([dirs.cscdatadir '\' iident(1:end-4) '_csc_tt*']);    
end
load([dirs.cscdatadir '\' EEGname(1).name],'Fs')
Filter_Design_For_Ripple=fdesign.bandpass(Ripple_Stop_Low, Ripple_Pass_Low, Ripple_Pass_High, Ripple_Stop_High, Stop_Band_Attenuation_One, Pass_Band, Stop_Band_Attenuation_Two, Fs);
Ripple_Filter=design(Filter_Design_For_Ripple,'butter');  %'equiripple' and 'cheby1' and 'cheby2' also work, but 'butter' was generally faster and gave similar results

if length(b)>=3
    touse = Tetrode(b(1:3));
else
    touse = Tetrode(b);
end
if strcmp(iident,'Buran_20090708_Run4.mat') || strcmp(iident,'Discovery_20111124_Run1.mat') ...
        || strcmp(iident,'Harpy_20100115_Run1.mat') || strcmp(iident,'Harpy_20100115_Run2.mat') || ...
        strcmp(iident,'Harpy_20100115_Run3.mat')
    touse = Tetrode(b([1 3:4])); %maybe a csc was unchecked
end

for icsc = 1:length(touse)
   if contains(iident,'cpp') || contains(iident,'sal')
         if touse(icsc)<10
             EEGname=[dirs.cscdatadir '\' iident(1:6) 'tt0' num2str(touse(icsc)) iident(6:end)];
         else
             EEGname=[dirs.cscdatadir '\' iident(1:6) 'tt' num2str(touse(icsc)) iident(6:end)];
         end 
    elseif contains(iident,'Janni_2010408_Run') && str2double(iident(end-4))>1
        EEGname=[dirs.cscdatadir '\Janni_tt' num2str(touse(icsc)) iident(6:end-4)];               
   else
        EEGname=[dirs.cscdatadir '\' iident(1:end-4) '_csc_tt' num2str(touse(icsc))];    
   end
   load([EEGname],'cscsamples','csctimes','Fs')
    LFPdata = [csctimes cscsamples];       
    RipDat=zeros(size(LFPdata,1),3);
    RipDat(:,1)=LFPdata(:,1);
    RipDat(:,2)=filter(Ripple_Filter,LFPdata(:,2));
    RipDat(:,2)=RipDat(end:-1:1,2);   
    RipDat(:,2)=filter(Ripple_Filter,RipDat(:,2));
    RipDat(:,2)=RipDat(end:-1:1,2);
    for M=1:2000000:size(RipDat,1)  % In case the program crashes if much more than about 2000000 samples are transformed at once
        RipDat(M:min([size(RipDat,1),M+2000000]),3)=hilbert(RipDat(M:min([size(RipDat,1),M+2000000]),2));
    end
    RipDat(:,3)=abs(RipDat(:,3));
    % The following gaussian filter has a sigma of 12 ms
    Ripple_Gaussian_Filter=fspecial('gaussian',[round(7*(12.5/((1/Fs)*1000))),1],round(12.5/((1/Fs)*1000)));
    RipDat(:,3)=filtfilt(Ripple_Gaussian_Filter,1,RipDat(:,3));
    clear Ripple_Gaussian_Filter;
    % The following z-scores the filtered trace
    RipDat(:,3)=zscore(RipDat(:,3));
    if icsc==1
        HP_Ripple = RipDat(:,[1 3]);
        HP_RippleRaw = [RipDat(:,1) LFPdata(:,2) RipDat(:,2:end)];
    else
        HP_Ripple(:,2) = HP_Ripple(:,2)+RipDat(:,3);
    end
    disp(['Done with Ripple ' num2str(icsc) ' in get_ripple_events'])
end
HP_Ripple(:,2) = HP_Ripple(:,2)./length(touse); % time, zscored ripple power


%rip_thresh is how many std above the mean
disp('Start make_RPCandEvents_ripple')
VelThresh = 5;
maxlength = 1; 
returnto = 1; 
minlength = .05; 
minsep = .05;
rip_thresh = 2;

load(iident,'pos','vel')


HP_Ripple(HP_Ripple(:,1)<min(pos(:,1)) | HP_Ripple(:,1)>max(pos(:,1)),:) = [];
[~,~,Index] = histcounts(pos(:,1),HP_Ripple(:,1));
[~,~,Index2] = histcounts(HP_Ripple(:,1),pos(:,1));
k = find(Index<=0);
Index(k(k<length(Index)/2))=1;
Index(k(k>length(Index)/2))=max(Index);
k = find(Index2<=0);
Index2(k(k<length(Index2)/2))=1;
Index2(k(k>length(Index2)/2))=max(Index2);



meanRip = mean(HP_Ripple(Index(vel<VelThresh),2));
stdRip = std(HP_Ripple(Index(vel<VelThresh),2));

[~,Locations]=findpeaks(HP_Ripple(:,2),'MinPeakHeight',meanRip+(stdRip*rip_thresh));

% This eliminates peaks that occurred when the rat was moving
Locations=Locations(vel(Index2(Locations))<VelThresh);

% This finds the start and end timepoints for each event.
Ripple_Events=zeros(length(Locations),6);
Ripple_Events(:,3)=HP_Ripple(Locations,1);
for N=1:size(Ripple_Events,1)
    L=Locations(N);
    Ripple_Events(N,6) = L;
    while HP_Ripple(L,2)>(meanRip+(stdRip*returnto)) && L>1 %this finds the closest timepoint prior to the current peak that crosses the mean
        L=L-1;
    end
    Ripple_Events(N,1)=HP_Ripple(L,1);
    Ripple_Events(N,4) = L;
    L=Locations(N);
    while HP_Ripple(L,2)>(meanRip+(stdRip*returnto)) && L<size(HP_Ripple,1) %this finds the closest timepoint after the current peak that crosses the mean
        L=L+1;
    end
    Ripple_Events(N,2)=HP_Ripple(L,1);
    Ripple_Events(N,5) = L;
end

%get rid of weirdly long ones
Ripple_Events((Ripple_Events(:,2)-Ripple_Events(:,1))>maxlength,:) = [];

%combine close ones
for ind = 3:-1:1
    while sum(diff(Ripple_Events(:,ind))<minsep)>0 %ABL added
        merge = [diff(Ripple_Events(:,ind))<minsep;false];
        Ripple_Events(merge,[2 5]) = Ripple_Events(find(merge)+1,[2 5]);
        [~,tallest] = max([HP_Ripple(Ripple_Events(merge,6),2),HP_Ripple(Ripple_Events(find(merge)+1,6),2)],[],2);
        Ripple_Events(merge,[3 6]) = Ripple_Events(find(merge)+tallest-1,[3 6]);
        Ripple_Events(find(merge)+1,:) = [];
    end
end

%get rid of short ones that couldn't be combined even to make longer ones
Ripple_Events((Ripple_Events(:,2)-Ripple_Events(:,1))<minlength,:) = [];

Ripple_CandEvents = Ripple_Events;

disp(['Number of RP Candidate Events: ' num2str(size(Ripple_CandEvents,1))])

save(iident,'Ripple_CandEvents','HP_Ripple','HP_RippleRaw','-append')
