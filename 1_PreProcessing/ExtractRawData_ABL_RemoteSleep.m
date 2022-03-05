%  bring Delia's multiple novel tracks with sleep in between in to look at whether remote replays 
%  replays get longer without experience


%extracts the raw data
%run RunReplay_ABL_RemoteSleep.m after

%% used to copy over xclust cbfiles
clearvars
daynum = 13;
hd2 = 'D:\test\RemoteSleepData\rawdata\';
cd(hd2)
homedir1 = dir('*0*'); %Atlantis_2009-06-16-s1e1s2';
homedir = [hd2 homedir1(daynum).name];
cd(homedir)
r1 = [homedir '\e1e2'];
cd(r1)
thefile = 'cbfile-e1e2';

r2 = [homedir '\e1'];
foldertomovefrom = r1;
foldtomoveto = r2;

cd(foldertomovefrom)
d2 = dir('*tt*');
for itt = 1:size(d2,1)
   cd(d2(itt).name)
   l = dir;
   if sum(strncmp(extractfield(l,'name'),thefile,length(thefile)))>0 && ...
       isfolder([foldtomoveto '\' d2(itt).name '\'])
       copyfile(thefile,[foldtomoveto '\' d2(itt).name '\'])
       disp([foldtomoveto '\' d2(itt).name '/'])
   end
   cd ../
end
%% used to bring raw neuralynx data into xclust format
clearvars
hd2 = 'D:\test\RemoteSleepData\';
% foldnum = [1 2];
cd(hd2)
homedir1 = dir('*0*'); %Atlantis_2009-06-16-s1e1s2'; dont use 7!!!
hd = [3 4 10];
for id = length(hd)    
    homedir = [hd2 homedir1(hd(id)).name];
    cd(homedir)
    if id == 1
        params.Run_Times = [31407741139 32463062832; 36789551698 37731398450]./1e6;
        params.Sleep_Times = [33540774202 36472196572]./1e6;
    elseif id == 2
        params.Run_Times = [3100876282 4017510042; 7028481844 7817761481]./1e6;
        params.Sleep_Times = [4364209343 6895265347]./1e6;
    elseif id == 3
        params.Run_Times = [7592221699 9992877855; 13543318694 15750989496]./1e6;
        params.Sleep_Times = [10819725740 13320945501]./1e6;
    end
    CONVERT_NEURALYNX_TO_XCLUST_ABL(1:40,'VT1.nvt',min(min(params.Run_Times)),max(max(params.Run_Times)))
    %no labels, just cbfile and cl-
end
%% run on all of them to pull out position data and spike data

clearvars
hd2 = 'D:\test\RemoteSleepData\rawdata\';
matdir = 'D:\test\RemoteSleepData\';

cd(hd2)
d2 = dir('*0*'); 
for id = 1:size(d2,1)
    if id==11
        continue %turned out it was the same track twice
    end
    Date = str2double(d2(id).name(end-7:end));
    Rat_Name = d2(id).name(1:end-9);
    params.Date = Date;
    params.Rat_Name = Rat_Name;
    if id==1        
        clusterstring = {'s1e1s2-'};
        params.Run_Times = [10948126679 21019309309;23447063842 25794618821]./1e6;
        params.Sleep_Times = [21154882134 23343563467]./1e6;
    elseif id == 2       % in two folders 
        clusterstring = {'s1e1s2-v2-';'s1e1s2-'};        
        params.Run_Times = [13039379317 14592475821; 17599765943 18700919423]./1e6;
        params.Sleep_Times = [14651225432 17122632663]./1e6;
    elseif id == 3
        clusterstring = {'cl-'};
        params.Run_Times = [31407741139 32463062832; 36789551698 37731398450]./1e6;
        params.Sleep_Times = [33540774202 36472196572]./1e6;
    elseif id == 4
        clusterstring = {'cl-'};
        params.Run_Times = [3100876282 4017510042; 7028481844 7817761481]./1e6;
        params.Sleep_Times = [4364209343 6895265347]./1e6;
    elseif id == 5        
        clusterstring = {'e1e2-v2-';'e1e2-'};
        params.Run_Times = [15964483927 17957190329; 21841493578 24449*1e6]./1e6;
        params.Sleep_Times = [18822758648 21448186927]./1e6;
    elseif id == 6        
%         params.Run_Times = [6458143874 6475916466;11066637403 12101358616]./1e6;
%         params.Sleep_Times = [7449704924 10919667009]./1e6;
        params.Run_Times = [11067*1e6 12097*1e6;14897*1e6 16449*1e6]./1e6;
        params.Sleep_Times = [12233*1e6 14781*1e6]./1e6;
    elseif id == 7        
        params.Run_Times = [3965730112 6913016117; 10060681351 12376357371]./1e6;
        params.Sleep_Times = [7854655517 10020852518]./1e6;
    elseif id == 8        
        params.Run_Times = [5434780280 7746651870; 11478080816 13425251369]./1e6;
        params.Sleep_Times = [8698629382 11398066322]./1e6;
    elseif id == 9
        clusterstring = {'cl-'};
        params.Run_Times = [7592221699 9992877855; 13543318694 15750989496]./1e6;
        params.Sleep_Times = [10819725740 13320945501]./1e6;
    elseif id == 10
        params.Run_Times = [3768355141 5882139292; 8731108524 10201094579]./1e6;
        params.Sleep_Times = [6053715349 8494897384]./1e6;
    elseif id == 11
        params.Run_Times = [4815321287 7193405613; 10835787290 13588498059]./1e6;
        params.Sleep_Times = [8340349531 10708594687]./1e6;
     elseif id == 12
%         params.Run_Times = [2572675747 4676157644; 7979580161 10492595424]./1e6;
%         params.Sleep_Times = [4739291788 6932856972]./1e6;
        params.Run_Times = [4760*1e6 6930*1e6; 10570*1e6 13156*1e6]./1e6;
        params.Sleep_Times = [7980*1e6 10490*1e6]./1e6;
    elseif id == 13
%         params.Run_Times = [4240349859 5019572964; 7294168476 10208600299]./1e6;
%         params.Sleep_Times = [5022940677 6452202155]./1e6;   
        params.Run_Times = [4240*1e6 6450*1e6; 10270*1e6 13400*1e6]./1e6;
        params.Sleep_Times = [7294*1e6 10206*1e6]./1e6;   
    elseif id == 14
        params.Run_Times = [13000*1e6 13715*1e6; 3800*1e6 4325*1e6]./1e6; %from brad's notes and the Cheetah file
        params.Sleep_Times = [213170979 3800*1e6]./1e6;
    end
    
    if id > 5 && id~=9
        clusterstring = {'e1e2'};
    end
        
    daydir = [Rat_Name '_' num2str(Date)];
    save([matdir daydir '.mat'],'params','-append')
    if id ~= 2 && id~=14
       rawpos = extract_nvt_v2(params,[hd2 daydir]);
    else
        cd([hd2 daydir])
        ddd = dir;
        ddd = ddd(3:end);
        rawpos = [];
        for iddd = 1:2 %length(ddd)
            [rawpos1] = extract_nvt_v2(params,[hd2 daydir '\' ddd(iddd).name]);
            rawpos = cat(1,rawpos,rawpos1);
        end
    end
    
   rawspikedata = extract_spikes_from_xclust_v2(params,clusterstring,[hd2 daydir]);
   
   if id==14
        indd = find(rawpos(:,5)==2,1,'last');        
        rawpos(indd+1:end,:) = [];
        tadd = max(rawpos(:,1));
        tind = rawpos(1,1);
        
        rawspikedata(rawspikedata(:,1)<=tind,1) = rawspikedata(rawspikedata(:,1)<=tind,1)+tadd;
        rawpos(rawpos(:,1)<=tind,1) = rawpos(rawpos(:,1)<=tind,1)+tadd;
        params.Run_Times(2,:) = params.Run_Times(2,:) + tadd;
        params.Sleep_Times = params.Sleep_Times + tadd;
    end
   
   save([matdir daydir '.mat'],'params','rawpos','rawspikedata','-append')
   disp(['Done ' daydir ' - ' num2str(id)])
end

% params.Run_Times = []./1e6;
% params.Sleep_Times = []./1e6;
%%
