%%%%%% Extracts position, spike, and lfp data

cd('D:/test/')
load('dirs_linear')
cd(dirs.paramdir)
d = dir;
sz = extractfield(d,'bytes')>0;
ftype = contains(extractfield(d,'name'),'.mat');
d2 = d(sz&ftype);

for iepoch = 1:size(d2,1)
    try
        cd(dirs.paramdir)
        load(d2(iepoch).name,'params')
        %extracts raw position and saves into each run, not processed yet
        %needs to go first, creates mat file that others append to
        params = extract_nvt(dirs,params);

        %extracts xclust spike data and saves into each run
        extract_spikes_from_xclust(dirs,params)

        %extract the cscs (only one per tetrode - theta/noise ratio)
        extract_cscs(dirs,params)
    catch ME
        disp(['ID: ' ME.identifier])    
        msgString = getReport(ME);
        disp(msgString)
        disp(iepoch)
        disp('Error Occured')
    end
end