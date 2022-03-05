function rs = extract_spikes_from_xclust_v2(params,clusterstring,daydir)
cd(daydir)
numtetrodes = 40;
    

%In Converted Folder, get tt folder names
% if no Converted Folder, go into each of the folders and pull out
d3 = dir;
if sum(contains(extractfield(d3,'name'),'Converted'))~=0
    ddirs = {'Converted'};
else
    touse1 = isfolder(extractfield(d3,'name'));
    touse1(1:2) = false;
    ddirs = extractfield(d3(touse1),'name');
end

rawspikedata = [];
for idir = 1:length(ddirs)
    cd([daydir '\' ddirs{idir}])
    celldat = [];
    for it = 1:numtetrodes            
        if isfolder(['tt' num2str(it)])
            cd(['tt' num2str(it)])
            d = dir;
            d2 = extractfield(d,'name');   

            ttstring = 'cl-';

            %in tt folder, get cluster names
            if isfield(params,'TTname')          
                if size(params.TTname,1)>1
                    celldr =d2(~cellfun('isempty',strfind(d2,ttstring,'ForceCellOutput',1)) & ...
                    (~cellfun('isempty',strfind(d2,params.TTname{1},'ForceCellOutput',1)) | ...
                    ~cellfun('isempty',strfind(d2,params.TTname{2},'ForceCellOutput',1))));
                else
                celldr =d2(~cellfun('isempty',strfind(d2,ttstring,'ForceCellOutput',1)) & ...
                    ~cellfun('isempty',strfind(d2,params.TTname,'ForceCellOutput',1)));
                end
            else
                celldr =d2(~cellfun('isempty',strfind(d2,ttstring,'ForceCellOutput',1)));
            end


            % dont load photos, mat files, or files labed as bad or
            % suspiciously long
            touse1 = ~contains(celldr,'bad') & ...
                    ~contains(celldr,'Bad') & ...
                    ~contains(celldr,'BAD') & ...
                    ~contains(celldr,'.PNG') & ...
                    ~contains(celldr,'Corrected') & ...
                    ~contains(celldr,'.mat') & ...
                    ~contains(celldr,'all') ;


            touse2 = contains(celldr,clusterstring{1});
            if sum(touse2)==0 && length(clusterstring)>1
                touse2 = contains(celldr,clusterstring{2});
            end
            touse = touse1 & touse2;
            
            celldr = celldr(touse);
            for ic = 1:length(celldr)
                %load each cluster
                tempdat2 = load(celldr{ic});  
%                 CN = CN+1;
                if ~isempty(tempdat2)
                    
                    CN = str2double(celldr{ic}(max(strfind(celldr{ic},'-'))+1:end))+it*100; % ttcl
                    
                    %add spike time, cell number, tetrode number, and spike widths
                    tempdat = [tempdat2(:,8) CN*ones(size(tempdat2,1),1) it*ones(size(tempdat2,1),1) tempdat2(:,6)]; 
                    clear tempdat2
                    
                    %add to cell data matrix
                    celldat = cat(1,celldat,tempdat);
                    clear tempdat
                end
            end   
            cd ../
        end
    end


    rawspikedata = cat(1,rawspikedata,celldat);
end

%sort by time
% rs = sort(rawspikedata,1);
[~,ind] = sort(rawspikedata(:,1));
rs = rawspikedata(ind,:);
% clear celldat

% take out clusters that aren't in all 3
st = (repmat(rs(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(rs(:,1))]);
nd = (repmat(rs(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(rs(:,1))]);
rundat = (st&nd)';
st1 = (repmat(rs(:,1),[1 size(params.Sleep_Times,1)]))'>=repmat(params.Sleep_Times(:,1),[1 length(rs(:,1))]);
nd1 = (repmat(rs(:,1),[1 size(params.Sleep_Times,1)]))'<=repmat(params.Sleep_Times(:,2),[1 length(rs(:,1))]);
sleepdat = (st1&nd1)';
clear st nd st1 nd1
e1_cells = unique(rs(rundat(:,1),2));
e2_cells = unique(rs(rundat(:,2),2));
s2_cells = unique(rs(sleepdat(:,1),2));
exclude = unique([setdiff(e1_cells,e2_cells); setdiff(e1_cells,s2_cells); setdiff(e2_cells,s2_cells)]);
rs(ismember(rs(:,2),exclude),:) = [];

%rename clusters to small numbers
[~,i] = histc(rs(:,2),unique(rs(:,2)));
rs(:,2) = i;
% figure; plot(rs(:,2),rs2(:,2),'.')
