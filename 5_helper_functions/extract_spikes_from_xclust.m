function extract_spikes_from_xclust(dirs,params)
cd(params.daydir)
numtetrodes = 40;
celldat = [];
CN = 0;

if ~strcmp(params.Rat_Name,'Chesapeake') && ~strcmp(params.Rat_Name,'Haw')

    cd(params.converted)

    %In Converted Folder, get tt folder names
    

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
            touse = ~contains(celldr,'bad') & ...
                    ~contains(celldr,'Bad') & ...
                    ~contains(celldr,'BAD') & ...
                    ~contains(celldr,'.PNG') & ...
                    ~contains(celldr,'Corrected') & ...
                    ~contains(celldr,'.mat') & ...
                    ~contains(celldr,'all') ;

            celldr = celldr(touse);
            for ic = 1:length(celldr)
                %load each cluster
                tempdat2 = load(celldr{ic});  
                CN = CN+1;
                if ~isempty(tempdat2)
                    
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


    %sort by time
    [~,ind] = sort(celldat(:,1));
    celld = celldat(ind,:);
    clear celldat

    %extract the run times
    st = (repmat(celld(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(celld(:,1))]);
    nd = (repmat(celld(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(celld(:,1))]);
    rundat = (st&nd)';
    clear st nd
    
    if 1
        rawspikedata = celld;
        save([dirs.spikedatadir '\' params.ident '.mat'],'rawspikedata','-append')
    else
    %save each run seperately
    for irun = 1:size(rundat,2)
        rawspikedata = celld(rundat(:,irun),:);
%         save([dirs.spikedatadir '\' params.Rat_Name '_' num2str(params.Date) '_Run' num2str(irun) '.mat'],'rawspikedata','-append')
        if isempty(rawspikedata)
            disp('Error')
        end
        save([dirs.spikedatadir '\' params.ident '.mat'],'rawspikedata','-append')
        clear rawspikedata
    end
    end

else
    Cdirs = [dir('*Converted');dir('*Converted_Run*')];    

    for irun = 1:size(params.Run_Times,1)
        celldat = [];
        CN = 0;

        for it = 1:numtetrodes

            if ismember(it,params.Othertt)
                cd Converted
            elseif ismember(it,params.HPtt)
                cd(Cdirs(irun+1).name)
            end

            if isfolder(['tt' num2str(it)])
                cd(['tt' num2str(it)])
                d = dir;
                d2 = extractfield(d,'name');   

                %in tt folder, get cluster names
                celldr =d2(~cellfun('isempty',strfind(d2,'cl-','ForceCellOutput',1)));

                % dont load photos, mat files, or files labed as bad
                touse = ~contains(celldr,'bad') & ...
                        ~contains(celldr,'.PNG') & ...
                        ~contains(celldr,'Corrected') & ...
                        ~contains(celldr,'.mat');

                celldr = celldr(touse);
                for ic = 1:length(celldr)
                    %load each cluster
                    tempdat2 = load(celldr{ic});      
                    CN = CN+1;
                    %add spike time, cell number, tetrode number, and spike widths
                    tempdat = [tempdat2(:,8) CN*ones(size(tempdat2,1),1) it*ones(size(tempdat2,1),1) tempdat2(:,6)]; 
                    clear tempdat2
                    %add to cell data matrix
                    celldat = cat(1,celldat,tempdat);
                    clear tempdat
                end   
                cd ../
            end
            cd ../
        end

        %sort by time
        [~,ind] = sort(celldat(:,1));
        celld = celldat(ind,:);
        clear celldat

         %extract the run times
        st = (repmat(celld(:,1),[1 size(params.Run_Times,1)]))'>=repmat(params.Run_Times(:,1),[1 length(celld(:,1))]);
        nd = (repmat(celld(:,1),[1 size(params.Run_Times,1)]))'<=repmat(params.Run_Times(:,2),[1 length(celld(:,1))]);
        rundat = (st&nd)';
        clear st nd

        rawspikedata = celld(rundat(:,irun),:);
        save([dirs.spikedatadir '\' params.ident '.mat'],'rawspikedata','-append')
        clear rawspikedata
    end
    
end
