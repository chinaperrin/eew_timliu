function profile = read_bosai_depth_profile(fileFullName,opt)
 
% There are two file types:
% . with vs-values at regular 1m-intervals
% . with vs-values at irregular intervals

zMaster = (1:30)';
nMaster = numel(zMaster);

fid = fopen(fileFullName);
if opt.plotVsProfile; StrRay = fscanf(fid,'%c')
else                  StrRay = fscanf(fid,'%c');
end
fclose(fid);
eolLim      = '\r?\n';
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

fieldArray = cell(nlines,1);
for iline = 1:nlines
    thisLine = thelines{iline};
    fields   = regexp(thisLine,'\ ','split');
    fields   = fields(cellfun(@(x) ~isempty(x), fields));    
    fieldArray{iline} = fields; 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE TYPE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try and detect the '1m' field. If it is found, file is of file-type 1
idx = find(cellfun(@(x) strcmp(x(1),'1m'), fieldArray));
if ~isempty(idx); 
    is      = idx;
    %nfields = cell2mat(cellfun(@(x) size(x,2), fieldArray(is:end),'uniformOutput',0));
    nfields = cell2mat(cellfun(@(x) size(x,2), fieldArray,'uniformOutput',0));
    ie     = nlines;
    %ie      = find(nfields>=5,1,'last');
    nz      = ie-is+1;
    
    %     ucmd = sprintf('sed -n %i,1000p %s |awk ''{print $1}''',idx,fileFullName);
    %     [status,uOut] = unix(ucmd);
    %     a             = regexp(uOut,LinePattern,'match');
    %     depth         = str2double(cellfun(@(x) strrep(x,'m',''),a,'uniformOutput',0)');
    
%     depth = cell(nz,1);
%     vp    = cell(nz,1);
%     vs    = cell(nz,1);
%     for iline=is:ie
    depth = cell(nlines,1);
    vp    = cell(nlines,1);
    vs    = cell(nlines,1);
    for iline=idx:nlines
        if size(fieldArray{iline},2)>=4;
            depth(iline) = fieldArray{iline}(1);
            vp   (iline) = fieldArray{iline}(3);
            vs   (iline) = fieldArray{iline}(4);
        end
    end
    
    % Cut out lines that do not correspond to depth 
    aintEmpty = cellfun(@(x) ~isempty(x)      , vs);
    aintNan   = cellfun(@(x) ~strcmp(x,'Nan') , vs);
    aintTrunc = nfields>=5;
    useMe     = logical(aintEmpty.*aintNan.*aintTrunc);
    %useMe  = cellfun(@(x) ~isempty(x), vs);
    %lgcNan = 
    %useMe(isnan(vs)) = 0;
    profile.vs = str2double(vs(useMe));
    profile.vp = str2double(vp(useMe));
    profile.z  = str2double(strrep(depth(useMe),'m',''));
    
    profile.fileType = 1;
    if numel(unique(diff(profile.z)))~=1;
        fprintf(1,sprintf('8ung: depth profile does not have regular 1m intervals as assumed\n'))
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE TYPE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if nlines~=1
        nz            = nlines-2;
        depth = zeros(nz,1);
        vp    = zeros(nz,1);
        vs    = zeros(nz,1);
        for ii=1:nz
            iline = ii+2;
            depth(ii) = str2double(strrep(fieldArray{iline}(3),',',''));
            vp   (ii) = str2double(strrep(fieldArray{iline}(4),',',''));
            vs   (ii) = str2double(strrep(fieldArray{iline}(5),',',''));
        end
        
        % Find first depth value that is larger than 30m
        idx30 = find(depth>30);
        if isempty(idx30); idx30=numel(depth); end
        vsVect = zeros(nMaster,1);
        vpVect = zeros(nMaster,1);
        for ii=idx30:-1:1
            vsVect(zMaster<=depth(ii))=vs(ii);
            vpVect(zMaster<=depth(ii))=vp(ii);
        end
        profile.vs = vsVect;
        profile.vp = vpVect;
        profile.z  = zMaster;
        profile.fileType = 2;
        %clf; plot(vsVect,zMaster,'-xk'); set(gca,'ydir','reverse')
    else
        fprintf(1,'Depth profile file is empty\n')
        profile.fileType = 0;
    end
end