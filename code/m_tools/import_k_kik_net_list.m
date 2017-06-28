function [stationList] = import_k_kik_net_list(ListFileName)

fid    = fopen(ListFileName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

% End of line
eolLim = '\r?\n';

% Line
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nst         = numel(thelines);

if (nst==0)
    error('Station list could not be read.')
end

stationList.name    = cell(nst,1);
stationList.stLat   = zeros(nst,1);
stationList.stLon   = zeros(nst,1);
stationList.stLat2  = zeros(nst,1);
stationList.stLon2  = zeros(nst,1);
stationList.stAlt   = zeros(nst,1);
stationList.stDepth = zeros(nst,1);
stationList.nw      = cell(nst,1);

for ist=1:nst
    
    ctLine  = thelines{ist};
    fields  = regexp(ctLine,'\s+', 'split');
    
    stationList.name{ist}    = fields{1};
    stationList.stLat(ist)   = str2double(fields{3});
    stationList.stLon(ist)   = str2double(fields{4});
    stationList.stLat2(ist)  = str2double(fields{8});
    stationList.stLon2(ist)  = str2double(fields{9});
    stationList.stAlt(ist)   = str2double(fields{5});
    stationList.stDepth(ist) = str2double(fields{6});
    stationList.nw{ist}      = 'jpsm';
    %stationList.nw{ist}      = fields{10};

end


% Clean out doubles from list
[~,sorted_idx,~]  = unique(stationList.name,'first');
ind_uniqueSt      = sort(sorted_idx);

stationList.name    = stationList.name(ind_uniqueSt);
stationList.stLat   = stationList.stLat(ind_uniqueSt);
stationList.stLon   = stationList.stLon(ind_uniqueSt);
stationList.stLat2  = stationList.stLat2(ind_uniqueSt);
stationList.stLon2  = stationList.stLon2(ind_uniqueSt);
stationList.stAlt   = stationList.stAlt(ind_uniqueSt);
stationList.stDepth = stationList.stDepth(ind_uniqueSt);
stationList.nw      = stationList.nw(ind_uniqueSt);

fprintf(1,[' ... ',num2str(numel(ind_uniqueSt)), ' out of ',num2str(nst),' stations in Japan      station list are unique\n'])