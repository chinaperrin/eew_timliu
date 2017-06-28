function ctlg = read_eq_catalog(catalogFullName,catalogType,opts)

%catalogFullName='/scratch/memeier/data/japan/jma/nps/mch/mch201601';
%catalogType    ='jmaFM';

% Catalog-type list
% jmaFM
% bombayBeachFM
% bombayBeachFM2

    
fid    = fopen(catalogFullName);
StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
fclose(fid);

eolLim      = '\r?\n';
LinePattern = ['[^\r\n]*', eolLim];
thelines    = regexp(StrRay,LinePattern,'match');
nlines      = numel(thelines);

if strcmp(catalogType,'jmaFM')
    
    % See format descriptions in 
    % .line 1: Hypocenter_record_format.pdf  
    % .line 2: Nodal_plane_solution_format.pdf
    neq = nlines/2;
    numFieldNames  = {'yr','mt','dy','hr','mn','sc','lat','lon','depth', ...
        'magnitude','strike','strike2','dip','dip2','rake','rake2'};
    for ifield =1:numel(numFieldNames);
        ctlg.(numFieldNames{ifield})=zeros(neq,1);
    end
    ctlg.mtype = cell(neq,1);
    
    for ieq=1:neq
        iline1 = 2*(ieq-1)+1;
        iline2 = 2*(ieq-1)+2;
        line1  = thelines{iline1};
        line2  = thelines{iline2};
        
        ctlg.yr(ieq) = str2double(line1(2:5));
        ctlg.mt(ieq) = str2double(line1(6:7));
        ctlg.dy(ieq) = str2double(line1(8:9));
        ctlg.hr(ieq) = str2double(line1(10:11));
        ctlg.mn(ieq) = str2double(line1(12:13));
        ctlg.sc(ieq) = str2double(line1(14:17))/100;
        
        ctlg.lat(ieq)   = str2double(line1(22:24)) + str2double(line1(25:28))/100/60;
        ctlg.lon(ieq)   = str2double(line1(33:36)) + str2double(line1(37:40))/100/60;
        ctlg.depth(ieq) = str2double(line1(45:49))/100 ;
        
        ctlg.magnitude(ieq) = str2double(line1(53:54))/10;
        ctlg.mtype{ieq}     = line1(55);
        
        ctlg.strike (ieq) = str2double(line2(41:44));
        ctlg.strike2(ieq) = str2double(line2(45:47));
        ctlg.dip    (ieq) = str2double(line2(48:52));
        ctlg.dip2   (ieq) = str2double(line2(53:56));
        ctlg.rake   (ieq) = str2double(line2(57:59));
        ctlg.rake2  (ieq) = str2double(line2(60:64));
    end
    
    
% EHs special format for the 2016 Bombay Beach Sequence
elseif strcmp(catalogType,'bombayBeachFM_fmt1')

    numFieldNames  = {'evid','yr','mt','dy','hr','mn','sc','lat','lon','depth', ...
        'magnitude','strike','dip','rake'};
    for ifield =1:numel(numFieldNames);
        ctlg.(numFieldNames{ifield})=zeros(nlines,1);
    end
    
    for iline = 1:nlines

        thisline = thelines{iline};
        
        ctlg.yr(iline)   = str2double(thisline(1:4));
        ctlg.mt(iline)   = str2double(thisline(5:6));
        ctlg.dy(iline)   = str2double(thisline(7:8));
        ctlg.hr(iline)   = str2double(thisline(10:11));
        ctlg.mn(iline)   = str2double(thisline(12:13));
        ctlg.sc(iline)   = str2double(thisline(15:19));
        
        latdeg            = str2double(thisline(21:22));
        latmin            = str2double(thisline(24:28));
        ctlg.lat(iline)   = latdeg+latmin/60;
        londeg            = -str2double(thisline(30:32));
        lonmin            =  str2double(thisline(34:38));
        ctlg.lon(iline)   = londeg-lonmin/60;
        ctlg.depth(iline) = str2double(thisline(41:45));
        
        ctlg.strike(iline) = str2double(thisline(113:115));
        ctlg.dip   (iline) = str2double(thisline(119:120));
        ctlg.rake  (iline) = str2double(thisline(122:125));
        
        ctlg.magnitude(iline) = str2double(thisline(49:52));
        ctlg.evid(iline)      = str2double(thisline(103:110));
    end
    
    
    
    % EHs special format for the 2001 & 2009 Bombay Beach Sequences
elseif strcmp(catalogType,'bombayBeachFM_fmt2')
    
    ctlg.head = thelines(1);
    thelines  = thelines(2:end);
    nlines    = numel(thelines);

    numFieldNames  = {'evid','yr','mt','dy','hr','mn','sc','lat','lon','depth', ...
        'magnitude','strike','dip','rake'};
    for ifield =1:numel(numFieldNames);
        ctlg.(numFieldNames{ifield})=zeros(nlines,1);
    end
    
    for iline = 1:nlines

        thisline = thelines{iline};
        
        ctlg.yr(iline)   = str2double(thisline(1:4));
        ctlg.mt(iline)   = str2double(thisline(5:6));
        ctlg.dy(iline)   = str2double(thisline(7:8));
        ctlg.hr(iline)   = str2double(thisline(10:11));
        ctlg.mn(iline)   = str2double(thisline(12:13));
        ctlg.sc(iline)   = str2double(thisline(15:19));
        
        latdeg            = str2double(thisline(21:22));
        latmin            = str2double(thisline(24:28));
        ctlg.lat(iline)   = latdeg+latmin/60;
        londeg            = -str2double(thisline(30:32));
        lonmin            =  str2double(thisline(34:38));
        ctlg.lon(iline)   = londeg-lonmin/60;
        ctlg.depth(iline) = str2double(thisline(41:45));
        
        %c---------------------------------------- Jeff Unruh -- Right hand rule ----
        %         if( (str_avg(i) .ge. 0.) .and. (str_avg(i) .le. 270.) )
        %      	idpdr1 = nint(str_avg(i)) + 90
        %         if( (str_avg(i) .gt. 270.) .and. (str_avg(i) .le. 360.) )
        %         idpdr1 =  nint(str_avg(i)) - 270
        dipdir = str2double(thisline(84:86));
        if     dipdir>=0   &dipdir<=270; strike=dipdir+90;
        elseif dipdir> 270 &dipdir<=360; strike=dipdir-270;
        else   fprintf(1,'sth must be wrong with strike. check.\n')
        end
        ctlg.strike(iline) = strike;
        ctlg.dip   (iline) = str2double(thisline(88:89));
        ctlg.rake  (iline) = str2double(thisline(90:93));
     
        ctlg.magnitude(iline) = str2double(thisline(49:52));
        ctlg.evid(iline)      = str2double(thisline(133:141));
    end
    
    
    
% EHs older special format for the 2016 Bombay Beach Sequence
% e.g. /scratch/memeier/data/socal/fm/bombaybeach/bb_2016_focal_mech.dat
elseif strcmp(catalogType,'bombayBeachFM')
    
    ctlg.head = thelines(1:30);
    thelines  = thelines(31:end);
    nlines    = numel(thelines);

    numFieldNames  = {'evid','yr','mt','dy','hr','mn','sc','lat','lon','depth', ...
        'magnitude','strike','dip','rake'};
    for ifield =1:numel(numFieldNames);
        ctlg.(numFieldNames{ifield})=zeros(nlines,1);
    end
    
    for iline = 1:nlines
        thisline = thelines{iline};
        
        ctlg.evid(iline) = str2double(thisline(1:8));
        ctlg.yr(iline)   = str2double(thisline(10:13));
        ctlg.mt(iline)   = str2double(thisline(15:16));
        ctlg.dy(iline)   = str2double(thisline(18:19));
        ctlg.hr(iline)   = str2double(thisline(21:22));
        ctlg.mn(iline)   = str2double(thisline(24:25));
        ctlg.sc(iline)   = str2double(thisline(27:31));
        
        ctlg.lat(iline)   = str2double(thisline(46:51));
        ctlg.lon(iline)   = str2double(thisline(53:60));
        ctlg.depth(iline) = str2double(thisline(65:68));
        ctlg.magnitude(iline) = str2double(thisline(38:41));
        
        ctlg.strike(iline) = str2double(thisline(113:115));
        ctlg.dip   (iline) = str2double(thisline(119:120));
        ctlg.rake  (iline) = str2double(thisline(122:125));
    end
    
    
% Relocated hypocenter catalog for SoCal: Hauksson, Yang & Shearer, 2011(?), BSSA(?)
elseif strcmp(catalogType,'hys2011')
    
    numFieldNames  = {'evid','yr','mt','dy','hr','mn','sc','lat','lon','depth','magnitude'};
    for ifield =1:numel(numFieldNames);
        ctlg.(numFieldNames{ifield})=zeros(nlines,1);
    end
    
    for iline = 1:nlines
        
        print_iteration_numbers(iline,nlines,'tenthousands')
        
        thisline = thelines{iline};
        
        ctlg.yr(iline)   = str2double(thisline(1:4));
        ctlg.mt(iline)   = str2double(thisline(6:7));
        ctlg.dy(iline)   = str2double(thisline(9:10));
        ctlg.hr(iline)   = str2double(thisline(12:13));
        ctlg.mn(iline)   = str2double(thisline(15:16));
        ctlg.sc(iline)   = str2double(thisline(18:23));
        
        ctlg.evid(iline) = str2double(thisline(25:33));
        
        ctlg.lat(iline)   = str2double(thisline(35:42));
        ctlg.lon(iline)   = str2double(thisline(44:53));
        ctlg.depth(iline) = str2double(thisline(55:61));
        ctlg.magnitude(iline) = str2double(thisline(64:67));
    end
    
% Relocated FM catalog for SoCal: Yang, Hauksson & Shearer, 2011, BSSA
elseif strcmp(catalogType,'yhs2011')

    numFieldNames  = {'evid','yr','mt','dy','hr','mn','sc','lat','lon','depth', ...
        'magnitude','strike','dip','rake'};
    for ifield =1:numel(numFieldNames);
        ctlg.(numFieldNames{ifield})=zeros(nlines,1);
    end
    
    for iline = 1:nlines
        
        print_iteration_numbers(iline,nlines,'tenthousands')

        thisline = thelines{iline};
        
        ctlg.yr(iline)   = str2double(thisline(1:4));
        ctlg.mt(iline)   = str2double(thisline(6:7));
        ctlg.dy(iline)   = str2double(thisline(9:10));
        ctlg.hr(iline)   = str2double(thisline(12:13));
        ctlg.mn(iline)   = str2double(thisline(15:16));
        ctlg.sc(iline)   = str2double(thisline(18:23));
        
        ctlg.evid(iline) = str2double(thisline(25:32));
        
        ctlg.lat(iline)   = str2double(thisline(34:42));
        ctlg.lon(iline)   = str2double(thisline(44:53));
        ctlg.depth(iline) = str2double(thisline(55:61));
        ctlg.magnitude(iline) = str2double(thisline(64:68));
        
        ctlg.strike (iline) = str2double(thisline(71:73));
        ctlg.dip    (iline) = str2double(thisline(76:77));
        ctlg.rake   (iline) = str2double(thisline(79:82));
    end

    
end