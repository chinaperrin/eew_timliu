function [CMT] = import_cmtList(cmtListFullName,latMin,latMax,lonMin,lonMax)

% Load full CMT list; go through all lines; if event is within bounding
% box, extract values of interest into CMT-structure. Save CMT Structure.
%
% menandrin@gmail.com, 150412

% latMin = 25;
% latMax = 50;
% lonMin = 125;
% lonMax = 150;

if nargin<5
    latMin =  -90;
    latMax =   90;
    lonMin = -360;
    lonMax =  360;
end

cmtMatName = sprintf('%s_lat%i_%i_lon%i_%i.mat',cmtListFullName(1:end-4),latMin,latMax,lonMin,lonMax);

if exist(cmtMatName,'file')
    fprintf(1,'\nLoading existing cmt-list ... ')
    load(cmtMatName)
    fprintf(1,'done.\n')
else
    
    
    tenthousands = linspace(1e4,1e7,1e3);
    
    fprintf(1,'\nImporting new cmt-list ... ')
    fid    = fopen(cmtListFullName);
    StrRay = fscanf(fid,'%c');        % Read the entire text file into an array
    fclose(fid);
    fprintf(1,'done.\n')
    
    % End of line
    eolLim = '\r?\n';
    
    % Line
    LinePattern = ['[^\r\n]*', eolLim];
    thelines    = regexp(StrRay,LinePattern,'match');
    nl          = numel(thelines);
    
    % All events have five lines
    fives = 1:5:nl-4;   % List of first line numbers
    ieq   = 0;          % Event count
    
    neqmax     = nl/5;
    CMT.lat    = zeros(neqmax,1);
    CMT.lon    = zeros(neqmax,1);
    CMT.date   = cell (neqmax,1);
    CMT.t0     = cell (neqmax,1);
    CMT.depth  = zeros(neqmax,1);
    CMT.mb     = zeros(neqmax,1);
    CMT.MS     = zeros(neqmax,1);
    CMT.M0     = zeros(neqmax,1);
    CMT.strike = zeros(neqmax,1);
    CMT.dip    = zeros(neqmax,1);
    CMT.rake   = zeros(neqmax,1);
    
    % Go through all lines
    fprintf(1,'Going through all lines of list ...\n')
    for il = 1:nl
        
        if ismember(il,tenthousands) fprintf(1,sprintf('%i / %i\n',il,nl)); end
        
        % If current line is a first line of an event: find latitude and check
        % if it fits inside box
        if ismember(il,fives);
            
            line1 = thelines{il};
            fields1 = regexp(line1,'\s+', 'split');
            
            lat    = str2num(fields1{4});
            lon    = str2num(fields1{5});
            
            if ( (lat>=latMin) &&(lat<latMax) &&(lon>=lonMin) &&(lon<lonMax) )
                
                ieq = ieq+1;
                
                CMT.lat(ieq)   = lat;
                CMT.lon(ieq)   = lon;
                CMT.date{ieq}  = fields1{2};
                CMT.t0{ieq}    = fields1{3};
                CMT.depth(ieq) = str2num(fields1{6});
                CMT.mb(ieq)    = str2num(fields1{7});
                CMT.MS(ieq)    = str2num(fields1{8});
                
                %line2 = thelines{il+1};
                %line3 = thelines{il+2};
                
                % Line 4
                line4           = thelines{il+3};
                fields4         = regexp(line4,'\s+', 'split');
                momentExponent  = str2num(fields4{1});
                
                % Line 5
                line5           = thelines{il+4};
                fields5         = regexp(line5,'\s+', 'split');
                scalarMoment    = str2num(fields5{11});
                
                M0dyncm         = scalarMoment*10^momentExponent; % Seismic moment in [dyne-cm] (hallelujah!)
                CMT.M0(ieq)     = M0dyncm*1e-7;                   % Seismic moment in [Nm]
                
                CMT.strike(ieq) = str2double(fields5{12});
                CMT.dip(ieq)    = str2double(fields5{13});
                CMT.rake(ieq)   = str2double(fields5{14});
                
            end
        end
    end
    
    % Remove empty fields
    neq = ieq;
    CMT.lat    = CMT.lat(1:neq);
    CMT.lon    = CMT.lon(1:neq);
    CMT.date   = CMT.date(1:neq);
    CMT.t0     = CMT.t0(1:neq);
    CMT.depth  = CMT.depth(1:neq);
    CMT.mb     = CMT.mb(1:neq);
    CMT.MS     = CMT.MS(1:neq);
    CMT.M0     = CMT.M0(1:neq);
    CMT.strike = CMT.strike(1:neq);
    CMT.dip    = CMT.dip(1:neq);
    CMT.rake   = CMT.rake(1:neq);
    
    % Compute Moment Magnitudes
    CMT.Mw = moment2magnitude(CMT.M0);
    
    fprintf(1,'Saving CMT as mat-file ...\n')
    save(cmtMatName)
    
    fprintf(1,'done.\n')
end
