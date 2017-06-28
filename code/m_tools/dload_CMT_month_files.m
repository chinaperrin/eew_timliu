
Note: does not fully work becuase not all months are labelled consistently... 
    June is sometimes jun03 and sometimes june03.
I have manually downloaded the files instead.

targetDir = '/scratch/memeier/data/cmt/monthly/';
sourceDir = 'http://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_MONTHLY/';
startYear = 5;
endYear   = 14;

monthNames = {'jan';'feb';'mar';'apr';'may';'june';'july';'aug';'sept';'oct';'nov';'dec';};
nm         = numel(monthNames);

wget_command_fileName = '/scratch/memeier/data/cmt/monthly/wget_commands.sh';
fid                   = fopen(wget_command_fileName,'w');


for iy = startYear:endYear
    
    for im = 1:nm
        
        fileName       = sprintf('%s%02d.ndk',monthNames{im},iy);
        fullName       = sprintf('%s20%02d/%s',sourceDir,iy,fileName);
        targetFullName = sprintf('%s%s',targetDir,fileName);
        
        % try downloading file
        ucmd = sprintf('wget -O%s %s',targetFullName,fullName);
        %unix(ucmd)
        fprintf(fid,sprintf('%s\n',ucmd))
        
    end
end
fclose(fid);
