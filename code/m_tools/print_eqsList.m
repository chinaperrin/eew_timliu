function print_eqsList(eqs,fileName)

fprintf(1,'Include times, m-Types, etc.\n')
%pause

fid = fopen(fileName,'w');

% Sort eqs
[~,idx] = sort(eqs.m,'descend');
neq     = numel(eqs.m);

fprintf(fid,'Date \t Magnitude \t Lat \t Lon \t Depth \t No. records \t Event name\n');
for ieq = 1:neq
    
    ii   = (idx(ieq));
    
    m      = eqs.m(ii);
    %mType  = eqs.mType{ii};
    ntr    = numel(eqs.traceId{ii});
    eqDate = eqs.date{ii};
    eqName = eqs.name{ii};
    t0     = eqs.t0{ii};
    lat    = eqs.lat(ii);
    lon    = eqs.lon(ii);
    z      = eqs.z(ii);
    %line   = sprintf('%s \t %s \t %3.1f \t %s \t %5.2f \t %5.2f \t %5.1f \t %i',eqDate,t0,m,mType,lat,lon,z,ntr);
    line   = sprintf('%s \t %3.1f \t %5.2f \t %5.2f \t %5.1f \t %i \t %s',eqDate,m,lat,lon,z,ntr,eqName);
    fprintf(fid,'%s\n',line);
end
fclose(fid);







% fileID = fopen('namalentescht.txt','w');
% %fprintf(fileID,'%3s %12s\n','M','nrec');
% fprintf(fileID,'%i\t%i\n',NTR,NTR);
% fclose(fileID);
% 
% A1 = [9.9, 9900];
% A2 = [8.8,  7.7 ; ...
%       8800, 7700];
% formatSpec = '%4.1f%i\n';
% fprintf(formatSpec,M,NTR)
% 
% formatSpec = 'X is %4.2f meters or %8.3f mm\n';
% fprintf(formatSpec,NTR,NTR)
% 
% whos NTR
% fileID = fopen('entescht.txt','w');
% formatSpec = 'X is %4.2f meters or %8.3f mm\n';
% for ieq = 1:neq
%     fprintf(formatSpec,[ieq ieq])
% end
% fclose(fileID)
% 
% T = cell2table(list,'VariableNames',{'Magnitude','No records','Date'});
% writetable(list,'tabledata.dat')
% %[m,idx] = sort(eqs.Mw,'descend');
% 
% 
% date    = {eqs.eqDate{idx}}';
% nwf_nf  = eqs.nwf_NF(idx);
% nwf_ff  = eqs.nwf_FF(idx);
% fNetIdx = eqs.fNetIdx(idx);
% % 
% % m       = eqs.Mw;
% % date    = eqs.eqDate;
% % nwf_nf  = eqs.nwf_NF;
% % nwf_ff  = eqs.nwf_FF;
% % fNetIdx = eqs.fNetIdx;
% 
% 
% msg = 'This files contains the number of NF and FF traces. Some of the NF traces will have been skipped at a later stage, e.g. because they have too low sampling rates.';
% fprintf(fid, '%s \n\n', msg);
% 
% hdr = 'Date / Mw / No. traces in near-field / No. traces in far-field / fNetList-index';
% fprintf(fid, '%s \n', hdr);
% 
% neq = size(list.Mw,1);
% for ieq = 1:neq
%     %    fprintf(fid, '%s %3.1f %d %d %d \n', list.eqDate{ieq}, list.Mw(ieq),list.nwf_NF(ieq), list.nwf_FF(ieq), list.fNetIdx(ieq));
%     %fprintf(fid, '%s %3.1f %3.1f %3.1f %3.1f \n', date{ieq}, m(ieq), nwf_nf(ieq), nwf_ff(ieq), fNetIdx(ieq));
%     fprintf(fid, '%22s\t %9.1f\t %9d\t %9d\t %9d\n', date{ieq}, m(ieq), nwf_nf(ieq), nwf_ff(ieq), fNetIdx(ieq));
% end
% fclose(fid);
% 
% 
% 
% 
% %     % Print record origins
%     %     printListName = sprintf('out/recordList_m%s_to_m%s.txt',strrep(num2str(mRanges(ir,1),'%3.1f'), ...
%     %         '.','p'),strrep(num2str(mRanges(ir,2),'%3.1f'),'.','p'));
%     %     tmpList       = trList.selectSubList(idx{ir});
%     % 
%     %     if ir==1; fprintf(1,'Printing record fullnames to txt-file\n'); end
%     %     print_ascii_tracelist(tmpList,[],printListName)
%     %     
%     %     rec.lat = [rec.lat; tmpList.stationLat];
%     %     rec.lon = [rec.lon; tmpList.stationLon];
%     %     rec.m   = [rec.m  ; tmpList.m    ];
%     %     
%     %     [~,ii]  = unique(tmpList.eqLat);
%     %     eqk.lat = [eqk.lat; tmpList.eqLat(ii)];
%     %     eqk.lon = [eqk.lon; tmpList.eqLon(ii)];
%     % 
%     
%     
%     fprintf(1,'Printing record list for plotting map\n');
% recordListFileName = 'out/recList.txt';
% eventListFileName  = 'out/eventList.txt';
% dlmwrite(recordListFileName,[rec.lon, rec.lat, rec.m],'delimiter','\t')
% dlmwrite(eventListFileName, [eqk.lon, eqk.lat, eqk.m],'delimiter','\t')
