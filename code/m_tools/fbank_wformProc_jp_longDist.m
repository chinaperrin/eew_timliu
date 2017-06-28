function fid=print_asciiTable_tPGx(printList,asciiTableFullName)

fid = fopen(asciiTableFullName,'w');
    %ntr = numel(zList.m)

    headerLine = sprintf('# Date\n# Time\n# Mag.\n# Hyp. Dist. [km]\n# Lat\n# Lon\n# PGA [m/s/s]\n# PGV [m/s]\n# PGD [m]\n# t(PGA)\n# t(PGV)\n# t(PGD)\n# Data Set\n');
    fprintf(fid,headerLine);
    
    printList = zList;
    ntr = numel(printList.m);
    for itr = 1:ntr
        
        print_iteration_numbers(itr,ntr,'thousands')
        
        eqDate = printList.eqDate{itr};
        t0     = printList.t0{itr};
        eqLat  = printList.eqLat(itr);
        eqLon  = printList.eqLon(itr);
        dsn    = printList.dataSetName{itr};
        oc     = printList.orntCode{itr};
        ic     = printList.instrCode{itr};
        m      = printList.m(itr);
        hd     = printList.hypDist(itr);
        ed     = printList.epiDist(itr);
        dsName = printList.dataSetName{itr};
        pga    = printList.pga(itr);
        pgv    = printList.pgv(itr);
        pgd    = printList.pgd(itr);
        %         tpga   = (printList.pgaIdx(itr)-printList.ppxIdx(itr))/printList.sRate(itr);
        %         tpgv   = (printList.pgvIdx(itr)-printList.ppxIdx(itr))/printList.sRate(itr);
        %         tpgd   = (printList.pgdIdx(itr)-printList.ppxIdx(itr))/printList.sRate(itr);
        tpga   = printList.pgaIdx(itr)/printList.sRate(itr);
        tpgv   = printList.pgvIdx(itr)/printList.sRate(itr);
        tpgd   = printList.pgdIdx(itr)/printList.sRate(itr);
        
        printLine = sprintf('%s\t%s\t%4.2f\t%6.2f\t%9.4f\t%9.4f\t%e\t%e\t%e\t%5.2f\t%5.2f\t%5.2f\t%s\n',eqDate,t0,m,hd,eqLat,eqLon,pga,pgv,pgd,tpga,tpgv,tpgd,dsName);
        %printLine = sprintf('%i\n',itr);
    
        %fprintf(1,printLine)
        fprintf(fid,printLine);
    end
    fclose(fid);