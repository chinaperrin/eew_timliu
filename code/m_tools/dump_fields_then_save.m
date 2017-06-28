overwrite_unneeded_fields(TraceList,{'cav'})
save(TraceListFullName ,'TraceList')

overwrite_unneeded_fields(SkipList,{'cav'})
save(SkipListFullName ,'SkipList')

%overwrite_unneeded_fields(SkipList,{'cav','nbpga','nbpgv','nbpgd','pgaIdx', ...
%    'pgvIdx','pgdIdx','mErr','rErr','tauP','tauC','nbnoise'})


% Finish log-file
unix(['echo "And run through it did." >> ' logFileName]);
unix(['echo "Other scripts that were used in this run: " >> ' logFileName]);
unix(['cat *.m >> ' logFileName]);
diary off
unix(['echo "########################## Matlab Output ###############################" >> ' logFileName]);
unix(['cat ',diaryFileName,' >> ' logFileName]);
unix(['rm ',diaryFileName]);
unix(['bzip2 -f ',logFileName]);
