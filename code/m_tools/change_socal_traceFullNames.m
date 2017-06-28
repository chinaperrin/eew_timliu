function trList = change_socal_traceFullNames(trList)
% Replace all traceFullName with the ones rel. to the new directory
% structure. If it is not found, remove trace from list.

ntr = numel(trList.m);
for itr = 1:ntr
    
    fullName = trList.fullName{itr};
    isSocal  = ~isempty(regexp(fullName,'scsn_','once'));
    
    if isSocal
    
        m = trList.m(itr);
        if     m>=4 &m<5; mDirName = 'M4';
        elseif m>=5 &m<6; mDirName = 'M5';
        elseif m>=6     ; mDirName = 'M6p';
        end
    
        slashIdx    = regexp(fullName,'/');
        midDirName  = fullName(slashIdx(5)+1:slashIdx(6)-1);
        newFullName = strrep(fullName,midDirName,mDirName);    
        
        trList.fullName{itr} = newFullName;
        %         ucmd      = sprintf('find /scratch/memeier/data/socal/M* -iname %s',traceName);
        %         [status,newFullName] = unix(ucmd);
        %         if status==0 &~strcmp(newFullName,''); trList.fullName{itr} = newFullName
        %         else                                   flg_skipMe(itr) = true;
        %         end
    end
end

