function basename = compare_basenames(basename_v,basename_h1,basename_h2)
% Find common basename. If there is no common part in strings, return empty.

[~,l1,commonName] = LCS(basename_v,basename_h1);
if (l1~=0)
    [~,l2,commonName2] = LCS(commonName,basename_h2);
    if (l2~=0)
        basename = commonName2;
    end
end

if ( (l1==0) || (l2=0) )
    basename = [];
end

