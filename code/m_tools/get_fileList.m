function fileList = get_fileList(dirName,searchPattern)

fileList   = dir(dirName);
fileList   = {fileList.name}';
if ~isempty(searchPattern); 
    hasPattern = cellfun(@(x) ~isempty(x), regexp(fileList,searchPattern));
    fileList   = fileList(hasPattern);
end


