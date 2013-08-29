function fileStruct = excludeFile(fileStruct, varargin)
    nfile = length(fileStruct);
    exclude = false(nfile, 1);
    fns = cell(nfile, 1);
    for i = 1:nfile
        fns{i} = fileStruct(i).name;
    end
    for i = 1:length(varargin)
        exclude = exclude | ~cellfun(@isempty, strfind(fns, varargin{i}));
    end
    fileStruct(exclude) = [];
end
