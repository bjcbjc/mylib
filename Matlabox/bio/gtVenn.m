function [h, setnum] = gtVenn(data, varargin)
    % data: coded genotype, numbers
    %
    para.missing = -1;
    para.homref = 0;
    
    fds = fieldnames(para);    
    idx = (find(ismember(varargin( cellfun(@ischar, varargin)), fds))-1) * 2 + 1;    
    tf = false(size(varargin));
    tf(idx) = true;
    if any(tf)
        tf(idx+1) = true;
    end    
    passpara = varargin(~tf);
    varargin = varargin(tf);
    para = assignpara(para, varargin{:});
    
    setflagMap = containers.Map(2:4, ...
        { arrayfun(@(x) str2double(x), dec2bin(1:3,2)), ...
        arrayfun(@(x) str2double(x), dec2bin(1:7,3)), ...
        arrayfun(@(x) str2double(x), dec2bin(1:15,4)) } );
    
    nSample = size(data, 2);
    if nSample > 4
        error('not supported');
    end
    valid = any(data ~= para.missing & data ~= para.homref, 2);
    data = data(valid, :);
    setflag = setflagMap(nSample);
    vennData = zeros(size(setflag, 1), 1);
    ugt = setdiff(unique(data), [para.missing, para.homref]);
    
    for gtIdx = 1:length(ugt)
        vennIdx = binvec2dec( data == ugt(gtIdx)  );
        [count, vennIdx] = eleCounts(vennIdx(vennIdx ~= 0), false);
        vennData(vennIdx) = vennData(vennIdx) + count;
    end
    
    [h, setnum] = fixvenn(vennData, 'precalculated', true, 'setflag', setflag, passpara{:}); 
end
    

    
