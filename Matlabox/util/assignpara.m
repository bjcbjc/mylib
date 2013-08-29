function [para, passthrough] = assignpara(para, varargin)
    %assign values to the existing fields of para
    %error if no such field in para
    %
    % 'PASSPARA': ignore assignment instead of error if option is unknown; 
    %       the passed arguments will be returned as passthrough
    %
    
    if isstruct(para)
        testfun = @isfield;
    elseif isobject(para)
        testfun = @isprop;
    else
        error('unknown data type for para');
    end
    
    sidx = find(cellfun(@ischar, varargin));
    ifpass = strcmpi(varargin(sidx), 'passpara');
    if any(ifpass)
        varargin(sidx(ifpass)) = [];
        ifpass = true;
    else
        ifpass = false;
    end
    passthrough = {};
    
    for i = 1:2:length(varargin)
        if testfun(para, varargin{i})
            para.(varargin{i}) = varargin{i+1};
        elseif ifpass
            passthrough = [passthrough; varargin(i:i+1)];
        else
            error('Unknown option %s',varargin{i});
        end
    end