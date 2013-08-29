function res = geneInMarker(markers, window, varargin)
    %markers: cell arrays of marker strings
    %window: window size, default = 10k
    %SGD: from SGD.mat
    %'mode': {'complete','partial'}; def='partial'
    %'name': {'orf','gene','index'}; def='orf'
    %'type': {'Verified','all'}; def='all'
    %
    %res: cell array, each contains a cell array of orfs in the marker
    %    
    if nargin < 2, window = 10000; end            
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'genelocdb'
                SGD = varargin{i+1};
            case 'mode'
                mode = varargin{i+1};
            case 'name'
                name = varargin{i+1};
            case 'type'
                type = varargin{i+1};
            otherwise
                error('unknown parameter %s\n',varargin{i});
        end
    end
    if ~exist('SGD','var')
        s = load('SGD.mat');
        SGD = s.SGD;
    end
    if ~exist('mode','var'), mode = 'partial'; end
    if ~exist('name','var'), name = 'orf'; end
    if ~exist('type','var'), type = 'all'; end
    
    if ~iscell(markers)
        markers = {markers};
    end
    nm = length(markers);
    res = cell(nm,1);
    for i = 1:nm
        [chrm ps pe] = parseMarker(markers{i});
        if strcmp(mode, 'complete')
            %genes that are "completely included" in the region !!!
            imask = (SGD.loc(:,1)==chrm & SGD.loc(:,2)>=ps-window & ...
                SGD.loc(:,3)<=pe+window);
        else
            %genes that are "included" in the region, including "partial gene"
            imask = (SGD.loc(:,1)==chrm & SGD.loc(:,2)<=pe+window & ...
                SGD.loc(:,3)>=ps-window);
        end
        if ~strcmp(type,'all')
            imask = imask & (~cellfun('isempty',strfind(SGD.type,type)));
        end
        if strcmp(name, 'index')
            res{i} = find(imask);
        else
            res{i} = SGD.(name)(imask);
        end
    end
    if nm == 1, res = res{1}; end
end

function [chrm ps pe] = parseMarker(marker)
    [chrm remain] = strtok(marker,'_');
    chrm = str2double(strrep(chrm,'M',''));
    [ps pe] = strtok(remain,'_');
    ps = str2double(ps);
    pe = str2double(strrep(pe,'_',''));
end