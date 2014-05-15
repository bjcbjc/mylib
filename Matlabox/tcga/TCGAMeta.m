classdef TCGAMeta < handle
    
    methods (Static)
        function data = readSDRF(fn, sample)
            if nargin < 2, sample = {}; end
                       
            keepProtocol = {'reverse_transcription', 'library_preparation', 'RNA_Sequencing'};
            t = parseText(fn, 'ncolname', 1, 'nrowname', 0);
            barcodeIdx = strcmp(t.colname, 'Comment [TCGA Barcode]');
            if ~isempty(sample)
                [~, si] = ismember(sample, t.text(:, barcodeIdx));
                data.sample = sample;
            else
                [data.sample, si] = unique( t.text(:, barcodeIdx));
            end
            data.material = t.text(si, strcmp(t.colname, 'Material Type'));
            n = length(unique(data.material));
            if n > 1
                fprintf('Material Type has %d unique values\n',n);
            end
            idx = find(strcmp(t.colname, 'Protocol REF'));
            for i = 1:length(keepProtocol)
                colIdx = ~cellfun(@isempty, strfind(t.text(1,idx), keepProtocol{i}));
                if any(colIdx)
                    data.(keepProtocol{i}) = t.text(si, idx(colIdx));
                    n = length(unique(data.(keepProtocol{i})));
                    if n > 1
                        fprintf('%s has %d unique values\n',keepProtocol{i},n);
                    end
                end
            end
            idx = strcmp(t.colname, 'Comment [TCGA Include for Analysis]');
            data.includeInAnalysis = true(size(data.sample));
            for i = 1:length(data.sample)
                if any(~strcmpi(t.text(strcmp(t.text(:,barcodeIdx), data.sample{i}), idx), 'yes'))
                    data.includeInAnalysis(i) = false;
                end
            end
            if ~all(data.includeInAnalysis)
                fprintf('%d samples are excluded from TCGA analysis\n',sum(data.includeInAnalysis));
            end
        end
    end
end