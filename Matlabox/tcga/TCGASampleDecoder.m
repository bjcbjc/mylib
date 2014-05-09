classdef TCGASampleDecoder < handle

    properties (Constant)
        indextb = {'tss', 7; 'participant', 12; 'patient', 12; ...
            'sample', 15; 'vial', 16; 'portion', 20; 'analyte', 21; ...
            'plate', 26; 'center', 29};
        sampleTypeTb = { ...
            1, 'Primary solid Tumor', 'TP'; ...
            2, 'Recurrent Solid Tumor', 'TR'; ...
            3, 'Primary Blood Derived Cancer - Peripheral Blood', 'TB'; ...
            4, 'Recurrent Blood Derived Cancer - Bone Marrow', 'TRBM'; ...
            5, 'Additional - New Primary', 'TAP'; ...
            6, 'Metastatic', 'TM'; ...
            7, 'Additional Metastatic', 'TAM'; ...
            8, 'Human Tumor Original Cells', 'THOC'; ...
            9, 'Primary Blood Derived Cancer - Bone Marrow', 'TBM'; ...
            10, 'Blood Derived Normal', 'NB'; ...
            11, 'Solid Tissue Normal', 'NT'; ...
            12, 'Buccal Cell Normal', 'NBC'; ...
            13, 'EBV Immortalized Normal', 'NEBV'; ...
            14, 'Bone Marrow Normal', 'NBM'; ...
            20, 'Control Analyte', 'CELLC'; ...
            40, 'Recurrent Blood Derived Cancer - Peripheral Blood', 'TRB'; ...
            50, 'Cell Lines', 'CELL'; ...
            60, 'Primary Xenograft Tissue', 'XP'; ...
            61, 'Cell Line Derived Xenograft Tissue', 'XCL'};
    end
 
    methods (Static)
        
        function sampleCollection = decode(sampleid, fulldecode)
            if nargin < 2, fulldecode = false; end
            
            if ~iscellstr(sampleid) && ischar(sampleid)
                sampleid = {sampleid};
            elseif ~iscellstr(sampleid)
                error('sample id must be arrays of strings or a string');
            end
            
            %format
            %TCGA-TSS-participant-XXY-WWZ-plate-center
            %TSS: tissue source site
            %XX: sample, normal: 10-19, tumor:01-09, control:20-29
            %Y: vial
            %WW: portion
            %Z: analyte
            %
            
            n = length(sampleid);
            s = textscan(sampleid{1}, '%s', 'delimiter', '-');
            sampstruct = TCGASampleDecoder.sampleStruct(s{1});
            for i = 2:n
                s = textscan(sampleid{i}, '%s', 'delimiter', '-');
                sampstruct(i,1) = TCGASampleDecoder.sampleStruct(s{1});
            end
            
            if length(unique({sampstruct.proj})) ~= 1                 
                sampleCollection = [];
                error('sample process error');
            end
            
            sampleCollection.proj = sampstruct(1).proj;
            sampleCollection.tss = {sampstruct.tss}';
            sampleCollection.participant = {sampstruct.participant}';
%             sampleCollection.sample = strcat(sampleCollection.tss, sampleCollection.participant);
            if isfield(sampstruct(1), 'vial')
                sampleCollection.vial = {sampstruct.vial}';
                sampleCollection.sample = [sampstruct.sample]';
                %sampleCollection.samplecode = {sampstruct.samplecode}';
                %sampleCollection.sampletype = [sampstruct.sampletype]';
            end
            if fulldecode
                if isfield(sampstruct(1), 'portion')
                    sampleCollection.portion = {sampstruct.portion}';
                end
                if isfield(sampstruct(1), 'plate')
                    sampleCollection.plate = {sampstruct.plate}';
                end
                if isfield(sampstruct(1), 'center')
                    sampleCollection.center = {sampstruct.center}';
                end
            end
        end
        
        function [sample, idx1, idx2] = intersect(sample1, sample2, reduceto)
            if nargin < 3, reduceto = []; end
            if isempty(reduceto)
                n1 = length(sample1{1});
                n2 = length(sample2{1});
                [n, mi] = min([n1, n2]);
                if mi == 1
                    sample2 = cellfun(@(x) x(1:n), sample2, 'uniformoutput', false);
                else
                    sample1 = cellfun(@(x) x(1:n), sample1, 'uniformoutput', false);
                end
            else
                sample1 = TCGASampleDecoder.reduceBarcode(sample1, reduceto);
                sample2 = TCGASampleDecoder.reduceBarcode(sample2, reduceto);
            end
            [sample, idx1, idx2] = intersect(sample1, sample2);            
        end
        
        function patientcode = reduceToParticipant(barcode)
            % barcode: TCGA barcodes
            % strip barcodes to participant level
            patientcode = TCGASampleDecoder.reduceBarcode(barcode, 'patient');            
        end
        
        function newcode = reduceBarcode(barcode, labelIndex)
            % Reduce barcode to indicated index (numeric index value or
            % TCGA barcode labels eg. Participant, Sample, Vial, Portion,
            % Analyte and etc. Return bardcode(1:k), k is labelIndex or the
            % index defined by the label
            %
            % labelIndex: index number or label (string), valid strings:
            % {'TSS', 'Participant'(or 'patient'), 'Sample', 'Vial',
            % 'Analyte', 'Plate', 'Center'}, CASE-INSENSITIVE
            %
            
            newcode = barcode;
            if ischar(labelIndex)
                [~, j] = ismember(lower(labelIndex), TCGASampleDecoder.indextb(:,1));
                if j == 0
                    fprintf('Cannot find label %s; return input barcodes\n', labelIndex);                    
                    return
                end
                labelIndex = TCGASampleDecoder.indextb{j, 2};
            end
            n = cellfun(@length, barcode);
            idx = n > labelIndex;
            if any( n < labelIndex)
                fprintf('Some barcodes are shorted than the requested length\n');
            end
            newcode(idx) = cellfun(@(x) x(1:labelIndex), barcode(idx), 'unif', 0);
        end
        
        function [matchcode, idx] = longestMatch(barcodelist, query)
            % return the longest matched barcode of query, compared to
            % barcodelist; query is a string
            %
            % barcodelist: list of barcodes to search in
            % query: a barcode for matching
            %
            % matchcode: the longest barcode that matches the prefix of
            %   query
            % idx: indices to barcodelist that contains the longest match
            %
            %
            indices = cell2mat(TCGASampleDecoder.indextb(:,2));
            matchcode = '';
            n = length(query);
            indices( indices > n ) = [];
            for i = length(indices):-1:1
                idx = find( cellfun(@(x) ismember(1,x), strfind(barcodelist, query(1:indices(i))) ) ); 
                if ~isempty(idx)
                    matchcode = query(1:indices(i));
                    return
                end
            end            
        end
        
        function [tf, idx] = ismember(sampleid1, sampleid2, reduceto)
            if nargin < 3, reduceto = []; end
            if ischar(sampleid1), sampleid1 = {sampleid1}; end
            if ischar(sampleid2), sampleid2 = {sampleid2}; end
            if isempty(reduceto)
                n1 = unique(cellfun(@length, sampleid1));
                n2 = unique(cellfun(@length, sampleid2));
                if length(n1) > 1 || length(n2) > 1
                    error('lengths of sample IDs are different in one of the array');
                end
                if n1 < n2
                    sampleid2 = TCGASampleDecoder.reduceBarcode(sampleid2, n1);
                elseif n1 > n2
                    sampleid1 = TCGASampleDecoder.reduceBarcode(sampleid1, n2);
                end
            else
                sampleid1 = TCGASampleDecoder.reduceBarcode(sampleid1, reduceto);
                sampleid2 = TCGASampleDecoder.reduceBarcode(sampleid2, reduceto);
            end
            [tf, idx] = ismember(sampleid1, sampleid2);
        end
        
        function sampletype = translateSampleType(sampleTypeCode, shortname)
            if nargin < 2, shortname = false; end
            if ischar(sampleTypeCode), sampleTypeCode = {sampleTypeCode}; end
            if iscellstr(sampleTypeCode)
                s = TCGASampleDecoder.decode(sampleTypeCode);
                sampleTypeCode = s.sample;
            end
            code = cell2mat(TCGASampleDecoder.sampleTypeTb(:,1));
            [~, i] = ismember(sampleTypeCode, code);
            sampletype = repmat({'NA'}, length(sampleTypeCode), 1);            
            sampletype(i~=0) = TCGASampleDecoder.sampleTypeTb(i(i~=0), shortname+2);
            if length(sampletype) == 1
                sampletype = sampletype{1};
            end
        end
        
        function pair = getSamplePair(barcode, refType)
            if nargin < 2, refType = 11; end
            sample = TCGASampleDecoder.decode(barcode);
            refIdx = sample.sample == refType;
            valid = ismember(sample.participant, sample.participant(refIdx));
            [count, upatient] = eleCounts(sample.participant(valid));
            validPatient = upatient(count > 1); %has pair
            nonRefIdx = find( ismember(sample.participant, validPatient) & ...
                sample.sample ~= refType);
            n = length(nonRefIdx);
            pair = cell(n ,2);
            pair(:,2) = barcode(nonRefIdx);
            refSample = barcode(refIdx);
            curLast = n;
            for i = 1:n
                [~, idx] = TCGASampleDecoder.longestMatch(refSample, pair{i,2});
                nMatch = length(idx);
                if nMatch == 1
                    pair{i,1} = refSample{idx};
                elseif nMatch > 1
                    pair{i,1} = refSample{idx};
                    pair(curLast+1:curLast+nMatch-1, 2) = pair(i,2);
                    pair(curLast+1:curLast+nMatch-1, 1) = refSample(idx);
                    curLast = curLast + nMatch - 1;
                end
            end
            rmIdx = cellfun(@isempty, pair(:,1));
            pair(rmIdx, :) = [];
        end
        
    end

    methods (Access = private, Static)
        function sample = sampleStruct(s)
            n = length(s);
            sample.proj = s{1};
            sample.tss = s{2};
            sample.participant = s{3};
            if n >= 4
                sample.vial = s{4}(end);
                sample.sample = str2double(s{4}(1:end-1));
%                 sample.samplevial = s{4};                
%                 c = str2double( s{4}(1:end-1) );
%                 sample.samplecode = s{4}(1:end-1);
%                 if c < 10 %tumor
%                     sample.sampletype = 2;
%                 elseif c < 20 %normal
%                     sample.sampletype = 1;
%                 else %control
%                     sample.sampletype = 0;
%                 end
            end
            if n >= 5
                sample.portion = s{5};
            end
            if n >= 6
                sample.plate = s{6};
            end
            if n >= 7
                sample.center = s{7};
            end
        end
    end
end
    
    