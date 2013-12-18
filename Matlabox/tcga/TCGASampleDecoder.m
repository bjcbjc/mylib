classdef tcgaSampleDecoder < handle

 
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
            sampstruct = tcgaSampleDecoder.sampleStruct(s{1});
            for i = 2:n
                s = textscan(sampleid{i}, '%s', 'delimiter', '-');
                sampstruct(i,1) = tcgaSampleDecoder.sampleStruct(s{1});
            end
            
            if length(unique({sampstruct.proj})) ~= 1                 
                sampleCollection = [];
                error('sample process error');
            end
            
            sampleCollection.proj = sampstruct(1).proj;
            sampleCollection.tss = {sampstruct.tss}';
            sampleCollection.participant = {sampstruct.participant}';
            sampleCollection.sample = strcat(sampleCollection.tss, sampleCollection.participant);
            if isfield(sampstruct(1), 'samplevial')
                sampleCollection.samplevial = {sampstruct.samplevial}';
                sampleCollection.samplecode = {sampstruct.samplecode}';
                sampleCollection.sampletype = [sampstruct.sampletype]';
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
        
        function [idx1, idx2] = matchSample(sample1, sample2)
            n1 = length(sample1{1});
            n2 = length(sample2{1});
            [n, mi] = min([n1, n2]);
            if mi == 1
                sample2 = cellfun(@(x) x(1:n), sample2, 'uniformoutput', false);
            else
                sample1 = cellfun(@(x) x(1:n), sample1, 'uniformoutput', false);
            end
            [~, idx1, idx2] = intersect(sample1, sample2);            
        end
        
        function [tf idx] = sampleIsMember(sampleid1, sampleid2)
            n1 = unique(cellfun(@length, sampleid1));
            n2 = unique(cellfun(@length, sampleid2));
            if length(n1) > 1 || length(n2) > 1
                error('lengths of sample IDs are different in one of the array');
            end
            if n1 < n2
                for i = 1:length(sampleid2)
                    sampleid2{i} = sampleid2{1}(1:n1);
                end
            elseif n1 > n2
                for i = 1:length(sampleid1)
                    sampleid1{i} = sampleid1{i}(1:n2);
                end
            end
            [tf idx] = ismember(sampleid1, sampleid2);
        end
    end

    methods (Access = private, Static)
        function sample = sampleStruct(s)
            n = length(s);
            sample.proj = s{1};
            sample.tss = s{2};
            sample.participant = s{3};
            if n >= 4
                sample.samplevial = s{4};                
                c = str2double( s{4}(1:end-1) );
                sample.samplecode = s{4}(1:end-1);
                if c < 10 %tumor
                    sample.sampletype = 2;
                elseif c < 20 %normal
                    sample.sampletype = 1;
                else %control
                    sample.sampletype = 0;
                end
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
    
    