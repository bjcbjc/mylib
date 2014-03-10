classdef RNASTATS < handle
    methods (Static)
        function geneBodyCov = readGeneBodyCoverage(pathPattern, sampleNamePattern)            
            if nargin < 2, sampleNamePattern = '/(Sample_[\w\-\_]+)/'; end
            geneBodyCov.fns = listfilename(pathPattern, true);
            geneBodyCov.sample = RNASTATS.getSampleNameFromPath(geneBodyCov.fns, sampleNamePattern);            
            geneBodyCov.count = NaN(101, length(geneBodyCov.fns));
            na = 0;
            for i = 1:length(geneBodyCov.fns)
                t = parseText(geneBodyCov.fns{i}, 'skip', 3, 'ncol',2, 'nrowname', 0, 'ncolname', 0, 'numeric', true);
                if size(t.text, 1) == 101
                    geneBodyCov.count(:, i) = t.text(:,2);
                else
                    na = na + 1;
                end
            end
            if na > 0
                fprintf('%d samples have no gene body coverage infomation\n',na);
            end
            
            geneBodyCov.distribution = bsxfun(@rdivide, geneBodyCov.count, sum(geneBodyCov.count,1));
            vi = ~all(isnan(geneBodyCov.count), 1);
            if sum(vi) > 20
                d = geneBodyCov.distribution(:, vi)';
                eva = evalclusters(d, 'kmeans', 'silhouette', 'klist', 2:5);
                kidx = kmeans(d, eva.OptimalK, 'replicates', 10, 'emptyaction', 'drop');
                geneBodyCov.kmeans = NaN(length(geneBodyCov.fns),1);
                geneBodyCov.kmeans(vi) = kidx;
            end
            RNASTATS.anyNaN(geneBodyCov.count, 'gene body cov');
        end
        
        function [inferExpStat, cmdStatus] = readInferExp(pathPattern, sampleNamePattern)
            %assume ls and grep have the same order of file names
            if nargin < 2, sampleNamePattern = '(Sample_[\w\-\_]+).'; end            
            inferExpStat.pattern = {'1++,1--,2+-,2-+', '1+-,1-+,2++,2--', 'other combination'};            
            
            npattern = length(inferExpStat.pattern);
            inferExpStat.data = cell(1, npattern);            
            valid = true(1, npattern);            
            for i = 1:npattern
                [cmdStatus, out] = system(sprintf('grep "%s" %s', inferExpStat.pattern{i}, pathPattern));                
                out = textscan(out, '%s', 'delimiter', sprintf('\n'));
                out = out{1};                
                if i == 1
                    inferExpStat.sample = regexp(out, sampleNamePattern, 'tokens', 'once');
                    inferExpStat.sample = cellfun(@(x) x{1}, inferExpStat.sample, 'unif', 0);
                end
                out = regexp(out, ': ([\d\.]+)$', 'tokens', 'once');
                out = cellfun(@(x) x{1}, out, 'unif', 0);
                inferExpStat.data{i} = str2double(out);                
                if all(inferExpStat.data{i}==0)
                    valid(i) = false;
                end
            end
            inferExpStat.pattern(~valid) = [];
            inferExpStat.data(~valid) = [];
            RNASTATS.anyNaN(inferExpStat.data, 'inferexp data');
        end
        
        function [spliceRate, cmdStatus] = readSplicingRate(pathPattern)
            cmdStatus = NaN(1,2);
            [cmdStatus(1), out] = system(sprintf('egrep -vh "Spliced_read" %s | cut -f 4', pathPattern));
            [cmdStatus(2), sample] = system(sprintf('egrep -vh "Spliced_read" %s | cut -f 1', pathPattern));
            out = textscan(out, '%s');
            spliceRate.data = str2double(out{1});
            sample = textscan(sample, '%s');
            spliceRate.sample = sample{1};
            RNASTATS.anyNaN(spliceRate.data, 'splicing rate');
        end
        
        function markDupStat = readMarkDuplicate(pathPattern, sampleNamePattern)
            if nargin < 2, sampleNamePattern = '/(Sample_[\w\-\_]+)/'; end
            markDupStat.fns = listfilename(pathPattern, true);
            markDupStat.sample = RNASTATS.getSampleNameFromPath(markDupStat.fns, sampleNamePattern);                        
            for i = 1:length(markDupStat.fns)
                f = fopen(markDupStat.fns{i}, 'r');
                line = fgetl(f);
                while ischar(line)
                    if ~isempty(strfind(line, sprintf('LIBRARY\tUNPAIRED')))
                        break
                    end
                    line = fgetl(f);
                end
                header = textscan( line, '%s', 'delimiter', sprintf('\t'));
                if i == 1
                    markDupStat.label = header{1}(2:end);
                    nstat = length(markDupStat.label);
                    markDupStat.data = NaN(nstat, length(markDupStat.sample));                    
                else
                    if any(~strcmp(markDupStat.label, header{1}(2:end)) )
                        error('inconsist #stats in %s\n', markDupStat.fns{i});
                    end
                end
                data = textscan(fgetl(f), ['%s ', repmat('%f ', 1, nstat)], 'delimiter', sprintf('\t'));
                if ~strcmp(data{1}, markDupStat.sample{i})
                    warning('markdup: %s -> %s',markDupStat.sample{i}, data{1});
                end
                markDupStat.data(:, i) = cell2mat( data(2:end) );
                fclose(f);
            end
            RNASTATS.anyNaN(markDupStat.data, 'markDup');
        end
        
        function readDistribution = readReadDistribution(pathPattern, sampleNamePattern)
            if nargin < 2, sampleNamePattern = '/(Sample_[\w\-\_]+)/'; end
            readDistribution.fns = listfilename(pathPattern, true);
            readDistribution.sample = RNASTATS.getSampleNameFromPath(readDistribution.fns, sampleNamePattern);                        
            for i = 1:length(readDistribution.fns)
                f = fopen(readDistribution.fns{i}, 'r');
                t = textscan(f, '%s %s %s %s', 'delimiter', ' ', 'multipledelimsasone',1, 'headerlines',4,'commentstyle','=');
                t = horzcat(t{:});
                t(1,:) = strrep(t(1,:), '/', '_');
                if i == 1
                    readDistribution.label = t(2:end,1);
                    nstat = length(readDistribution.label);
                    for j = 2:size(t,2)
                        readDistribution.(t{1,j}) = NaN(nstat, length(readDistribution.sample));
                    end                       
                else
                    if any(~strcmp(readDistribution.label, t(2:end,1)) )
                        error('inconsist #stats in %s\n', readDistribution.fns{i});
                    end
                end
                data = str2double( t(2:end, 2:end) );
                for j = 2:size(t,2)
                    readDistribution.(t{1,j})(:, i) = data(:, j-1);
                end
                fclose(f);
            end
            for j = 2:size(t,2)
                RNASTATS.anyNaN(readDistribution.(t{1,j}), sprintf('read distribution, %s', t{1,j}));
            end
        end
        
        function fcSummary = readFeatureCountSummary(pathPattern, sampleNamePattern)
            if nargin < 2, sampleNamePattern = '/(Sample_[\w\-\_]+)/'; end
            fcSummary.fns = listfilename(pathPattern, true);
            fcSummary.sample = RNASTATS.getSampleNameFromPath(fcSummary.fns, sampleNamePattern);     
            for i = 1:length(fcSummary.fns)
                t = parseText(fcSummary.fns{i}, 'skip', 1, 'nrowname', 1, 'ncolname', 0, 'numeric', true, 'ncol', 2);
                if i == 1
                    fcSummary.label = t.rowname;
                    nstat = length(fcSummary.label);
                    fcSummary.data = NaN(nstat, length(fcSummary.sample));
                else
                    if nnz(~strcmp(fcSummary.label, t.rowname)) > 0
                        error('inconsist FC summary labels, %s', fcSummary.fns{i});
                    end
                end
                fcSummary.data(:,i) = t.text;
            end
            RNASTATS.anyNaN(fcSummary.data, 'featureCountsSummary');
        end
        
        function featureCount = readFeatureCount(pathPattern, sampleNamePattern, multlen)
            if nargin < 2, sampleNamePattern = '/(Sample_[\w\-\_]+)/'; end
            if nargin < 3, multlen = false; end
            featureCount.fns = listfilename(pathPattern, true);
            featureCount.sample = RNASTATS.getSampleNameFromPath(featureCount.fns, sampleNamePattern);     
            nsample = length(featureCount.sample);
            for i = 1:nsample
                data = RNASTATS.readFeatureCountSingleFile(featureCount.fns{i});
                if i == 1
                    featureCount.geneId = data.geneId;
                    featureCount.chrm = data.chrm;
                    featureCount.strand = data.strand;
                    if ~multlen
                        featureCount.featureLength = data.featureLength;
                    else
                        featureCount.featureLength = NaN(length(featureCount.geneId), nsample);
                    end
                    featureCount.count = NaN(length(featureCount.geneId), nsample);
                else
                    if nnz(~strcmp(featureCount.geneId, data.geneId)) > 0
                        error('inconsist FC gene IDs, %s', featureCount.fns{i});
                    end
                    if nnz(featureCount.featureLength ~= data.featureLength) > 0
                        error('feature lengths are different between files, %s; set ''multlen'' true', featureCount.fns{i});
                    end
                end
                if multlen
                    featureCount.feaureLength(:,i) = data.featureLength;
                end
                featureCount.count(:,i) = data.count;
            end
            RNASTATS.anyNaN(featureCount.count, 'featureCounts');
        end
        
        function featureCount = readFeatureCountSingleFile(fn)
            %read data from a single file
            t = parseText(fn, 'skip', 1, 'nrowname',0,'ncolname',1,'numericcol',[6,7]);
            featureCount.geneId = t.text(:,strcmpi(t.colname, 'geneid'));
            featureCount.chrm = regexp(t.text(:,strcmpi(t.colname, 'chr')), '^(\w+)', 'match', 'once');
            featureCount.strand = regexp(t.text(:,strcmpi(t.colname, 'strand')), '([\+\-])', 'match', 'once');
            featureCount.featureLength = t.numtext(:, strcmpi(t.numcolname, 'length'));
            featureCount.count = t.numtext(:, ~strcmpi(t.numcolname, 'length'));
        end
        
        function GC = readGCDistribution(fn)
            %fn is the matrix text format from Heather's R data
            %GC_read_distribution.txt
            t = parseText(fn, 'nrowname', 1, 'ncolname', 1, 'numeric', true);
            GC.sample = t.colname;
            GC.GCContent = str2double(t.rowname);
            GC.data = t.text;
            GC.note = sprintf('data is read count per million reads in each sample');
            RNASTATS.anyNaN(GC.data);
        end
        
        function innerDistance = readInnerDistanceDistribution(fn, distanceBin)
            %fn is the matrix text format from Heather's R data
            %inner_distance_read_distribution.txt
            %setting from RSEQC: -250:5:250 for bins of inner distance
            if nargin < 2, distanceBin = -250:5:250; end
            t = parseText(fn, 'nrowname', 1, 'ncolname', 1, 'numeric', true);            
            innerDistance.sample = t.colname;
            innerDistance.distance = (distanceBin(1:end-1)+distanceBin(2:end))./2;
            innerDistance.data = t.text;
            innerDistance.note = sprintf('data is read count per million reads in each sample');
            RNASTATS.anyNaN(innerDistance.data);
        end
        
        function report = readReportStats(reportpath, varargin)
            para.fns = {'QC_statistics.txt', ...                
                'Read_distribution_statistics.txt'};
                        
            para = assignpara(para, varargin{:});
            if reportpath(end) ~= '/', reportpath = [reportpath '/']; end
            
            for i = 1:length(para.fns)
                t = parseText([reportpath para.fns{i}], 'nrowname', 1, 'ncolname', 1, 'numeric', true);
                if i == 1
                    report.sample = t.rowname;
                    nlabel = length(t.colname);
                    report.label = cell(nlabel+50, 1);
                    report.data = NaN(nlabel+50, length(report.sample));
                    report.label(1:nlabel) = strrep(t.colname, '(%)', '_PCT');
                    report.data(1:nlabel, :) = t.text';
                else
                    if nnz(~strcmp(t.rowname, report.sample)) > 0
                        error('sample names are not the same, %s', para.fns{i});
                    end
                    addnlabel = length(t.colname);
                    report.label(nlabel+1:nlabel+addnlabel) = strrep(t.colname, '(%)', '_PCT');
                    report.data(nlabel+1:nlabel+addnlabel, :) = t.text';
                    nlabel = nlabel + addnlabel;
                end
            end
            report.label(nlabel+1:end) = [];
            report.data(nlabel+1:end, :) = [];
            RNASTATS.anyNaN(report.data, 'readReportStats');
        end
        
        function readcount = readReportReadCount(reportpath, varargin)
            para.fns = {'DESeq_normalized_count_matrix.txt', ...
                'featureCounts_count_matrix.txt'};
            para.stripstr = '_count_matrix.txt';
            
            para = assignpara(para, varargin{:});
            if reportpath(end) ~= '/', reportpath = [reportpath '/']; end
            
            fdname = strrep(para.fns, para.stripstr, '');
            for i = 1:length(para.fns)
                if ~isempty(strfind(lower(para.fns{i}), 'featurecounts'))
                    t = parseText([reportpath para.fns{i}], 'nrowname', 2, 'ncolname', 1, 'numeric', true);
                    rmi = strcmpi(t.colname, 'Length');
                    readcount.geneLength = t.text(:, rmi);
                    t.colname(rmi) = [];
                    t.text(:, rmi) = [];
                    rmi = strcmpi(t.rownamelabel, 'Chromosome');
                    t.rownamelabel(rmi) = [];
                    t.rowname(:, rmi) = [];
                else
                    t = parseText([reportpath para.fns{i}], 'nrowname', 1, 'ncolname', 1, 'numeric', true);
                end
                if i == 1
                    readcount.sample = t.colname;                    
                    readcount.geneId = t.rowname;                    
                else
                    if nnz(~strcmp(t.colname, readcount.sample)) > 0
                        error('sample names are not the same, %s', para.fns{i});
                    end
                    if nnz(~strcmp(t.rowname, readcount.geneId)) > 0
                        error('gene names are not the same, %s', para.fns{i});
                    end
                    
                end
                readcount.(fdname{i}) = t.text;
                RNASTATS.anyNaN(readcount.(fdname{i}), sprintf('readReportReadCount, %s', fdname{i}));
            end            
        end
        
        function alleleCount = readAlleleCount(pathPattern, sampleNamePattern, locfile, stranded)
            alleleCount.fns = listfilename(pathPattern, true);
            alleleCount.sample = RNASTATS.getSampleNameFromPath(alleleCount.fns, sampleNamePattern);            
            
            t = parseText(locfile, 'ncolname',0,'nrowname',1,'numeric',true,'ncol',2);
            alleleCount.locidx = gloc2index( numericchrm( t.rowname), t.text);
            
            for i = 1:length(alleleCount.fns)
                data = AlleleCountData.readTableFormatOutput(alleleCount.fns{i}, stranded);
                if i == 1
                    alleleCount.ntBase = data.ntbase;
                    alleleCount.count = NaN(length(alleleCount.locidx), length(data.ntbase), length(alleleCount.fns));
                    alleleCount.indel = cell(length(alleleCount.locidx), length(alleleCount.fns));
                    alleleCount.indelCount = cell(length(alleleCount.locidx), length(alleleCount.fns));
                end
                [~, idx] = ismember(data.locidx, alleleCount.locidx);
                alleleCount.count(idx, :, i) = data.count;
                alleleCount.indel(idx,i) = data.indel;
                alleleCount.indelCount(idx,i) = data.indelcount;
            end
        end
        
        function sample = getSampleNameFromPath(filepaths, samplepattern)
            sample = regexp(filepaths, samplepattern, 'tokens', 'once');
            sample = cellfun(@(x) x{1}, sample, 'unif', 0);
        end
                                
    end
    methods (Access=private, Static)
        function anyNaN(data, msg)
            if any(isnan(data(:)))
                warning('NaN data, %s', msg);
            end
        end
    end
end