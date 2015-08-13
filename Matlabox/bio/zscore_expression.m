function zscore_expression(cohortFn, countFn, outPrefix, cohortName, validSampleFn)
if nargin < 5, validSampleFn = ''; end

%for gene name mapping
GENCODE = loadStructData('/nethome/bjchen/DATA/GenomicInfo/GENCODE/GENCODE.gene.v18.mat');

%read cohort count data
cohort = parseText(cohortFn, 'nrowname', 2, 'ncolname', 1, 'numeric', true);
rmi = strcmpi(cohort.colname, 'Length');
cohort.geneLength = cohort.text(:, rmi);
cohort.colname(rmi) = [];
cohort.text(:, rmi) = [];
rmi = strcmpi(cohort.rownamelabel, 'Chromosome') | strcmpi(cohort.rownamelabel, 'Chrm');
cohort.rownamelabel(rmi) = [];
cohort.rowname(:, rmi) = [];
cohort.count = cohort.text;
cohort.geneId = cohort.rowname;
cohort.sample = cohort.colname;
cohort = rmfield(cohort, {'rowname', 'colname', 'text', 'rownamelabel'});

%filter samples from cohort?
if ~isempty(validSampleFn)
    f = fopen(validSampleFn);
    validSample = textscan(f, '%s');
    validSample = validSample{1};
    fclose(f);
    cohort = DATASTRUCTFUNC.orderDataStruct(cohort, 'sample', validSample);%save memory
end

cohort.DESeq_normalized = RNACOUNT.getDESeqNormalizedCount(cohort.count);
cohort.logread = log2( cohort.DESeq_normalized + 1);
cohort.logTPM = log2( RNACOUNT.getTPM( cohort.DESeq_normalized, cohort.geneLength, 0.9));


%read sample count data
ReadCount = RNASTATS.readFeatureCountSingleFile(countFn);
ReadCount = DATASTRUCTFUNC.orderDataStruct(ReadCount, 'geneId', cohort.geneId);
[~, idx] = ismember(ReadCount.geneId, GENCODE.id);
ReadCount.geneName = GENCODE.name(idx);


%normalization
ReadCount.DESeq_normalized = RNACOUNT.getDESeqNormalizedCount( ...
    cohort.count, ReadCount.count);
ReadCount.TPM = RNACOUNT.getTPM(ReadCount.DESeq_normalized, ReadCount.featureLength, 0);
ReadCount.logTPM = log2( RNACOUNT.getTPM( ...
    ReadCount.DESeq_normalized, ReadCount.featureLength, 0.9));
ReadCount.logread = log2( ReadCount.DESeq_normalized + 1);
ReadCount.zBy_count = modzscore(ReadCount.logread, cohort.logread, 2);
ReadCount.zBy_TPM = modzscore(ReadCount.logTPM, cohort.logTPM, 2);
ReadCount.logRatio_count = ReadCount.logread - nanmedian(cohort.logread,2);
ReadCount.logRatio_TPM = ReadCount.logTPM - nanmedian(cohort.logTPM,2);


sample = strrep(regexp(countFn, 'Sample_[\w\-\_]+', 'match', 'once'), 'Sample_', '');
% tag = {'count', 'TPM', ...
%     'zBy_count','zBy_TPM', ...
%     'logRatio_count', 'logRatio_TPM', ...
%     };
tag = {'count', ...
    'zBy_count', ...
    'logRatio_count', ...
    };
header = {'gene', 'gencodeId', sample};    
for i = 1:length(tag)
    if ~isempty(strfind(tag{i}, 'zBy'))
        outFn = [outPrefix, '.', strrep(tag{i}, '_', sprintf('%s_', cohortName)), '.txt'];
    else
        outFn = [outPrefix, '.', strrep(tag{i}, 'Ratio', sprintf('Ratio2%s', cohortName)), '.txt'];
    end
    table = cell(length(ReadCount.geneId), 3);
    table(:,1) = ReadCount.geneName;
    table(:,2) = ReadCount.geneId;
    table(:,3) = numarray2strarray(ReadCount.(tag{i}));    
    tabwrite( outFn, [header; table]);
end

