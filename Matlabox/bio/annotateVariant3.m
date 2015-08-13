function annotateVariant3(variantFn, cohortReadCountFn, sampleReadCountFn, rnaAlleleCountFn, outputFn, varargin)

    para.stranded = false;
    para.locColumnName = {'Chromosome', 'Start_Position', 'End_Position'};
    para.refColumnName = 'Reference_Allele';
    para.altColumnName = 'Tumor_Seq_Allele2';
    para.refCountColumnName = 't_ref_count';
    para.altCountColumnName = 't_alt_count';

    %ASE with purity
    para.minAltReadForASE = 6;
    para.minTotalReadForASE = 6;
    para.minAEFForASE = 0.02;
    para.minNumVarForASE = 30;
    para.runTestASE = true;
    
    para.GTEx = '/nethome/bjchen/DATA/GTEx/GTEx_Gene.mat'; %can pre-load to speed up multiple samples
    para.GTExSample = '';
%     para.GTExTissueKey = 'brain';
    para.cohortSubtypeInfoFile = '';
    para.rnaGtfFn = '/nethome/bjchen/DATA/GenomicInfo/GENCODE/gencode.v18.annotation.gene.gtf';
    para.gtexGtfFn = '/nethome/bjchen/DATA/GTEx/20140117/gtex.gencode.v18.genes.removedNotInGTEx.gtf';
    
    para = assignpara(para, varargin{:});    
    if para.stranded
        error('stranded not supported yet');
    end

    if isempty(para.GTExSample)
        error('GTExSample is empty; need list of GTEx samples');
    end
    f = fopen(para.GTExSample);
    para.GTExSample = textscan(f, '%s');
    para.GTExSample = para.GTExSample{1};
    fclose(f);
    
    %read vaiant
    t = parseText(variantFn, 'ncolname', 1, 'nrowname', 0, 'numeric', false, 'commentstyle','#');
    [~, colidx] = ismember(para.locColumnName, t.colname);
    if any(colidx==0)
        error('cannot find columns for SNV locations');
    end
    variant.loc = [numericchrm( t.text(:, colidx(1))), str2double(t.text(:, colidx(2:3)))];
    variant.locidx = gloc2index(variant.loc(:,1:2));
    [~, colidx] = ismember({para.refColumnName, para.altColumnName}, t.colname);
    variant.ref = t.text(:, colidx(1));
    variant.alt = t.text(:, colidx(2));
    [~, colidx] = ismember({para.refCountColumnName, para.altCountColumnName}, t.colname);
    variant.dnaCount = str2double(t.text(:, colidx));
    
    %sample readcount
    sample = RNASTATS.readFeatureCountSingleFile(sampleReadCountFn);
    
    %get expression annotation
    [expDataHeader, expData] = getExpression(cohortReadCountFn, sample, variant, para);    
    %read RNA allele count
    [acDataHeader, acData] = getRnaAlleleCount(rnaAlleleCountFn, variant, para);
    %test ASE
    [aseDataHeader, aseData] = runTestASE(acData, variant.dnaCount, para);
    %convert to string
    acData = numarray2strarray(acData);
    %get GTEx expression
    [gtexDataHeader, gtexData] = getGtexData(sample, variant, para);
    
    %output
    f = fopen(variantFn, 'r');
    fout = fopen(outputFn, 'w');
    fpos = ftell(f);
    line = fgetl(f);
    offset = fpos - ftell(f);
    while ~isempty( regexp(line, '^#', 'once'))        
        fprintf(fout, '%s\n',line);
        fpos = ftell(f);
        line = fgetl(f);
        offset = fpos - ftell(f);                
    end
    fseek(f, offset, 'cof');        
    
    if isempty(gtexData)
        fprintf(fout, '%s\t%s\t%s\t%s\n', fgetl(f), strjoin(expDataHeader, sprintf('\t')), ...
            strjoin(acDataHeader, sprintf('\t')), strjoin(aseDataHeader, sprintf('\t')) );
        line = fgetl(f);
        lineIdx = 1;
        while ischar(line)
            fprintf(fout, '%s\t%s\t%s\t%s\n', line, strjoin(expData(lineIdx,:), sprintf('\t')), ...
                strjoin(acData(lineIdx,:), sprintf('\t')), strjoin(aseData(lineIdx,:), sprintf('\t')) );
            line = fgetl(f);
            lineIdx = lineIdx + 1;
        end
    else
        fprintf(fout, '%s\t%s\t%s\t%s\t%s\n', fgetl(f), strjoin(expDataHeader, sprintf('\t')), ...
            strjoin(acDataHeader, sprintf('\t')), strjoin(aseDataHeader, sprintf('\t')), ...
            strjoin(gtexDataHeader, sprintf('\t')) );
        line = fgetl(f);
        lineIdx = 1;
        while ischar(line)
            fprintf(fout, '%s\t%s\t%s\t%s\t%s\n', line, strjoin(expData(lineIdx,:), sprintf('\t')), ...
                strjoin(acData(lineIdx,:), sprintf('\t')), strjoin(aseData(lineIdx,:), sprintf('\t')), ...
                strjoin(gtexData(lineIdx,:), sprintf('\t')) );
            line = fgetl(f);
            lineIdx = lineIdx + 1;
        end
    end
    fclose(f);
    fclose(fout);
end


function [expDataHeader, expData] = getExpression(cohortReadCountFn, sample, variant, para)
    if ~isempty(para.cohortSubtypeInfoFile)
        cohortSubtypeInfo = readSubtypeInfo(para.cohortSubtypeInfoFile);
    else
        cohortSubtypeInfo = [];
    end
    
    %read expression
    if ischar(cohortReadCountFn)
        cohort = loadStructData(cohortReadCountFn);
    else
        cohort = cohortReadCountFn;
    end
    if isfield(cohort, 'highRRnaSample')
        validSampleForZ = ~ismember(cohort.sample, cohort.highRRnaSample);
%         rmi = ismember(cohort.sample, cohort.highRRnaSample);
%         cohort.sample(rmi) = [];
%         cohort.DESeq_normalized(:, rmi) = [];
%         cohort.featureCounts(:, rmi) = [];
    else
        validSampleForZ  = true(size(cohort.sample));
    end
    cohort.logread = log2( cohort.DESeq_normalized(:, validSampleForZ) + 1);
    cohort.logmedian = nanmedian( cohort.logread, 2);

%     cohort.logread = log2( cohort.DESeq_normalized + 1);
%     cohort.logmedian = nanmedian( cohort.logread, 2);    
    
    sample.geneName = GENOMEFUNC.idQuery(cohort, sample.geneId, 'geneName', 'geneId');
    if ~strcmp(sample.geneId, cohort.geneId)
        error('inconsistent geneId between sampleReadCount and cohortReadCount');
    end
    %calculate expression: 1. raw read counts for genes, 2. TPM, 3. log2(normalized
    % count+1) - median(log2(normalized cohort count+1)), 4. modified z score (x-med)/MAD
    normcount = RNACOUNT.getDESeqNormalizedCount(cohort.featureCounts, sample.count);
    sample.tpm = RNACOUNT.getTPM(normcount, sample.featureLength, 0);
%     sample.tpm = RNACOUNT.getTPM(normcount, sample.featureLength, 0.9);
    sample.lognormcount = log2(normcount+1);
    sample.logRatio2Cohort = sample.lognormcount - cohort.logmedian;
    sample.z = modzscore( sample.lognormcount, cohort.logread, 2);
    clear t

%     expDataHeader = {'gene_name_expression', 'raw_read_count', ...
%         'transcript_per_million', 'log2_ratio_to_cohort', 'zscore_cohort'};
    expDataHeader = {'gene_name_expression', 'raw_read_count', ...
        'transcript_per_million', 'log2_ratio_to_TCGA', 'zscore_TCGA'};
    if ~isempty(cohortSubtypeInfo) %prepare subtype data
        subtypes = fieldnames(cohortSubtypeInfo);
        valid = false(1, length(subtypes));
        for i = 1:length(subtypes)
            if length(cohortSubtypeInfo.(subtypes{i})) > 10
                valid(i) = true;
                expDataHeader{end+1} = sprintf('log2_ratio_to_%s',subtypes{i});
                expDataHeader{end+1} = sprintf('zscore_%s',subtypes{i});
                validCohortSampIdx = ismember(cohort.sample, cohortSubtypeInfo.(subtypes{i}));
                sample.(['log' subtypes{i}]) = sample.lognormcount - nanmedian(cohort.logread(:,validCohortSampIdx), 2);
                sample.(['z' subtypes{i}]) = modzscore(sample.lognormcount, cohort.logread(:,validCohortSampIdx), 2);
            end
        end
        subtypes = subtypes(valid);
    else
        subtypes = {};
    end

    %find the gene(s) first, need for loop in case more than one gene
    nvar = length(variant.locidx);
    expData = cell(nvar, length(expDataHeader)); %first column for gene name
    expData(:) = {'NA'};
        
    %use bedtools to query overlap genes     
    matchGene = GENOMEFUNC.searchAllOverlapGTF( variant.loc, para.rnaGtfFn, ...
        'padChr', true, 'mt', 'M', 'geneIdFormat', 'ENSG[\w\.]+');
    [matchCount, uniqMatch] = eleCounts(matchGene.locidx_start);
    uniqIdx = ismember(variant.locidx, uniqMatch(matchCount == 1));
    [~, uniqMatchIdx] = ismember(variant.locidx(uniqIdx), matchGene.locidx_start);
    multiIdx = find( ismember(variant.locidx, uniqMatch(matchCount > 1)) );    
    [~, gencodeIdx2Sample] = ismember(matchGene.geneId, sample.geneId);
    
    %variant mapped to only one gene
    gidx = gencodeIdx2Sample(uniqMatchIdx);
    expData(uniqIdx,1) = sample.geneName(gidx);
    expData(uniqIdx,2) = numarray2strarray(sample.count(gidx));
    expData(uniqIdx,3) = numarray2strarray(sample.tpm(gidx));
    expData(uniqIdx,4) = numarray2strarray(sample.logRatio2Cohort(gidx));
    expData(uniqIdx,5) = numarray2strarray(sample.z(gidx));
    for subtypeIdx = 1:length(subtypes)
        expData(uniqIdx, 5+2*subtypeIdx-1) = numarray2strarray(sample.(['log' subtypes{subtypeIdx}])(gidx));
        expData(uniqIdx, 5+2*subtypeIdx) = numarray2strarray(sample.(['z' subtypes{subtypeIdx}])(gidx));
    end
    
    %variant mapped to multiple genes
    for i = 1:length(multiIdx)
        gidx = gencodeIdx2Sample( matchGene.locidx_start == variant.locidx(multiIdx(i)) );        
        gidx(gidx==0) = [];
        
        if ~isempty(gidx)
            expData{multiIdx(i), 1} = strjoin(sample.geneName(gidx)', ', ');
            expData{multiIdx(i), 2} = strjoin(numarray2strarray(sample.count(gidx))', ', ');
            expData{multiIdx(i), 3} = strjoin(numarray2strarray(sample.tpm(gidx))', ', ');
            expData{multiIdx(i), 4} = strjoin(numarray2strarray(sample.logRatio2Cohort(gidx))', ', ');
            expData{multiIdx(i), 5} = strjoin(numarray2strarray(sample.z(gidx))', ', ');
            for subtypeIdx = 1:length(subtypes)
                expData{multiIdx(i), 5+2*subtypeIdx-1} = strjoin(numarray2strarray(sample.(['log' subtypes{subtypeIdx}])(gidx))', ', ');
                expData{multiIdx(i), 5+2*subtypeIdx} = strjoin(numarray2strarray(sample.(['z' subtypes{subtypeIdx}])(gidx))', ', ');
            end        
        end
    end
    expData(:,2:end) = strrep(expData(:, 2:end), 'NaN', 'NA');
end

function [acDataHeader, acData] = getRnaAlleleCount(rnaAlleleCountFn, variant, para)
    alleleCount = AlleleCountData.readTableFormatOutput(rnaAlleleCountFn, para.stranded);    
    nread = sum(alleleCount.count(:, ismember(alleleCount.ntbase, {'A', 'T','C','G','N'})), 2);
    nskip = sum(alleleCount.count(:, ismember(alleleCount.ntbase, {'>', '*'})), 2);
    ntotal = nread + nskip;
    mtxsize = size(alleleCount.count);
    [~, refIdx] = ismember(variant.ref, alleleCount.ntbase);
    [~, altIdx] = ismember(variant.alt, alleleCount.ntbase);        
    [~, varIdx] = ismember(variant.locidx, alleleCount.locidx);
    valid = refIdx ~= 0 & altIdx ~= 0 & varIdx ~= 0;
    
    nvar = length(variant.locidx);
    %get RNA allele counts    
    acData = zeros( nvar, 4); %ref, alt, skip, other
    acDataHeader = {'ref_RNA_count', 'alt_RNA_count', 'skip_RNA_count', 'other_RNA_count'};
    acData(valid,1:3) = [ alleleCount.count( sub2ind(mtxsize, varIdx(valid), refIdx(valid)) ), ...
        alleleCount.count( sub2ind(mtxsize, varIdx(valid), altIdx(valid)) ), ...
        nskip(varIdx(valid))];
    acData(valid, 4) = ntotal(varIdx(valid)) - sum(acData(valid, 1:3),2);
    acData(refIdx == 0 | altIdx == 0,:) = NaN;    
end

function [gtexDataHeader, gtexData] = getGtexData(sample, variant, para)
    if ~isempty(para.GTEx)
        nvar = length(variant.locidx);
        if ischar(para.GTEx)
            GTEx = loadStructData(para.GTEx);
        else
            GTEx = para.GTEx;
        end

        %limit to subset of samples (eg. samples that are highly correlated
        %with TCGA cohort)
        validSample = ismember(GTEx.sample, para.GTExSample);
%         validSample = ~cellfun(@isempty, strfind(lower(GTEx.SMTSD), lower(para.GTExTissueKey)));
        GTEx.pctSampleExpressed = NaN(length(GTEx.gene), 3);

        [~, sampGeneIdx, gtexGeneIdx] = GENOMEFUNC.intersectGencode(sample.geneId, GTEx.gencodeId);
        gtexLogRead = log2( RNACOUNT.getDESeqNormalizedCount( GTEx.readcount(gtexGeneIdx, validSample))+1);
        gtexLogMedian = nanmedian(gtexLogRead , 2);
        
        sampleLogCount = log2(RNACOUNT.getDESeqNormalizedCount(GTEx.readcount(gtexGeneIdx,validSample), sample.count(sampGeneIdx,:))+1);
        sample.logRatio2GTEx = NaN(length(sample.geneId),2); %logratio and median of log readcount of GTEx samples
        sample.logRatio2GTEx(sampGeneIdx,:) = [sampleLogCount - gtexLogMedian, gtexLogMedian];
        sample.zGTEx = NaN(length(sample.geneId), 1);
        sample.zGTEx(sampGeneIdx) = modzscore(sampleLogCount, gtexLogRead, 2);
        
        for i = 1:3
            GTEx.pctSampleExpressed(:,i) = sum(GTEx.discreteExpLevel(:,validSample)==i,2) ./ sum(validSample);
        end
        GTEx.pctSampleExpressed = GTEx.pctSampleExpressed * 100;
        gtexDataHeader = {'GTEx_gene', '%GTEx_noExp', '%GTEx_poorExp', '%GTEx_Exp', 'log2_ratio_to_GTEx', 'GTEx_medianLogCount', 'zscore_GTEx'};
        gtexData = repmat({'NA'}, nvar, length(gtexDataHeader));
        
        %use bedtools to query overlap genes
        matchGene = GENOMEFUNC.searchAllOverlapGTF( variant.loc, para.gtexGtfFn, ...
            'padChr', false, 'mt', 'MT', 'geneIdFormat', 'ENSG[\w\.]+');
        [matchCount, uniqMatch] = eleCounts(matchGene.locidx_start);
        uniqIdx = ismember(variant.locidx, uniqMatch(matchCount == 1));
        [~, uniqMatchIdx] = ismember(variant.locidx(uniqIdx), matchGene.locidx_start);
        multiIdx = find( ismember(variant.locidx, uniqMatch(matchCount > 1)) );
        [~, gencodeIdx2Sample] = GENOMEFUNC.ismemberGencode(matchGene.geneId, sample.geneId);
        [~, gencodeIdx2GTEx] = ismember(matchGene.geneId, GTEx.gencodeId);
        
        %variant mapped to only one gene
        gidx = gencodeIdx2GTEx(uniqMatchIdx);        
        gtexData(uniqIdx, 1) = GTEx.gene(gidx);
        for expi = 1:3
            gtexData(uniqIdx, expi+1) = numarray2strarray(GTEx.pctSampleExpressed(gidx,expi));
        end
        
        gidx = gencodeIdx2Sample(uniqMatchIdx);        
        exp2Gtex = NaN(length(gidx),1);
        exp2Gtex(gidx~=0) = sample.logRatio2GTEx(gidx(gidx~=0), 1);
        gtexData(uniqIdx, 5) = numarray2strarray(exp2Gtex);
        exp2Gtex(gidx~=0) = sample.logRatio2GTEx(gidx(gidx~=0), 2);
        gtexData(uniqIdx, 6) = numarray2strarray(exp2Gtex);
        exp2Gtex(gidx~=0) = sample.zGTEx(gidx(gidx~=0));
        gtexData(uniqIdx, 7) = numarray2strarray(exp2Gtex);
        
        %variant mapped to multiple genes
        for i = 1:length(multiIdx)
            gidx = gencodeIdx2GTEx( matchGene.locidx_start == variant.locidx(multiIdx(i)) );
            gtexData{multiIdx(i), 1} = strjoin(GTEx.gene(gidx)', ', ');
            for expi = 1:3
                gtexData{multiIdx(i), expi+1} = strjoin(numarray2strarray(GTEx.pctSampleExpressed(gidx,expi))', ', ');
            end
            
            gidx = gencodeIdx2Sample( matchGene.locidx_start == variant.locidx(multiIdx(i)) );
            if ~isempty(gidx)
                exp2Gtex = NaN(length(gidx),1);
                exp2Gtex(gidx~=0) = sample.logRatio2GTEx(gidx(gidx~=0), 1);
                gtexData{multiIdx(i),5} = strjoin(numarray2strarray(exp2Gtex)', ', ');
                exp2Gtex(gidx~=0) = sample.logRatio2GTEx(gidx(gidx~=0), 2);
                gtexData{multiIdx(i),6} = strjoin(numarray2strarray(exp2Gtex)', ', ');
                exp2Gtex(gidx~=0) = sample.zGTEx(gidx(gidx~=0));
                gtexData{multiIdx(i),7} = strjoin(numarray2strarray(exp2Gtex)', ', ');
            end
        end
        gtexData(:, 2:end) = strrep(gtexData(:, 2:end), 'NaN', 'NA');
    else
        gtexData = {};
    end    
end


function [aseDataHeader, aseData] = runTestASE(data, dnaData, para)
    nvar = size(data,1);
%     testIdx = ~any(isnan(data),2) & data(:,2) ~= 0; %remove not transcribed alt
            
%     testRes.AEF = data(:, 2) ./ sum(data(:, [1,2,4]), 2);
    testRes.AEF = data(:, 2) ./ sum(data(:, [1,2]), 2);

    %test with purity estimate
    res = testAseWithPurityEstimate(data(:,2), sum(data(:, 1:2), 2), ...
        'minAltRead', para.minAltReadForASE, 'minRead', para.minTotalReadForASE, ...
        'minAEF', para.minAEFForASE, 'minN', para.minNumVarForASE);
    testRes.purExpectedAEF = repmat( res.expectedAEF, nvar, 1);
    testRes.purPvalAltOverExp = res.purityPvalAltOverExp;
    
    %test with comparison to DNA data
    res = testAseWithDna( data(:, 1:2), dnaData, ...
        'minRnaRead', para.minTotalReadForASE); 
    testRes.vafDiff = res.vafDiff;
    testRes.vafPvalAltOverExpBino = res.vafDnaPvalAltOverExp;
    testRes.vafPvalAltOverExpNorm = res.vafDiffPvalAltOverExp;

%     lookup = {'AEF', 'allele_expressed_frequency(AEF)'; ...
%         'expectedAEF', 'estimated_expected_AEF'; ...
%         'pvalue', 'pvalue_allelic_expression'; ...
%         'evalue', 'evalue_allelic_expression'};
    lookup = {'AEF', 'RNA_VAF'; ...
        'purExpectedAEF', 'estimated_expected_RNA_VAF'; ...
        'purPvalAltOverExp', 'pval_alt_overexp_by_RNA_VAF'; ...
        'vafDiff', 'VAF_diff_DNA_vs_RNA'; ...
        'vafPvalAltOverExpBino', 'pval_alt_overexp_vs_DNAVAF_bino'; ...
        'vafPvalAltOverExpNorm', 'pval_alt_overexp_vs_DNAVAF_norm'};
    fds = fieldnames(testRes);
    fds = fds(ismember(fds, lookup(:,1)));
    aseData = NaN(nvar, length(fds));
    [~, i] = ismember(fds, lookup(:,1));
    aseDataHeader = lookup(i, 2)';
    for i = 1:length(fds)
        aseData(:, i) = testRes.(fds{i});
    end
    aseData = strrep(numarray2strarray(aseData), 'NaN', 'NA');
end

function subtypeinfo = readSubtypeInfo(subtypeFn)
    %return a struct: subtype.subtype1 = {sample1, sample2...}
    f = fopen(subtypeFn, 'r');
    line = fgetl(f);
    while ischar(line)
        t = textscan(line, '%s');
        t = t{1};
        subtypeinfo.(strrep(t{1},'-','_')) = t(2:end);
        line = fgetl(f);
    end
    fclose(f);
end
