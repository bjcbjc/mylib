
load data/TCGA_UCEC.mat
%%

CnHighSample = UCEC.clinical.sample(strcmp(UCEC.clinical.IntegrativeCluster, 'CN high'));

%% RTK, limit to CN high

gene = {'FGFR2', 'FGFR1', 'ERBB2', 'ERBB3', 'ERBB4', 'EGFR', 'MET', ...
    'KIT', 'IGF1R', 'PDGFRA', 'PDGFRB', 'KDR', 'FLT4', 'ALK'};

validGene = ismember(UCEC.RTK.geneName, gene);
%test (a) amp (b) amp or overexp (c ) amp and overexp (d) overexp ...
%(e) mut or amp  (f) mut or amp or overexp
expThres = 3;
testLabel = {'amp', any(UCEC.RTK.altamp(validGene,:), 1); ...    
    'ampAndOverexp', any(UCEC.RTK.altamp(validGene,:) & UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    'ampOrOverexp', any(UCEC.RTK.altamp(validGene,:) | UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    'overexp', any(UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    };

%     'mutOrAmp', any(UCEC.RTK.mut(validGene,:) | UCEC.RTK.altamp(validGene,:), 1); ...
%     'mutOrAmpOrOverexp', any(UCEC.RTK.mut(validGene,:) | ...
%     UCEC.RTK.altamp(validGene,:) | UCEC.RTK.rnaZ(validGene,:) >= expThres, 1); ...
%     'mutOrAmpAndOverexp', any(UCEC.RTK.mut(validGene,:) | ...
%     (UCEC.RTK.altamp(validGene,:) & UCEC.RTK.rnaZ(validGene,:) >= expThres),1) ...
%     };

for testIdx = 1:size(testLabel,1)
    validSample = ismember(UCEC.RTK.sample', CnHighSample);
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'mut'))
        validSample = validSample & UCEC.RTK.validMutSample;
    end
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'amp'))
        validSample = validSample & UCEC.RTK.validCnSample;
    end
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'exp'))
        validSample = validSample & UCEC.RTK.validRnaSample;
    end
        
    alteredSample = UCEC.RTK.sample(testLabel{testIdx, 2} & validSample);
    nonAlteredSample = UCEC.RTK.sample( ~testLabel{testIdx, 2} & validSample);

    maxMonth = 120; %[60, 120, 250];
    fig = offFigure(1200, 450);
    for i = 1:length(maxMonth)
        %OS
        subplot(1,2,1);
        survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
            strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
        idx = survivalData(:,1) > maxMonth(i);
        survivalData(idx,1) = maxMonth(i);
        survivalData(idx,2) = 1;
        [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
            alteredSample, nonAlteredSample, 'survival');

        %DFS
        subplot(1,2,2);
        survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
            strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
        survivalData(idx,1) = maxMonth(i);
        survivalData(idx,2) = 1;
        [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
            alteredSample, nonAlteredSample, 'disease free');
        saveas(fig, sprintf('figure/survival/suvival_RTK_%s_max_%dMonth.CNHigh.eps',...
            testLabel{testIdx,1},maxMonth(i)), 'epsc');
    end
    close(fig);
end

%% RTK, all samples

gene = {'FGFR2', 'FGFR1', 'ERBB2', 'ERBB3', 'ERBB4', 'EGFR', 'MET', ...
    'KIT', 'IGF1R', 'PDGFRA', 'PDGFRB', 'KDR', 'FLT4', 'ALK'};

validGene = ismember(UCEC.RTK.geneName, gene);

expThres = 3;
testLabel = {'amp', any(UCEC.RTK.altamp(validGene,:), 1); ...    
    'ampAndOverexp', any(UCEC.RTK.altamp(validGene,:) & UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    'mut', any(UCEC.RTK.mut(validGene,:), 1); ...
    };

%     'mutOrAmp', any(UCEC.RTK.mut(validGene,:) | UCEC.RTK.altamp(validGene,:), 1); ...
%     'mutOrAmpOrOverexp', any(UCEC.RTK.mut(validGene,:) | ...
%     UCEC.RTK.altamp(validGene,:) | UCEC.RTK.rnaZ(validGene,:) >= expThres, 1); ...
%     'mutOrAmpAndOverexp', any(UCEC.RTK.mut(validGene,:) | ...
%     (UCEC.RTK.altamp(validGene,:) & UCEC.RTK.rnaZ(validGene,:) >= expThres),1) ...
%     };

for testIdx = 1:size(testLabel,1)
    validSample = true(1, length(UCEC.RTK.sample));
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'mut'))
        validSample = validSample & UCEC.RTK.validMutSample;
    end
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'amp'))
        validSample = validSample & UCEC.RTK.validCnSample;
    end
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'exp'))
        validSample = validSample & UCEC.RTK.validRnaSample;
    end
        
    alteredSample = UCEC.RTK.sample(testLabel{testIdx, 2} & validSample);
    nonAlteredSample = UCEC.RTK.sample( ~testLabel{testIdx, 2} & validSample);

    maxMonth = 120; %[60, 120, 250];
    fig = offFigure(1200, 450);
    for i = 1:length(maxMonth)
        %OS
        subplot(1,2,1);
        survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
            strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
        idx = survivalData(:,1) > maxMonth(i);
        survivalData(idx,1) = maxMonth(i);
        survivalData(idx,2) = 1;
        [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
            alteredSample, nonAlteredSample, 'survival');

        %DFS
        subplot(1,2,2);
        survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
            strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
        survivalData(idx,1) = maxMonth(i);
        survivalData(idx,2) = 1;
        [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
            alteredSample, nonAlteredSample, 'disease free');
        saveas(fig, sprintf('figure/survival/suvival_RTK_%s_max_%dMonth.all.eps',...
            testLabel{testIdx,1},maxMonth(i)), 'epsc');
    end
    close(fig);
end

%% RTK, each subtype

gene = {'FGFR2', 'FGFR1', 'ERBB2', 'ERBB3', 'ERBB4', 'EGFR', 'MET', ...
    'KIT', 'IGF1R', 'PDGFRA', 'PDGFRB', 'KDR', 'FLT4', 'ALK'};

validGene = ismember(UCEC.RTK.geneName, gene);

expThres = 3;
testLabel = { ...
    'mut', any(UCEC.RTK.mut(validGene,:), 1); ...
    };
subtypes = {'CN high', 'MSI', 'CN low', 'POLE'};
for testIdx = 1:size(testLabel,1)
    for subtypeIdx = 1:length(subtypes)
        validSample = ismember(UCEC.RTK.sample', ...
            UCEC.clinical.sample(strcmp(UCEC.clinical.IntegrativeCluster, subtypes{subtypeIdx})));
        if ~isempty(strfind(lower(testLabel{testIdx,1}), 'mut'))
            validSample = validSample & UCEC.RTK.validMutSample;
        end
        if ~isempty(strfind(lower(testLabel{testIdx,1}), 'amp'))
            validSample = validSample & UCEC.RTK.validCnSample;
        end
        if ~isempty(strfind(lower(testLabel{testIdx,1}), 'exp'))
            validSample = validSample & UCEC.RTK.validRnaSample;
        end
        
        alteredSample = UCEC.RTK.sample(testLabel{testIdx, 2} & validSample);
        nonAlteredSample = UCEC.RTK.sample( ~testLabel{testIdx, 2} & validSample);
        
        maxMonth = 120; %[60, 120, 250];
        fig = offFigure(1200, 450);
        for i = 1:length(maxMonth)
            %OS
            subplot(1,2,1);
            survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
                strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
            idx = survivalData(:,1) > maxMonth(i);
            survivalData(idx,1) = maxMonth(i);
            survivalData(idx,2) = 1;
            [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
                alteredSample, nonAlteredSample, 'survival');
            
            %DFS
            subplot(1,2,2);
            survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
                strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
            survivalData(idx,1) = maxMonth(i);
            survivalData(idx,2) = 1;
            [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
                alteredSample, nonAlteredSample, 'disease free');
            saveas(fig, sprintf('figure/survival/suvival_RTK_%s_max_%dMonth.%s.eps',...
                testLabel{testIdx,1},maxMonth(i), strrep(subtypes{subtypeIdx},' ','')), 'epsc');
        end
        close(fig);
    end
end


%% ERBB2 

gene = {'ERBB2'};
validGene = ismember(UCEC.RTK.geneName, gene);
%test (a) amp (b) amp or overexp (c ) amp and overexp (d) overexp ...
%(e) mut or amp  (f) mut or amp or overexp
expThres = 3;
testLabel = {'amp', any(UCEC.RTK.altamp(validGene,:), 1); ...
    'ampOrOverexp', any(UCEC.RTK.altamp(validGene,:) | UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    'ampAndOverexp', any(UCEC.RTK.altamp(validGene,:) & UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    'overexp', any(UCEC.RTK.rnaZ(validGene,:) >= expThres,1); ...
    'mutOrAmp', any(UCEC.RTK.mut(validGene,:) | UCEC.RTK.altamp(validGene,:), 1); ...
    'mutOrAmpOrOverexp', any(UCEC.RTK.mut(validGene,:) | ...
    UCEC.RTK.altamp(validGene,:) | UCEC.RTK.rnaZ(validGene,:) >= expThres, 1); ...
    };

for testIdx = 1:size(testLabel,1)
    validSample = true(1, length(UCEC.RTK.sample));
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'mut'))
        validSample = validSample & UCEC.RTK.validMutSample;
    end
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'amp'))
        validSample = validSample & UCEC.RTK.validCnSample;
    end
    if ~isempty(strfind(lower(testLabel{testIdx,1}), 'exp'))
        validSample = validSample & UCEC.RTK.validRnaSample;
    end
        
    alteredSample = UCEC.RTK.sample(testLabel{testIdx, 2} & validSample);
    nonAlteredSample = UCEC.RTK.sample( ~testLabel{testIdx, 2} & validSample);

    maxMonth = [60, 120, 250];
    fig = offFigure(1200, 450);
    for i = 1:length(maxMonth)
        %OS
        subplot(1,2,1);
        survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
            strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
        idx = survivalData(:,1) > maxMonth(i);
        survivalData(idx,1) = maxMonth(i);
        survivalData(idx,2) = 1;
        [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
            alteredSample, nonAlteredSample, 'survival');

        %DFS
        subplot(1,2,2);
        survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
            strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
        survivalData(idx,1) = maxMonth(i);
        survivalData(idx,2) = 1;
        [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
            alteredSample, nonAlteredSample, 'disease free');
        saveas(fig, sprintf('figure/survival/suvival_ERBB2_%s_max_%dMonth.eps',...
            testLabel{testIdx,1},maxMonth(i)), 'epsc');
    end
    close(fig);
end

%% ERBB2 RPPA
addpath('cgds_toolbox/')
cgdsURL = 'http://www.cbioportal.org/public-portal';
profileData = getprofiledata(cgdsURL, 'ucec_tcga_pub_all',...
    {'ucec_tcga_pub_RPPA_protein_level' }, ...
    {'ERBB2','ERBB2_PY1248'},true);
expThres = 1.5; 
alteredSample = profileData.caseId( ...
    any(profileData.data >= expThres, 1) );
nonAlteredSample = profileData.caseId( ...
    ~any(profileData.data >= expThres, 1) & ~any(isnan(profileData.data),1) );

maxMonth = [60, 120, 500];
fig = offFigure(1200, 450);
for i = 1:length(maxMonth)
    %OS
    subplot(1,2,1);
    survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
    idx = survivalData(:,1) > maxMonth(i);
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'survival');
    
    %DFS
    subplot(1,2,2);
    survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'disease free');
    saveas(fig, sprintf('figure/survival/suvival_ERBB2_overRPPA%0.1f_max_%dMonth.eps',expThres,maxMonth(i)), 'epsc');
end
close(fig);

%% JNK
[~,gidx] = ismember({'MAPK8', 'MAPK9','MAPK10'}, UCEC.cbio.geneName);

expThres = 3; 
alteredSample = UCEC.cbio.sample( ...
     any(UCEC.cbio.rnaz(gidx,:) >= expThres,1) );
nonAlteredSample = UCEC.cbio.sample( ...
     all(UCEC.cbio.rnaz(gidx,:) < expThres,1) );

maxMonth = [60, 120, 500];
fig = offFigure(1200, 450);
for i = 1:length(maxMonth)
    %OS
    subplot(1,2,1);
    survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
    idx = survivalData(:,1) > maxMonth(i);
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'survival');
    
    %DFS
    subplot(1,2,2);
    survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'disease free');
    saveas(fig, sprintf('figure/survival/suvival_JNK_overexp%0.1f_max_%dMonth.eps',expThres,maxMonth(i)), 'epsc');
end
close(fig);


%% JNK RPPA

profileData = getprofiledata(cgdsURL, 'ucec_tcga_pub_all',...
    {'ucec_tcga_pub_RPPA_protein_level', }, ...
    {'MAPK9'},true);
expThres = 1.5; 
alteredSample = profileData.caseId( ...
    any(profileData.data >= expThres, 1) );
nonAlteredSample = profileData.caseId( ...
    ~any(profileData.data >= expThres, 1) & ~any(isnan(profileData.data),1) );

maxMonth = [60, 120, 500];
fig = offFigure(1200, 450);
for i = 1:length(maxMonth)
    %OS
    subplot(1,2,1);
    survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
    idx = survivalData(:,1) > maxMonth(i);
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'survival');
    
    %DFS
    subplot(1,2,2);
    survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'disease free');
    saveas(fig, sprintf('figure/survival/suvival_JNK_overRPPA%0.1f_max_%dMonth.eps',expThres,maxMonth(i)), 'epsc');
end
close(fig);

%%
rppaIdx = find(strcmp(UCEC.RppaRbn.protein,'JNK2'));
z = zeromean_univar_normalization(UCEC.RppaRbn.data(rppaIdx,:),2);

expThres = 2; 
alteredSample = TCGASampleDecoder.reduceToParticipant(UCEC.RppaRbn.sample( z>expThres ));
nonAlteredSample = TCGASampleDecoder.reduceToParticipant(UCEC.RppaRbn.sample( z<expThres ));

maxMonth = [60, 120, 500];
fig = offFigure(1200, 450);
for i = 1:length(maxMonth)
    %OS
    subplot(1,2,1);
    survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
    idx = survivalData(:,1) > maxMonth(i);
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'survival');
    
    %DFS
    subplot(1,2,2);
    survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'disease free');
    saveas(fig, sprintf('figure/suvival_JNK_MD_RPPA%0.1f_max_%dMonth.eps',expThres,maxMonth(i)), 'epsc');
end
close(fig);


%% RTK RPPA from cbio
addpath('cgds_toolbox/')
cgdsURL = 'http://www.cbioportal.org/public-portal';
proteinId = {'ERBB2', 'ERBB2_PY1248', 'ERBB3', 'ERBB3_PY1298', ...
    'KIT', 'MET' , 'MET_PY1235' , 'IGF1R' , 'KDR' , ...
    'EGFR' , 'EGFR_PY1068' , 'EGFR_PY1173' , 'EGFR_PY992' };
profileData = getprofiledata(cgdsURL, 'ucec_tcga_pub_all',...
    {'ucec_tcga_pub_RPPA_protein_level', }, ...
    proteinId,true);
expThres = 3.5; 
alteredSample = profileData.caseId( ...
    any(profileData.data >= expThres, 1) );
nonAlteredSample = profileData.caseId( ...
    ~any(profileData.data >= expThres, 1) & ~any(isnan(profileData.data),1) );

maxMonth = [60, 120, 500];
fig = offFigure(1200, 450);
for i = 1:length(maxMonth)
    %OS
    subplot(1,2,1);
    survivalData = [str2double(UCEC.cbio.clinical.OS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.OS_STATUS,'living')];
    idx = survivalData(:,1) > maxMonth(i);
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'survival');
    
    %DFS
    subplot(1,2,2);
    survivalData = [str2double(UCEC.cbio.clinical.DFS_MONTHS), ...
        strcmpi(UCEC.cbio.clinical.DFS_STATUS,'diseasefree')];
    survivalData(idx,1) = maxMonth(i);
    survivalData(idx,2) = 1;
    [curvehandle, pval] = plotKM(survivalData, UCEC.cbio.sample, ...
        alteredSample, nonAlteredSample, 'disease free');
    saveas(fig, sprintf('figure/survival/suvival_RTK_overRPPA%0.1f_max_%dMonth.eps',expThres,maxMonth(i)), 'epsc');
end
close(fig);

