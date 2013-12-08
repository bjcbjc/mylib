classdef TCGAclinical < handle

    properties (Constant)
        folder = './data/clinical/';
        pafn = 'clinical_patient_all_brca.txt';
        sdfn = 'clinical_slide_all_brca.txt';
    end
    
    methods (Static)
        function s = readAllClinical()
            s.patient = TCGAclinical.readPatients();
            s.slide = TCGAclinical.readSlides();
        end
        
        function patientinfo = readPatients(fn)
            if nargin < 1
                fn = [TCGAclinical.folder TCGAclinical.pafn];
            end
            t = parseText(fn, 'nrowname',1,'ncolname',1);
            t.colname(1) = [];

            %discard columns w/o information or lack of variation within tumor or
            %normal samples
            remove = all(strcmpi('null', t.text), 1);
            remove = remove | all(~cellfun(@isempty, strfind(t.text, ...
                'Not Available')), 1);    
            t.colname(remove) = [];
            t.text(:,remove) = [];

            remove = strcmp(t.colname, 'patient_id') | ...
                strcmp(t.colname, 'bcr_patient_uuid') | ...
                strcmp(t.colname, 'tissue_source_site');
            t.colname(remove) = [];
            t.text(:, remove) = [];

            t.text = strrep(t.text, '[Not Available]', 'NaN');
            t.text = strrep(t.text, '[Not Applicable]', 'NaN');
            t.text = strrep(t.text, '[Completed]', 'NaN');
            t.text = strrep(t.text, '[Discrepancy]', 'NaN');

            remove = [];
            for i = 1:size(t.text, 2)
                if all(strcmp(t.text(:,i), t.text{1, i})) || ...
                        length(unique(setdiff(t.text(:,i), 'NaN'))) == 1
                    remove = [remove i];
                end
            end
            t.colname(remove) = [];
            t.text(:, remove) = [];

            sample = tcgaSampleDecoder.decode(t.rowname);

            patientinfo.sample = t.rowname;
            patientinfo.tss = sample.tss;
            patientinfo.participant = sample.participant; 
            %separate attributes into ref, categorical and numerical
            %ref: most data are not available or censored (can't be used for test)
            %categorical: categorical data
            %numerical
            %
            missratio = mean(strcmp(t.text, 'NaN'), 1);
            numeric = all(isnumericstring(t.text),1);
            daysattr = ~cellfun(@isempty, strfind(t.colname, 'days_'))';
            patientinfo.attr_num = ...
                t.colname((numeric & missratio <= 0.5) | daysattr);
            patientinfo.attr_ref = setdiff(t.colname(missratio > 0.5), patientinfo.attr_num);
            patientinfo.attr_cat = setdiff(t.colname, union(patientinfo.attr_ref, patientinfo.attr_num));

            [~,i] = ismember(patientinfo.attr_num, t.colname);
            patientinfo.data_num = str2double( t.text(:, i) );
            [~,i] = ismember(patientinfo.attr_ref, t.colname);
            patientinfo.data_ref = t.text(:, i);
            [~,i] = ismember(patientinfo.attr_cat, t.colname);
            patientinfo.data_cat = t.text(:, i);

            %stage features, big bins
            stage_attr = {'breast_tumor_clinical_m_stage', ...
                'breast_tumor_pathologic_grouping_stage', ...
                'breast_tumor_pathologic_n_stage', ...
                'breast_tumor_pathologic_t_stage'};
            patientinfo.attr_num(end+1:end+4) = strcat(stage_attr, '_num');

            [~,i] = ismember(stage_attr, t.colname);
            tmp = t.text(:, i);

            %anything with X -> NaN
            i = ~cellfun(@isempty, strfind(tmp, 'X'));
            tmp(i) = {'NaN'};
            stagetype = {'M', '', 'N', 'T'};
            for i = [1 3 4]
                for j = 0:4
                    k = ~cellfun(@isempty, ...
                        strfind(tmp(:,i), sprintf('%s%d',stagetype{i},j)));
                    tmp(k,i) = {sprintf('%d',j)};
                end
            end
            group = {'IV', 'III', 'II', 'I'}; %reverse search!
            for i = 1:length(group)
                k = ~cellfun(@isempty, ...
                    strfind(tmp(:,2), sprintf('Stage %s',group{i})));
                tmp(k,2) = {sprintf('%d',5-i)};
            end
            patientinfo.data_num(:, end+1:end+4) = str2double(tmp);
        end

        function slideinfo = readSlides(fn)
            if nargin < 1
                fn = [TCGAclinical.folder TCGAclinical.sdfn];
            end
            t = parseText(fn, 'nrowname',1,'ncolname',1);
            t.colname(1) = [];

            %discard columns w/o information or lack of variation within tumor or
            %normal samples
            remove = all(strcmpi('null', t.text), 1);
            remove = remove | all(~cellfun(@isempty, strfind(t.text, ...
                'Not Available')), 1);    
            t.colname(remove) = [];
            t.text(:,remove) = [];

            sample = tcgaSampleDecoder.decode(t.rowname, false);
            remove = [];
            normsamp = find(sample.sampletype == 1);
            tumorsamp = find(sample.sampletype == 2);
            for i = 1:size(t.text,2)
                if all(strcmp(t.text(normsamp,i), t.text{normsamp(1),i})) && ...
                        all(strcmp(t.text(tumorsamp,i), t.text{tumorsamp(1),i}))
                    remove = [remove i];
                end
            end
            t.colname(remove) = [];
            t.text(:,remove) = [];

            remove = strcmp(t.colname, 'bcr_slide_uuid');
            t.colname(remove) = [];
            t.text(:,remove) = [];

            remove = strcmp(t.colname, 'bcr_slide_barcode');
            t.colname(remove) = [];
            t.text(:,remove) = [];

            loc = t.text(:, strcmp(t.colname, 'section_location'));
            t.text(:, strcmp(t.colname, 'section_location')) = [];
            t.colname(strcmp(t.colname, 'section_location')) = [];

            %for the rest columns, create three matrices (top, bottom, mean) for
            %tumor and normal samples
            %
            %normindex = sample.sampletype == 1;

            sampfd = {'_norm', ''};
            copiedfield = {'tss', 'participant', 'samplevial', 'sample'};
            for idx = 1:2
                sampidx = sample.sampletype == idx;
                for cfi = 1:length(copiedfield)
                    if ~strcmp(copiedfield{cfi}, 'sample')
                        slideinfo.([copiedfield{cfi} sampfd{idx}]) = sample.(copiedfield{cfi})(sampidx);
                    else
                        slideinfo.([copiedfield{cfi} sampfd{idx}]) = t.rowname(sampidx);
                    end
                end
                [~, ui] = uniqueCellRows( [ slideinfo.([copiedfield{1} sampfd{idx}]) ...
                    slideinfo.([copiedfield{2} sampfd{idx}])  slideinfo.([copiedfield{3} sampfd{idx}]) ] );
                for cfi = 1:length(copiedfield)
                    slideinfo.([copiedfield{cfi} sampfd{idx}]) = ...
                        slideinfo.([copiedfield{cfi} sampfd{idx}])(ui);
                end
                %slideinfo.sample_norm = uniqueCellRows(slideinfo.sample_norm);
            end
                
            
%             slideinfo.sample = [sample.tss(~normindex) ...
%                 sample.participant(~normindex) sample.samplevial(~normindex)];
%             slideinfo.sample_norm = [sample.tss(normindex) ...
%                 sample.participant(normindex) sample.samplevial(normindex)];


            slideinfo.slide_attr = t.colname;

            %top, bottom, mean
            slideinfo.slide_data = NaN(length(t.colname), length(slideinfo.sample), 3);
            slideinfo.slide_data_norm = NaN(length(t.colname), length(slideinfo.sample_norm), 3);

            %top, bottom will match participants, samples and vials
            %but mean will only match particpants and sampletype (tumor or normal)

            
            dimfilter = [strcmpi(loc, 'top'), strcmpi(loc, 'bottom'), true(size(loc))];

            for sampi = 1:length(sampfd)
                smp = ['sample' sampfd{sampi}];
                d = ['slide_data' sampfd{sampi}];
                for dim = 1:3
                    for i = 1:length(slideinfo.(smp))
                        rowi = strcmp(sample.participant, slideinfo.(['participant' sampfd{sampi}]){i} ) & ...
                            strcmp(sample.samplevial, slideinfo.(['samplevial' sampfd{sampi}]){i} ) & ...
                            dimfilter(:, dim);
%                         rowi = strcmp(sample.participant, slideinfo.(smp){i,2}) & ...
%                             strcmp(sample.samplevial, slideinfo.(smp){i,3}) & dimfilter(:,dim);
                        if any(rowi)
                            slideinfo.(d)(:,i,dim) = nanmean( str2double(t.text(rowi,:)), 1);
                        end
                    end
                end
            end    

        end
        
        function res = testClinicalAssociation(clinicalData, testSample, testData)            
            %testSample: cell array
            %testData: #testSample x #features for test
            %clinicalData: patientinfo 
            %
            %
            [csampi, tsampi] = tcgaSampleDecoder.matchSample(clinicalData.sample, testSample);            
            clinicalData.data_num = clinicalData.data_num(csampi, :);
            clinicalData.data_ref = clinicalData.data_ref(csampi, :);
            clinicalData.data_cat = clinicalData.data_cat(csampi, :);
            testData = testData(tsampi, :);
            ntestVar = size(testData, 2);
            %test data type to determine tests: 1: numeric, 2: binary, n:
            %#categories
            testDataType = zeros(ntestVar, 1);
            if iscell(testData)
                for i = 1:ntestVar
                    vi = ~strcmp(testData(:, i), 'NaN') & ~cellfun(@isempty, testData(:,i));
                    testDataType(i) = length(unique( testData(vi, i) ) );
                end
                testDataType( testDataType == 1 ) = 0; %no test
            else
                testDataType(:) = 1;
            end
            %numeric clinical data
            res.pval_num = NaN(ntestVar, length(clinicalData.attr_num));
            if any(testDataType == 1)
                res.rho_num = NaN(ntestVar, length(clinicalData.attr_num));
                testi = testDataType == 1;
            
                [res.rho_num(testi, :) res.pval_num(testi,:)] = corr( ...
                    testData(:, testi),  clinicalData.data_num, ...
                    'rows', 'pairwise');
            end
            testi = find(testDataType >= 2);
            for i = 1:length(testi)
                for j = 1:length(clinicalData.attr_num)
                    vi = find( ~strcmp(testData(:, testi(i)), 'NaN') & ...
                        ~strcmp(testData(:, testi(i)), '') & ...
                        ~isnan(clinicalData.data_num(:, j)) );
                    if testDataType(testi(i)) == 2
                        uval = unique(testData(vi, testi(i)));
                        group1 = testData(vi, testi(i)) == uval(1);
                        res.pval_num(testi(i), j) = ranksum( ...
                            clinicalData.data_num( vi( group1), j), ...
                            clinicalData.data_num( vi( ~group1), j) );
                    else
                        res.pval_num(testi(i), j) = kruskalwallis( ...
                            clinicalData.data_num( vi, j ), ...
                            testData( vi, testi(i) ), 'off' );
                    end
                end
            end
                        
            %ref: categorical
            res.pval_ref = NaN(ntestVar, length(clinicalData.attr_ref) );
            if all(testDataType == 1) %kruskal-wallis
                for i = 1:ntestVar
                    vi1 = ~isnan(testData(:, i));
                    for j = 1:length(clinicalData.attr_ref)
                        vi = vi1 & ~strcmp(clinicalData.data_ref(:, j), 'NaN') & ...
                            ~strcmp(clinicalData.data_ref(:, j), '') ;
                        res.pval_ref(i,j) = kruskalwallis( testData(vi, i), ...
                            clinicalData.data_ref(vi, j), 'off');
                    end
                end
            else
                for i = 1:ntestVar
                    vi1 = ~strcmp(testData(:, i), 'NaN') & ...
                        ~strcmp(testData(:, i), '');
                    for j = 1:length(clinicalData.attr_ref)
                        vi = vi1 & ~strcmp(clinicalData.data_ref(:,j), 'NaN') & ...
                            ~strcmp(clinicalData.data_ref(:,j), '');
                        [~, ~, res.pval_ref(i,j)] = ...
                            crosstab(testData(vi, i), clinicalData.data_ref(vi, j));
                    end
                end
            end
            
            
            %cat
            res.pval_cat = NaN(ntestVar, length(clinicalData.attr_cat) );
            if all(testDataType == 1) %kruskal-wallis
                for i = 1:ntestVar
                    vi1 = ~isnan(testData(:, i));
                    for j = 1:length(clinicalData.attr_cat)
                        vi = vi1 & ~strcmp(clinicalData.data_cat(:, j), 'NaN') & ...
                            ~strcmp(clinicalData.data_cat(:, j), '') ;
                        res.pval_cat(i,j) = kruskalwallis( testData(vi, i), ...
                            clinicalData.data_cat(vi, j), 'off');
                    end
                end
            else
                for i = 1:ntestVar
                    vi1 = ~strcmp(testData(:, i), 'NaN') & ...
                        ~strcmp(testData(:, i), '');
                    for j = 1:length(clinicalData.attr_cat)
                        vi = vi1 & ~strcmp(clinicalData.data_cat(:,j), 'NaN') & ...
                            ~strcmp(clinicalData.data_cat(:,j), '');
                        [~, ~, res.pval_cat(i,j)] = ...
                            crosstab(testData(vi, i), clinicalData.data_cat(vi, j));
                    end
                end
            end
                
        end
        
        function res = testSlideAssociation(slideData, testSample, testData)
            %testSample: cell array
            %testData: #testSample x #features for test
            %slideData: slide
            %
            %
            [csampi, tsampi] = tcgaSampleDecoder.matchSample(slideData.sample, testSample);            
            slideData.slide_data = slideData.slide_data(:, csampi, :);            
            testData = testData(tsampi, :);
            ntestVar = size(testData, 2);
            %test data type to determine tests: 1: numeric, 2: binary, n:
            %#categories
            testDataType = zeros(ntestVar, 1);
            if iscell(testData)
                for i = 1:ntestVar
                    vi = ~strcmp(testData(:, i), 'NaN') & ~cellfun(@isempty, testData(:,i));
                    testDataType(i) = length(unique( testData(vi, i) ) );
                end
                testDataType( testDataType == 1 ) = 0; %no test
            else
                testDataType(:) = 1;
            end
            
            res.pval_num = NaN(ntestVar, length(slideData.slide_attr), size(slideData.slide_data,3));
            if all(testDataType == 1)
                res.rho_num = NaN(ntestVar, length(slideData.slide_attr), size(slideData.slide_data,3));
            end
            testi = testDataType == 1;
            for slidei = 1:size(slideData.slide_data, 3)
                [res.rho_num(testi, :, slidei) res.pval_num(testi,:, slidei)] = ...
                    corr( testData(:, testi),  squeeze(slideData.slide_data(:,:,i))', ...
                    'rows', 'pairwise');
            end
            testi = find(testDataType >= 2);
            for slidei = 1:size(slideData.slide_data, 3)
                for i = 1:length(testi)
                    for j = 1:length(slideData.slide_attr)
                        vi = find( ~strcmp(testData(:, testi(i)), 'NaN') & ...
                            ~strcmp(testData(:, testi(i)), '') & ...
                            ~isnan(squeeze(slideData.slide_data(j, :, slidei))') );
                        if testDataType(testi(i)) == 2
                            uval = unique(testData(vi, testi(i)));
                            group1 = testData(vi, testi(i)) == uval(1);
                            res.pval_num(testi(i), j, slidei) = ranksum( ...
                                squeeze(slideData.slide_data( j, vi( group1), slidei)), ...
                                squeeze(slideData.slide_data( j, vi( ~group1), slidei)) );
                        else
                            res.pval_num(testi(i), j) = kruskalwallis( ...
                                squeeze(slideData.slide_data( j, vi, slidei )), ...
                                testData( vi, testi(i) ), 'off' );
                        end
                    end
                end
            end
        end
        
    end
    
end






