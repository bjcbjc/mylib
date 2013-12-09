classdef TCGA < handle
    
    properties (Constant) 
        dirhead = '/Users/bjc/Lab/Projects/Cancer/breast/data/';
        %default filenames
        cna = 'cn/BreastMay12_all_cn.txt'; 
        exp = 'exp/Aug11Breast_Agilent244kExp_wnormals.txt'; 
        rna = 'RNASeq/BreastMay12_all_exp_rpkm.txt';
        mircna = 'cn/miRNA.continuous.matrix'; 
        mirexp = 'RNASeq/Aug11BreastMirna.txt'; 
        mut = 'mut/Dec11BreastMutationMatrix.txt';
        met = ''; 
        patient = 'clinical/clinical_patient_all_brca.txt'; 
        slide = 'clinical/clinical_slide_all_brca.txt';
        batch = 'BreastAug11_samplebatch.txt';
        subtype = {'BreastAug10AllHer2SampleId.txt', ...
            'BreastAug10BasalSampleId.txt', ...
            'BreastAug10Her2SampleId.txt', ...
            'BreastAug10LuminalSampleId.txt', ...
            'BreastAug10TripleNegSampleId.txt'};
    end
    
    methods (Access = public, Static)
        
        function data = getData(data, update)
            %read all data and build data structure
            %
            %data: old data structure; default is empty (update all data)
            %cancer: cancer type, also tells where to get the data
            %update: cell array indicating which types of data to update; default
            %   is to update all: {'cna','exp','mircna','mirexp','mut','met','cln'}
            %
            %            
            if nargin < 2
                update = {'cna', 'rna', 'mut', 'patient', 'slide'};%, ...
                    %'exp', 'mut', 'mircna','mirexp','met'};
            end
        
            %data includes the following fields:
            %sample (sample ID)
            %cnaid, cna (#cna x #sample)
            %expid, exp
            %mirid, mir
            %mutid, mut
            %metid, met
            %clnid1, cln (numeric)
            %clnid2, cln (categorical)
            %
            
            update = lower(update);
            for di = 1:length(update)
                if strcmp(update{di}, 'cna')
                    tmp.cna = TCGA.getCna();                    
                elseif strcmp(update{di}, 'exp')
                    tmp.exp = TCGA.getExp();
                elseif strcmp(update{di}, 'rna')
                    tmp.rna = TCGA.getRNA();
                elseif strcmp(update{di}, 'patient')
                    tmp.patient = TCGA.getPatients();
                elseif strcmp(update{di}, 'slide')
                    tmp.slide = TCGA.getSlides();
                elseif strcmp(update{di}, 'mut')
                    tmp.mut = TCGA.getMutation();
                elseif strcmp(update{di}, 'mircna')
                    tmp.mircna = TCGA.getMircna();
                elseif strcmp(update{di}, 'mirexp')
                    tmp.mirexp = TCGA.getMirexp();
                elseif strcmp(update{di}, 'met')
                    tmp.met = TCGA.getMethylation();
                else
                    error('unknown data type %s', update{di});
                end
            end

            GQ = genequery();
            samplecolumn = {};
            
            if isempty(data) %create a new structure
                allpatient = {};
                for di = 1:length(update)
                    allpatient = union(allpatient, tmp.(update{di}).sample(:,1));
                    
                    if any(tmp.(update{di}).sampletype == 2)
                        samplecolumn = union(samplecolumn, update{di});
                    end
                    if any(tmp.(update{di}).sampletype == 1)
                        samplecolumn = union(samplecolumn, [update{di}, '_normal']);
                    end                    
                end
                
                samplemtx = zeros(length(allpatient), length(samplecolumn));
                for coli = 1:length(samplecolumn)
                    fields = textscan(samplecolumn{coli}, '%s', 'delimiter', '_');
                    fields = fields{1};
                    if length(fields) == 2 %normal
                        subidx = find(tmp.(fields{1}).sampletype == 1);
                        [tf, idx] = ismember(allpatient, ...
                            tmp.(fields{1}).sample(subidx,1));                        
                    else
                        subidx = find(tmp.(fields{1}).sampletype == 2);
                        [tf, idx] = ismember(allpatient, ...
                            tmp.(fields{1}).sample(subidx,1));
                    end
                    samplemtx(tf, coli) = subidx(idx(tf));
                end
                
                % now collect all data
                data.allsample = allpatient;
                data.type = samplecolumn;
                data.typematrix = samplemtx;
                for di = 1:length(update)
                    data.([update{di} '_sample']) = tmp.(update{di}).sample;
                    if strcmp(update{di}, 'cna') || strcmp(update{di}, 'exp') ...
                            || strcmp(update{di}, 'mut') || strcmp(update{di}, 'rna')
                        loc = GQ.idquery(tmp.(update{di}).attr, {'chrom', 'start', 'end'});
                        [~,sortidx] = sortrows(loc);
                        data.([update{di} '_id']) = tmp.(update{di}).attr(sortidx);
                        data.(update{di}) = tmp.(update{di}).data(sortidx, :);                        
                    elseif strcmp(update{di}, 'patient')
                        subfds = fieldnames(tmp.(update{di}));
                        for sfi = 1:length(subfds)
                            if ~isempty(strfind(subfds{sfi}, 'data')) || ...
                                    ~isempty(strfind(subfds{sfi}, 'attr'))
                                data.([update{di} '_' subfds{sfi}]) = tmp.(update{di}).(subfds{sfi});
                            end
                        end                        
                    else
                        data.([update{di} '_attr']) = tmp.(update{di}).attr;
                        data.(update{di}) = tmp.(update{di}).data;
                    end
                end
            else %update data
                error('Not implemented yet. HaHaHa\n');
            end
            
        end
                
        function data = sampleMatch(alldata, matchnormal, types)
            % match samples between different data types
            % 
            % 
            if nargin < 3
                types = {'cna', 'exp', 'rna', 'mut', 'slide', 'patient'};
            end
            if nargin < 2
                matchnormal = false;
            end
            
            types(~ismember(types, alldata.type)) = [];
                        
            tumorcolumn = ismember(alldata.type, types);
            rowi = all(alldata.typematrix(:, tumorcolumn)~=0, 2);
            
            fds = fieldnames(alldata);
            
            for i = 1:length(types)                
                [~, typeidx] = ismember(types{i}, alldata.type);
                if i == 1
                    data.sample = ...
                        alldata.([types{i} '_sample'])(alldata.typematrix(rowi, typeidx), 1);
                    samplecopy = data.sample;
                    %if normal samples are required, sort them to be on top of
                    %sample list
                    if matchnormal
                        rowi = find(rowi);
                        normcolumn = ismember(alldata.type, strcat(types, '_normal'));
                        normsample = all(alldata.typematrix(rowi, normcolumn)~=0, 2);
                        [~, reordersample] = sort(normsample, 'descend');
                        data.sample = data.sample(reordersample);                        
                        normrowi = rowi(normsample);                    
                    end
                else %double check on samples
                    if nnz(~strcmp(samplecopy, ...
                        alldata.([types{i} '_sample'])(alldata.typematrix(rowi, typeidx),1))) > 0
                        error('some samples do not match, typematrix is corrupted');
                    end
                end
                
                %attribute names, id (labels of rows)
                cpfield = fds(~cellfun(@isempty, strfind(fds, types{i})));
                cpfield = setdiff(cpfield, [types{i} '_sample']);
                cpfield = setdiff(cpfield, types{i});
                cpfield = setdiff(cpfield, cpfield(~cellfun(@isempty, strfind(cpfield, 'data'))));
                for ci = 1:length(cpfield)
                    data.(cpfield{ci}) = alldata.(cpfield{ci});
                end
                
                %patient data fields
                pafields = fds(~cellfun(@isempty, strfind(fds, 'patient_data')));
                %sort samples and copy data
                sampleidx = alldata.typematrix(rowi, typeidx);
                %normal data
                if matchnormal
                    sampleidx = sampleidx(reordersample);
                    [~, normaltypeidx] = ismember([types{i} '_normal'], alldata.type);
                    if normaltypeidx ~=0 
                        normsampleidx = alldata.typematrix(normrowi, normaltypeidx);
                        if strcmp(types{i}, 'patient')
                            for pai = 1:length(pafields)
                                data.([pafields{pai} '_normal']) = alldata.(pafields{pai})(:, normsampleidx);
                            end
                        else
                            if length(size(alldata.(types{i}))) == 3
                                data.([types{i} '_normal']) = alldata.(types{i})(:, normsampleidx, :);
                            else
                                data.([types{i} '_normal']) = alldata.(types{i})(:, normsampleidx);
                            end
                        end
                    end
                end
                %tumor data
                if strcmp(types{i}, 'patient')
                    for pai = 1:length(pafields)
                        data.(pafields{pai}) = alldata.(pafields{pai})(:, sampleidx);
                    end
                else
                    if length(size(alldata.(types{i}))) == 3
                        data.(types{i}) = alldata.(types{i})(:, sampleidx, :);
                    else
                        data.(types{i}) = alldata.(types{i})(:, sampleidx);
                    end
                end
            end    
            %add batch
            b = TCGA.getBatch();
            [~, i] = ismember(data.sample, b.sample(:,1));
            data.batch = NaN(size(data.sample));
            data.batch(i~=0) = b.batch(i(i~=0));
        end
        
        function s = getCna(fn)
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.cna];
            end
            s = TCGA.readtable(fn, 'cna');
        end
        
        function s = getExp(fn)
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.exp];
            end
            s = TCGA.readtable(fn, 'exp');
            
        end
        
        function s = getRNA(fn)
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.rna];
            end
            s = TCGA.readtable(fn, 'rna');
        end
        
        function s = getMircna(fn)
            s = [];
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.mircna];
            end
        end
        
        function s = getMirexp(fn)
            s = [];
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [ TCGA.dirhead TCGA.mirexp];
            end
        end
                
        function s = getMethylation(fn)
            s = [];
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.met];
            end
        end
        
        function s = getMutation(fn)
            s = [];
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.mut];
            end
            s = TCGA.readtable(fn, 'mut');
        end
                
        function s = readAllClinical()
            s.patient = TCGA.getPatients();
            s.slide = TCGA.getSlides();
        end
        
        function patientinfo = getPatients(fn)
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.patient];
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

            patientinfo.sample = sample.sample;
            patientinfo.sampletype = ones(length(patientinfo.sample),1) * 2;
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
            patientinfo.data_num = str2double( t.text(:, i) )';
            [~,i] = ismember(patientinfo.attr_ref, t.colname);
            patientinfo.data_ref = t.text(:, i)';
            [~,i] = ismember(patientinfo.attr_cat, t.colname);
            patientinfo.data_cat = t.text(:, i)';

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
            patientinfo.data_num(end+1:end+4,:) = str2double(tmp)';
        end

        function slideinfo = getSlides(fn)
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.slide];
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
%             normindex = sample.sampletype == 1;

            slideinfo.sample = [sample.sample sample.samplecode sample.samplevial];
            slideinfo.sampletype = sample.sampletype;
                 
%            slideinfo.sample = uniqueCellRows(slideinfo.sample);
%             slideinfo.sample_normal = [sample.sample(normindex) ...
%                 sample.samplevial(normindex)];

%             slideinfo.sample_normal = uniqueCellRows(slideinfo.sample_normal);

            slideinfo.attr = t.colname;

            %top, bottom, mean
            slideinfo.data = NaN(length(t.colname), length(slideinfo.sample), 3);
%             slideinfo.slide_data_normal = NaN(length(t.colname), length(slideinfo.sample_normal), 3);

            %top, bottom will match participants, samples and vials
            %but mean will only match particpants and sampletype (tumor or normal)

            sampfd = {''};%, '_normal'};
            dimfilter = [strcmpi(loc, 'top'), strcmpi(loc, 'bottom'), true(size(loc))];

            for sampi = 1:length(sampfd)
                smp = ['sample' sampfd{sampi}];
                d = ['data' sampfd{sampi}];
                for dim = 1:3
                    for i = 1:length(slideinfo.(smp))                        
                        rowi = strcmp(sample.sample, slideinfo.(smp){i,1}) & ...
                            strcmp(sample.samplecode, slideinfo.(smp){i,2}) & ...
                            strcmp(sample.samplevial, slideinfo.(smp){i,3}) & dimfilter(:,dim);
                        if any(rowi)
                            slideinfo.(d)(:,i,dim) = nanmean( str2double(t.text(rowi,:)), 1);
                        end
                    end
                end
            end    

        end
        
        function batch = getBatch(fn)
            if nargin < 1
                fn = '';
            end
            if isempty(fn)
                fn = [TCGA.dirhead TCGA.batch];
            end
            t = parseText(fn, 'nrowname',1,'ncolname',1,'numeric',true);
            samp = tcgaSampleDecoder.decode(t.rowname);
            batch.sample = [samp.sample samp.samplecode samp.samplevial];
            batch.batch = t.text;
        end
        
        function matchdata = getSubtype(matchdata, fns)
           if nargin < 2
               fns = TCGA.subtype;
           end
           if isempty(fns)
               fns = TCGA.subtype;
           end
           if ~isfield(matchdata, 'subtype')
               matchdata.subtype = {};
               matchdata.subtypemtx = [];
           end
           nf = length(fns);
           ns = length(matchdata.sample);
           for fi = 1:nf
               f = fopen([TCGA.dirhead fns{fi}]);
               t = textscan(f, '%s');               
               fclose(f);
               t = t{1};
               t(strcmpi(t, 'sampleid')) = [];
               fd = strrep(strrep(strrep(fns{fi}, '.txt', ''), ...
                   'SampleId', ''), ...
                   'BreastAug10', '');
                   
               tmp = tcgaSampleDecoder.decode(t);
               if ~ismember(fd, matchdata.subtype)
                   matchdata.subtype{end+1,1} = fd;
                   matchdata.subtypemtx(end+1,:) = zeros(1, ns);
               end
               [~, rowi] = ismember(fd, matchdata.subtype);
               [~, coli] = ismember(tmp.sample, matchdata.sample);
               coli(coli==0) = [];
               %in case this call is to update the existing column
               matchdata.subtypemtx(rowi, :) = 0; 
               matchdata.subtypemtx(rowi, coli) = 1;
           end
           
        end
    end

    
    methods (Access = private, Static)
        function tb = readtable(fn, datatype)
            t = parseText(fn, 'ncolname', 1, 'nrowname', 1, 'numeric', true);
            if ismember(datatype, {'cna', 'exp', 'mut'})
                tb.attr = str2double(t.rowname);
                tb.sample = t.colname(2:end);
                tb.data = t.text;
                if strcmp(datatype, 'mut')
                    tb.data = sparse(tb.data);
                end
            elseif ismember(datatype, {'mircna', 'mirexp'})
                tb.attr = t.rowname;
                tb.sample = t.colname(2:end);
                tb.data = t.text;
            elseif ismember(datatype, {'rna'})
                %RNA-seq data, take log2
                tb.attr = str2double(t.rowname);
                tb.sample = t.colname(2:end);
                tb.data = log2(t.text + 1);
            else
                error('unknown datatype %s', datatype);
            end
                        
            
            %check if there are replicates in the samples or genes
            if length(unique(tb.sample)) ~= length(tb.sample)
                fprintf('%s, replicated samples; take mean\n', datatype);
                [count uid] = eleCounts(tb.sample);
                uid(count<2) = [];
                remove = [];
                for ui = 1:length(uid)
                    coli = find(strcmp(tb.sample, uid{ui}));
                    tb.data(:, coli(1)) = nanmean(tb.data(:, coli), 2);
                    remove = [remove; coli(2:end)];
                end
                tb.sample(remove) = [];
                tb.data(:, remove) = [];
            end
            
            if length(unique(tb.attr)) ~= length(tb.attr)
                fprintf('%s, replicated data (rows); take mean\n', datatype);
                [count uid] = eleCounts(tb.attr);
                uid(count<2) = [];
                remove = [];
                for ui = 1:length(uid)
                    if isnumeric(tb.attr)
                        rowi = find(tb.attr == uid(ui));
                    else
                        rowi = find(strcmp(tb.attr, uid{ui}));
                    end
                    tb.data(rowi(1), :) = nanmean(tb.data(rowi, :), 1);
                    remove = [remove; rowi(2:end)];
                end
                tb.attr(remove) = [];
                tb.data(remove, :) = [];
            end
            
            sample = tcgaSampleDecoder.decode(tb.sample);

            tb.sample = [sample.sample sample.samplecode sample.samplevial];
            tb.sampletype = sample.sampletype;
%             if any(tb.sampletype(:,2) == 1) %normal
%                 normcol = tb.sampletype(:,2) == 1;
%                 tb.sample_normal = tb.sample(normcol);
%                 tb.([datatype '_normal']) = tb.data(:, normcol);
%                 tb.sample(normcol) = [];
%                 tb.data(:, normcol) = [];
%             end
        end
    end
    
end

