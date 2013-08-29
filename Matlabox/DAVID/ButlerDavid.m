classdef ButlerDavid < handle
    
    properties
        davidfile = '/Users/bjc/Lab/Projects/DATA/DAVID/DAVID.mat';
        DAVID;                
        numAntMin = 15;
        numAntMax = Inf;   
        useSet;
        useGene;
        davidwebservice = [];
        davidwebserviceDefaultCat = {'BBID', 'BIOCARTA', 'COG_ONTOLOGY', ...
            'GOTERM_BP_FAT', 'GOTERM_CC_FAT', 'GOTERM_MF_FAT', 'KEGG_PATHWAY', ...
            'OMIM_DISEASE', 'SP_PIR_KEYWORDS', 'UP_SEQ_FEATURE'};
    end
    
    methods
        function obj = ButlerDavid(varargin)
            if nargin > 0                
                obj = assignpara(obj, varargin{:});
            end  
            obj.loadData();
        end
        
        function delete(obj)
            delete(obj.davidwebservice);
            obj.davidwebservice = [];
        end
        
        function filterSet(obj, varargin)
            para.and = [];
            para.or = [];
            para.numAntMin = [];
            para.numAntMax = [];
            para.dbname = [];
            para.dbidx = [];
            para.set = [];
            
            if nargin == 1 %reset
                obj.useSet(:) = true;
            else
                para = assignpara(para, varargin{:});
                fds = fieldnames(para);
                for i = 1:length(fds)
                    if isempty(para.(fds{i})), continue; end
                    if ismember(fds{i}, {'and', 'or'})
                        func = str2func(fds{i});
                        obj.useSet = func(obj.useSet, para.(fds{i}));
                    elseif ismember(fds{i}, {'dbname'})
                        obj.useSet = obj.useSet & ...
                            ismember(obj.DAVID.dbname(obj.DAVID.dbidx), para.dbname);
                    elseif ismember(fds{i}, {'dbidx', 'set'})
                        obj.useSet = obj.useSet & ...
                            ismember(obj.DAVID.(fds{i}), para.(fds{i}));
                    elseif ismember(fds{i}, {'numAntMin', 'numAntMax'})
                        if strcmp(fds{i}, 'numAntMin')
                            func = @ge;
                        else
                            func = @le;
                        end
                        obj.useSet = obj.useSet & ...
                            func(sum( obj.DAVID.mtx( obj.useGene, :), 1)', para.(fds{i}));
                    end
                end
            end
            
        end
        
        function filterGene(obj, varargin)
            para.and = [];
            para.or = [];
            para.gid = [];
            
            if nargin == 1 %reset
                obj.useGene(:) = true;
            else
                para = assignpara(para, varargin{:});
                fds = fieldnames(para);
                for i = 1:length(fds)
                    if isempty(para.(fds{i})), continue; end
                    if ismember(fds{i}, {'and', 'or'})
                        func = str2func(fds{i});
                        obj.useGene = func(obj.useGene, para.(fds{i}));
                    elseif ismember(fds{i}, {'gid'})
                        obj.useGene = obj.useGene & ...
                            ismember(obj.DAVID.(fds{i}), para.(fds{i}));
                    end
                end
            end
        end
       
        function numGene = numAntGene(obj)
            numGene = sum(obj.DAVID.mtx(obj.useGene, obj.useSet), 1);
        end
        
        function gene = validGene(obj)
            gene = obj.DAVID.gid(obj.useGene);
        end
        
        function [annotation, datasource] = validSet(obj)
            annotation = obj.DAVID.set(obj.useSet);
            datasource = obj.DAVID.dbname(obj.DAVID.dbidx(obj.useSet));
        end
        
        function [pval, numTestAnt, annotation, source] = enrichmentTest(obj, testGene, testSet)
            %testGene: vector of gene id
            %testSet: #gene x #set, binary
            [validTestGeneIdx, validTestSetIdx, testSet] = ...
                obj.setupEnrichTest(testGene, testSet);
            ns = size(testSet, 2);
            ng = length(validTestGeneIdx);            
            pval = NaN(ns, sum(obj.useSet));
            numTestAnt = NaN(size(pval));
            if sum(validTestSetIdx) < sum(obj.useSet)
                validTestSetIdx = find(validTestSetIdx);
                na = obj.numAntGene();
                for i = 1:length(validTestSetIdx)
                    numTestAnt(validTestSetIdx(i), :) = ...
                        sum( ...
                        bsxfun(@and, testSet(:, validTestSetIdx(i)), ...
                        obj.DAVID.mtx(obj.useGene, obj.useSet)), 1);
                    pval(validTestSetIdx(i), :) = ...
                        1 - hygecdf(max(0, numTestAnt(validTestSetIdx(i), :) - 1), ...
                        ng, sum(testSet(:, validTestSetIdx(i)), 1), na);
                end
            else
                runSet = find(obj.useSet);
                na = sum( testSet(:, validTestSetIdx), 1);
                for i = 1:length(runSet)
                    numTestAnt(validTestSetIdx, i) = ...
                        sum( ...
                        bsxfun(@and, testSet(:, validTestSetIdx), ...
                        obj.DAVID.mtx(obj.useGene, runSet(i))), 1);
                    pval(validTestSetIdx, i) = ...
                        1 - hygecdf(max(0, numTestAnt(validTestSetIdx, i) - 1), ...
                        ng, na, sum(obj.DAVID.mtx(obj.useGene, runSet(i)),1) );
                end
            end
            [annotation, source] = obj.validSet();
        end
        
        function [pval, numTestAnt, numTestNonAnt] = gseaTest(obj, testGene, testSet, varargin)
            %testGene: vector of gene id
            %testSet: #gene x #set, continuous scores
            para.score_func = @abs;
            para.P = 1;
            para.num_permutations = 1000;
            if ~isempty(varargin)
                para = assignpara(para, varargin{:});
            end            
            
            [~, validTestSetIdx, testSet] = ...
                obj.setupEnrichTest(testGene, testSet);
            
            pval = NaN( size(testSet,2), length(obj.useSet), 2);
            numTestAnt = NaN( size( pval ) );
            numTestNonAnt = NaN( size( pval ) );
            validTestSetIdx = find(validTestSetIdx);
            for i = 1:length(validTestSetIdx)                
                [p, ~, ~, nAnt, nNonAnt] = ...
                    GSEA(testSet(:, validTestSetIdx(i)), ...
                    obj.DAVID.mtx(obj.useGene, obj.useSet), ...
                    para.score_func, para.P, para.num_permutations);
                pval(validTestSetIdx(i), obj.useSet, :) = p';
                numTestAnt(validTestSetIdx(i), obj.useTest, :) = nAnt';
                numTestNonAnt(validTestSetIdx(i), obj.useTest, :) = nNonAnt';
            end
        end
        
        function setupDavidWebService(obj, backgroundGene, testGene, testCategories)
            if nargin < 4, testCategories = {}; end;
            if nargin < 3, testGene = []; end
            if nargin < 2, backgroundGene = []; end
            if isempty(obj.davidwebservice)
                obj.davidwebservice = DAVIDWebService();
                a = authenticate(obj.davidwebservice, 'bc2252@columbia.edu');
                if ~a
                    error('cannot authenticate account for DAVID service');
                end                
                setCategories(obj.davidwebservice, cellarray2str(obj.davidwebserviceDefaultCat, ','));
            end
            if ~isempty(testCategories)
                setCategories(obj.davidwebservice, cellarray2str(testCategories, ','));            
            end
            if ~isempty(backgroundGene)
                backgroundGene = backgroundGene(:)';
                backgroundGene = mat2str(backgroundGene);
                backgroundGene = strrep(backgroundGene, ' ', ',');
                backgroundGene([1 end]) = [];
                addList(obj.davidwebservice, backgroundGene, 'ENTREZ_GENE_ID', 'background', 1);
            end
            if ~isempty(testGene)
                testGene = testGene(:)';
                testGene = mat2str(testGene);
                testGene = strrep(testGene, ' ', ',');
                testGene([1 end]) = [];
                addList(obj.davidwebservice, testGene, 'ENTREZ_GENE_ID', 'list', 0);
            end
        end
        
        function report = runDavidGeneCluster(obj, backgroundGene, testGene, testCategories, varargin)
            para.overlap = 4;
            para.initGroup = 4;
            para.finalGroup = 4;
            para.linkage = 0.5;
            para.kappa = 20; %lowest stringency, more aggresive merging
            if nargin < 4, testCategories = {}; end;
            if ~isempty(varargin)
                para = assignpara(para, varargin{:});
            end
            
            setupDavidWebService(obj, backgroundGene, testGene, testCategories);
            
            report = getGeneClusterReport(obj.davidwebservice, para.overlap, ...
                para.initGroup, para.finalGroup, para.linkage, para.kappa);
            
        end
        
        function [cluster, score, report] = runDavidTermCluster(obj, backgroundGene, testGene, testCategories, varargin)
            para.overlap = 3;
            para.initGroup = 3;
            para.finalGroup = 3;
            para.linkage = 0.5;
            para.kappa = 20; %lowest stringency
            if nargin < 4, testCategories = {}; end;            
            if ~isempty(varargin)
                para = assignpara(para, varargin{:});
            end
            
            setupDavidWebService(obj, backgroundGene, testGene, testCategories);
            
            report = getTermClusterReport(obj.davidwebservice, para.overlap, ...
                para.initGroup, para.finalGroup, para.linkage, para.kappa);            
            ncluster = length(report);
            cluster = cell(ncluster,1);
            score = NaN(ncluster, 1);
            for i = 1:ncluster
                cluster{i} = report(i).name;
                score(i) = str2double( report(i).score );
            end
        end
        
    end
    
    methods (Access = private)
        function loadData(obj)
            d = load(obj.davidfile);
            fds = fieldnames(d);
            obj.DAVID = d.(fds{1});
            obj.useSet = true( length(obj.DAVID.set), 1);
            obj.useGene = true( length(obj.DAVID.gid), 1);
        end
        
        function [validTestGeneIdx, validTestSetIdx, testSet] = setupEnrichTest(obj, testGene, testSet)
            ng = length(testGene);
            [ng1, ns] = size(testSet);
            if ng ~= ng1
                error('length(testGene) ~= size(testSet,1)');
            end
            if size(testGene, 2) > size(testGene, 1)
                testGene = testGene';
            end
            validAntGene = ismember(obj.DAVID.gid, testGene);            
            if any(validAntGene ~= obj.useGene)
                fprintf('change .useGene index by overlapping testGene\n'); 
                obj.filterGene('and', validAntGene);                                
            end
            fprintf('change .useSet index by numAntMin/Max\n');
            obj.filterSet('numAntMin', obj.numAntMin, 'numAntMax', obj.numAntMax);
            [~, validTestGeneIdx] = ismember(obj.DAVID.gid(obj.useGene), testGene);
            %reorder the data
            testSet = testSet(validTestGeneIdx, :);
            ngInSet = sum(testSet, 1);
            validTestSetIdx = ngInSet >= obj.numAntMin & ngInSet <= obj.numAntMax;
        end        
    end
end