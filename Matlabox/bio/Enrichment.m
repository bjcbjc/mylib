classdef Enrichment < handle
    
    
    methods (Static)
        function res = enrichmentTest(ANT, testGene, background, varargin)
            %ANT: struct, contain gene, name, mtx
            %testGene: vector of gene, or cell ({gene_list, genelist}), or
            %struct (testGene.gene, testGene.mtx (binary) )
            %background: vector of gene
            %
            para.numAntMin = 30;
            para.numAntMax = 500;
            if nargin < 3, background = []; end
            para = assignpara(para, varargin{:});
            
            ANT.gene = double(ANT.gene);
            if isempty(background)
                background = ANT.gene;
            end
            ANT = Enrichment.filterAnnotation(ANT, background, para);
                        
            nSet = length(ANT.name);
            nBackgroundGene = length(background);
%             nBackgroundGene = length(ANT.gene);
            nAnnotatedGene = sum(ANT.mtx~=0, 1)';
            
            if isstruct(testGene)
                testMtx = false(length(ANT.gene), size(testGene.mtx,2));
                [~, i] = ismember(testGene.gene, ANT.gene);
                testMtx(i(i~=0), :) = testGene.mtx(i~=0, :);
            elseif iscell(testGene)
                if ~iscell(testGene{1}) %multiple lists of genes
                    testMtx = false(length(ANT.gene), length(testGene));
                    for i = 1:length(testGene)
                        testMtx(:, i) = ismember(ANT.gene, testGene{i});
                    end
                else
                    testMtx = ismember(ANT.gene, testGene); %single list of gene (names)
                end
            else
                testMtx = ismember(ANT.gene, testGene); %single list of gene (id)
            end
            
            nTestGene = sum(testMtx, 1)';
            nTestSet = size(testMtx, 2);
            pval = NaN(nSet, nTestSet);
            nTestAnt = NaN(size(pval));
            if nTestSet < nSet
                for i = 1:nTestSet
                    nTestAnt(:, i) = sum(bsxfun(@and, testMtx(:,i), ANT.mtx), 1);
%                     pval(:, i) = ...
%                         sum(hygecdf( nTestAnt(:, i):min(nTestGene(i), nAnnotatedGene) , ...
%                         nBackgroundGene, nTestGene(i), nAnnotatedGene));
                    pval(:, i) = ...
                        hygecdf( nTestAnt(:, i)-1 , ...
                        nBackgroundGene, nTestGene(i), nAnnotatedGene, 'upper');
                end
            else
                for i = 1:nSet
                    nTestAnt(i, :) = sum(bsxfun(@and, testMtx, ANT.mtx(:,i)), 1);
                    pval(i, :) = ...
                        hygecdf( nTestAnt , ...
                        nBackgroundGene, nTestGene, nAnnotatedGene(i), 'upper');
                end
            end        
            res.name = ANT.name;
            res.nBackgroundGene = nBackgroundGene;
            res.nAnnotatedGene = nAnnotatedGene;
            res.nAnnotatedTestGene = nTestAnt;
            res.nTestGene = nTestGene;
            res.pval = pval;
        end
        
        function ANT = filterAnnotation(ANT, background, para)
            validGene = ismember(ANT.gene, background);
            ANT.gene = ANT.gene(validGene);
            ANT.mtx = ANT.mtx(validGene, :);
            
            nAnnotatedGene = sum(ANT.mtx~=0, 1);
            validSet = nAnnotatedGene <= para.numAntMax & nAnnotatedGene >= para.numAntMin;
            ANT.name = ANT.name(validSet);
            ANT.mtx = ANT.mtx(:, validSet);
        end
        
        
        
    end
    
    
end