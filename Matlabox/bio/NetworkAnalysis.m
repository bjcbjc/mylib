classdef NetworkAnalysis < handle
    
    methods (Static)
        function ws = weightedScore(gene, score, weight, network)
            % gene, score, weight: #gene x 1
            % network:
            %   network.gene
            %   network.mtx: network.gene x network.gene
            %   if network is annt (annt.mtx: annt.gene x annt.name), it
            %   will be converted to network
            %
            score = score(:)';
            weight = weight(:)'; %weight need to be row for bsxfun
            valid = (~isnan(score) & ~isnan(weight))';
%             ws = score;
            ws = NaN(size(score));
            network.gene = double(network.gene);
            if isfield(network, 'name')
                network = NetworkAnalysis.annt2net(network);
            end
            
            %if orphan genes, don't get a score
            [~, validSourceIdx] = ismember(gene, network.gene);
            validSource = ismember(gene, network.gene(any(network.mtx,2)));
            [validTarget, validTargetIdx] = ismember(gene, network.gene);            
            validSource = validSource & valid;            
            validTarget = validTarget & valid;
            w = bsxfun(@times, weight(validTarget), network.mtx(validSourceIdx(validSource), validTargetIdx(validTarget)));
            ws(validSource) = sum( bsxfun(@times, w, score(validTarget)), 2) ./ sum(w,2);
            
%             for i = 1:length(validIdx)
%                 anntIdx = network.mtx(geneIdx(validIdx(i)),:) ~= 0;
%                 networkGene = network.gene( any(network.mtx(:, anntIdx), 2) );
%                 idx = ismember(gene, networkGene);
%                 ws(validIdx(i)) = sum(weight(idx).*score(idx)) / sum(weight(idx));
%             end
        end
        
        function res = neighborStats(gene, score, network, varargin)
            % gene, score: #gene x sample
            % network:
            %   network.gene
            %   network.mtx: network.gene x network.gene
            %   if network is annt (annt.mtx: annt.gene x annt.name), it
            %   will be converted to network
            %
            para.minTop = 10; %minimum neighbor for mean top 50% (that is, if n < 20, take top 10)
            para.topPct = 0.5;
            
            para = assignpara(para, varargin{:});
             
            if nnz(isnan(score)) > 0
                error('score has NaN');
            end
            nSample = size(score,2);
            %stats: [num neighbor, std, mean (abs), #score>=1, #score<=-1,
            %#score>=2, #score<=-2, mean top 50%]
            res.stats = NaN(length(gene), 8, nSample);
            res.name = {'numNeighbor', 'sd', 'meanAbs', 'numPos1', 'numNeg1', ...
                'numPos2', 'numNeg2', 'meanAbsTopHalf'};
            network.gene = double(network.gene);
            if isfield(network, 'name')
                network = NetworkAnalysis.annt2net(network);
            end            
            
            [~, scoreGeneIdx, networkGeneIdx] = intersect(gene, network.gene);            
            score = score(scoreGeneIdx, :);
            network.mtx = network.mtx(networkGeneIdx, :);
            network.mtx = network.mtx(:, networkGeneIdx);
            network.gene = network.gene(networkGeneIdx);
            
            res.stats(:, 1) = 0;
            for i = 1:length(scoreGeneIdx)
                if ~any(network.mtx(i, :)), continue; end
                neighbor = network.mtx(i,:);
                nn = sum(neighbor);
                res.stats(scoreGeneIdx(i), 1, :) = nn;
                res.stats(scoreGeneIdx(i), 3, :) = mean(abs(score(neighbor,:)),1);
                res.stats(scoreGeneIdx(i), 4, :) = sum(score(neighbor,:) >= 1, 1);
                res.stats(scoreGeneIdx(i), 5, :) = sum(score(neighbor,:) <= -1, 1);
                res.stats(scoreGeneIdx(i), 6, :) = sum(score(neighbor,:) >= 2, 1);
                res.stats(scoreGeneIdx(i), 7, :) = sum(score(neighbor,:) <= -2, 1);
                
                if nn >= 5
                    res.stats(scoreGeneIdx(i), 2, :) = std(score(neighbor,:), 0, 1);
                    ntop = min(nn, max(para.minTop, floor(nn*para.topPct)));
                    mtx = sort(abs(score(neighbor,:)), 1, 'descend');                    
                    res.stats(scoreGeneIdx(i), 8, :) = mean(mtx(1:ntop,:),1);
                end
            end                        
        end
        
        
        function ds = neightborScore(gene, score, network, querygene)
            % get neighbor scores for the queried genes
            % score: gene x sample
            % gene: gene x 1
            % querygene: query x 1
            % network: .gene, .mtx: .gene x .gene
            %
            % ds: query x 1 cell, each cell is downstream x sample
            %
            [~, gi, gj] = intersect(gene, network.gene);
            gene = gene(gi);
            score = score(gi,:);
            network.gene = network.gene(gj);
            network.mtx = network.mtx(gj,gj);
            
            n = length(querygene);
            ds = cell(n, 1);
            [~, idx] = ismember(querygene, gene);
            for i = 1:n
                if idx(i) == 0, continue; end
                ds{i} = score(network.mtx(idx(i),:), :);
            end
        end
        
        function network = annt2net(annt)
            % annt.gene, annt.name
            % annt.mtx: annt.gene x annt.name
            % return:
            %   network.gene
            %   network.mtx: network.gene x network.gene
            %
            [ng, na] = size(annt.mtx);            
            network.gene = annt.gene;
            network.mtx = false(ng, ng);
            if na < ng
                for i = 1:na
                    gi = annt.mtx(:,i);
                    network.mtx(gi, gi) = true;
                end
            else
                for i = 1:ng
                    ai = annt.mtx(i,:);
                    gi = any(annt.mtx(:, ai), 2);
                    network.mtx(i, gi) = true;
                    network.mtx(gi, i) = true;
                end
            end
        end
        
        function [reachedMtxCombined, reachedMtxList] = ...
                downstreamAcrossNetwork(networkMtxList, varargin)
            % for each srcIdx (row index), traverse each network (assumed
            % directed, row -> col) to get all downstream nodes, including
            % feedback loops; then all downstream nodes are aggregated
            % across multiple networks in the list
            %
            para.directed = true;
            para.loopToSelf = true; %include self if reachable via loop
            para.sparse = false;
            para = assignpara(para, varargin{:});
            
            paracopy = para;
            paracopy.sparse = false;
            paracopy = paraPairArray(paracopy);
            ng = size(networkMtxList{1},1);            
            nNetwork = length(networkMtxList);
            reachedMtxCombined = false(ng, ng);
            if nargout > 1                
                reachedMtxList = cell(nNetwork, 1);
                for i = 1:nNetwork                    
                    reachedMtxList{i} = NetworkAnalysis.downstream( ...
                        networkMtxList{i}, paracopy{:} );
                    reachedMtxCombined = reachedMtxCombined | reachedMtxList{i};
                    if para.sparse && nnz(reachedMtxList{i}) < 0.5*ng*ng
                        reachedMtxList{i} = sparse(reachedMtxList{i});
                    end
                end
            else                
                for i = 1:nNetwork
                    reachedMtxList = NetworkAnalysis.downstream( ...
                        networkMtxList{i}, paracopy{:} );
                    reachedMtxCombined = reachedMtxCombined | reachedMtxList;                    
                end
            end
        end
        
        function reachedMtx = downstream2(networkMtx, varargin)
            % for each srcIdx (row index), traverse the network (assumed
            % directed, row -> col) to get all downstream nodes, including
            % feedback loops
            %
            para.directed = true;
            para.loopToSelf = true; %include self if reachable via loop
            para.sparse = false;
            para = assignpara(para, varargin{:});
            
            n = size(networkMtx, 1);            
            if nnz(networkMtx) == 0
                reachedMtx = sparse(n, n);
                return;
            end
            reachedMtx = all_shortest_paths(double(networkMtx));
            reachedMtx = reachedMtx~=0 & ~isinf(reachedMtx);
            if para.loopToSelf
                [i, j] = find(reachedMtx & reachedMtx');
                reachedMtx( sub2ind([n, n], i,i) ) = true;
            end
            if nnz(reachedMtx) < 0.5*n*n && para.sparse
                reachedMtx = sparse(reachedMtx);
            end
        end
        
        function reachedMtx = downstream(networkMtx, varargin)
            % for each srcIdx (row index), traverse the network (assumed
            % directed, row -> col) to get all downstream nodes, including
            % feedback loops
            %
            para.directed = true;
            para.loopToSelf = true; %include self if reachable via loop
            para.sparse = false;
            para = assignpara(para, varargin{:});
            
            n = size(networkMtx, 1);            
            
            if nnz(networkMtx) == 0
                reachedMtx = sparse(false(n, n));
                return;
            end
            
            nEdge = nnz(networkMtx);
            if nEdge > 6000
                reachedMtx = all_shortest_paths(double(networkMtx));
                reachedMtx = reachedMtx~=0 & ~isinf(reachedMtx);
            else            
                reachedMtx = false(n, n);                            
                for i = 1:n
                    %the source node is always returned as the first node, but
                    %we don't know if there is a loop back to the source
                    if ~any(networkMtx(i,:)), continue; end
                    reachedIdx = graphtraverse(networkMtx, i, 'directed', para.directed);
                    reachedMtx(i, reachedIdx(2:end)) = true;
                end
            end
            if para.loopToSelf
                [i, j] = find(reachedMtx & reachedMtx');
                reachedMtx( sub2ind([n, n], i,i) ) = true;
            end
            if nnz(reachedMtx) < 0.5*n*n && para.sparse
                reachedMtx = sparse(reachedMtx);
            end
        end
        
        function network = getNetwork(databaseName, varargin)
            para.path = '/nethome/bjchen/DATA/pathway/';
            para.convertToEntrez = true;
            
            available = {'biocarta', 'kegg', 'nci', 'reactome', 'spike'};
            
            para = assignpara(para, varargin{:});
            databaseName = lower(databaseName);
            
            if ~ismember(databaseName, available)
                error('unknown database %s', databaseName);
            end
            if para.convertToEntrez
                fnname = [para.path 'Network.' databaseName '.entrez.mat'];
            else
                fnname = [para.path 'Network.' databaseName '.mat'];
            end
            if exist(fnname, 'file')
                network = loadStructData(fnname);
            else
                network = NetworkAnalysis.readNetwork([para.path, databaseName, '.txt'], para.convertToEntrez, para.path);
            end
            
        end
        function network = readNetwork(fn, useEntrez, entrezPath)
            if nargin < 2, useEntrez = false; end            
            
            t = parseText(fn, 'nrowname', 0, 'ncolname', 0);
            %name, source, dest, directed, edge-type
            a = regexp(t.text(:,2:3), '\w+:', 'match', 'once');
            a = unique(a(:));
            if length(a) == 1
                t.text(:,2:3) = strrep(t.text(:,2:3), a{1}, '');
            end
            network.identifier = strrep(a{1}, ':', '');
            
            if useEntrez && ~strcmpi(network.identifier, 'entrezgene')
                IDMap = loadStructData([entrezPath, 'UniProtIDMap.mat']);
                invalid = false(size(t.text,1), 1);
                for colIdx = 2:3
                    entrez = GENOMEFUNC.idQuery(IDMap, t.text(:,colIdx), 'Entrez',  'UniProtKB_AC');
                    invalid = invalid | isnan(str2double(entrez));                    
                    t.text(:,colIdx) = entrez;
                end
                %remove nodes that cannot be uniquelly mapped to entrez
                t.text(invalid,:) = [];
                network.identifier = 'entrezgene';
            end
            
            [network.flat.pathwayName, ~, network.flat.pathwayId] = ...
                unique(t.text(:,1));
            network.flat.src = t.text(:,2);
            network.flat.dst = t.text(:,3);
            network.flat.directed = strcmpi(t.text(:,4), 'directed');
            network.flat.edgeType = t.text(:,5);
            if strcmpi(network.identifier, 'entrezgene')
                network.flat.src = str2double(network.flat.src);
                network.flat.dst = str2double(network.flat.dst);
            end
            network.pathwayName = network.flat.pathwayName;
            np = length(network.pathwayName);
            network.pathwayMtxDirected = cell(np, 1);
            network.pathwayMtxUndirected = cell(np, 1);
            network.gene = union(network.flat.src, network.flat.dst);
            ng = length(network.gene);
            [~, rowIdx] = ismember(network.flat.src, network.gene);
            [~, colIdx] = ismember(network.flat.dst, network.gene);
            %there are duplicate edges (eg. different effects etc)
            %hence, contruct double sparse and then convert to logical
            for pIdx = 1:np
                edgeIdx = network.flat.pathwayId == pIdx & network.flat.directed;
                network.pathwayMtxDirected{pIdx} = logical(sparse( ...
                    rowIdx(edgeIdx), colIdx(edgeIdx), 1, ng, ng)); 
                edgeIdx = network.flat.pathwayId == pIdx;                
                network.pathwayMtxUndirected{pIdx} = logical(sparse( ...
                    [rowIdx(edgeIdx); colIdx(edgeIdx)], ...
                    [colIdx(edgeIdx); rowIdx(edgeIdx)], 1, ng, ng)); %symmetric
            end
        end        
    end
end