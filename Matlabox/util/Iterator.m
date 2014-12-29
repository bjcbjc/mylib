classdef Iterator < handle
    properties
        data;
        indice;
        n;
        total;
        anchorIdx;
        curIdx;
        done;
        N;
    end
    methods
        function obj = Iterator(data)
            obj.data = data;
            obj.indice = ones(1, length(obj.data));
            obj.n = cellfun(@length, obj.data);
            obj.total = prod(obj.n);
            obj.anchorIdx = 1;
            obj.curIdx = 1;
            obj.done = false;
            obj.N = length(obj.data);
        end
        function item = getNext(obj)
            if obj.done
                item = {};
            else
                item = cell(length(obj.data), 1);
                for i = 1:length(obj.data)                    
                    item{i} = obj.data{i}{obj.indice(i)};
                end
                obj.updateIndice();
            end
        end
        
        function updateIndice(obj)            
            if obj.indice(obj.curIdx) == obj.n(obj.curIdx) 
                loopIdx = obj.curIdx + 1;
                while 1                    
                    if obj.n(loopIdx) > obj.indice(loopIdx)
                        obj.indice(loopIdx) = obj.indice(loopIdx) + 1;
                        obj.indice(1:loopIdx-1) = 1;
                        obj.curIdx = 1;
                        break
                    end
                    loopIdx = loopIdx + 1;
                    if loopIdx > obj.N
                        obj.done = true;
                        break;
                    end                    
                end
            else            
                obj.indice(obj.curIdx) = obj.indice(obj.curIdx) + 1;
            end
%             disp(obj.indice);
        end
    end
end
