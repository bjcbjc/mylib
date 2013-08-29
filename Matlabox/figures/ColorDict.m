classdef ColorDict < handle
    properties
        key;
        cmap;        
    end
    
    methods
        function obj = ColorDict(color, key)
            if nargin < 2
                key = {};
            end
            if iscellstr(color)
                obj.cmap = ColorDict.translate(color);
            else
                obj.cmap = color;
            end
            
            obj.key = key(:);    
        end
        
        function add(obj, addcolor, addkey)
            if nargin < 3
                if ~isempty(obj.key)
                    error('require name/key for the color');
                end
                addkey = {};
            end
            if iscellstr(addcolor)
                obj.cmap = [obj.cmap; ColorDict.translate(addcolor)];
            else
                obj.cmap = [obj.cmap; addcolor];
            end
            if ~isempty(addkey)
                obj.key = [obj.key; addkey(:)];
            end
        end
        
        function clr = color(obj, name)
            [tf, i] = ismember(name, obj.key);
            if tf
                clr = obj.cmap(i,:);
            else
                error('unknown name/key %s to retrieve color', name);
            end
        end
        
        function k = keys(obj)
            k = obj.key;
        end
        
        function [tf, clr] = hasColor(obj, querykey)
            [tf, i] = ismember(querykey, obj.key);
            clr = obj.cmap(i,:);
        end
    end
    
    methods (Static = true) 
        function color = translate(colorchar)
            preset = {'r', [0.5, 0.1, 0.1]; ...
                'g', [0.1, 0.1, 0.6]; ...
                'b', [0.02, 0.19, 0.4]; ...
                'y', [0.8, 0.8, 0.1]; ...
                'm', [0.8, 0.1, 0.8]; ...
                'c', [0.1, 0.8, 0.8]; ...
                'w', [1, 1, 1]; ...
                'k', [0, 0, 0];
                'o', [251 170 39]/255; ...
                'p', [0.5 0.1 0.5];
                'n', [0.7, 0.7, 0.7]};
            [~, i] = ismember(colorchar, preset(:,1));
            if any(i==0)
                u = colorchar(i==0);
                error('unknown color %s', strcat(u{:}));
            end
            color = cell2mat(preset(i,2));
        end
    end
end