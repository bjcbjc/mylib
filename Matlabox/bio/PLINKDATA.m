classdef PLINKDATA < handle
    
    
    methods (Static)
        function [dataTable] = readGenome(fn)            
            f = fopen(fn);            
            line = fgetl(f);
            t = textscan(line, '%s');
            header = strrep(t{1}, '-', '_');
            t = textscan(f, ['%s %s %s %s %s ', repmat('%f ', 1, length(header)-5)], ...
                'treatasempty', {'NA'});
            fclose(f);            
            for i = 1:length(header)
                dataTable.(header{i}) = t{i};
            end
        end
        function dataTable = readLD(fn)
            f = fopen(fn);
            
            line = fgetl(f);
            t = textscan(line, '%s');
            header = t{1};            
            t = textscan(f, '%s %f %s %s %f %s %f', ...
                'treatasempty', {'NA'});            
            fclose(f); 
            for i = 1:length(header)
                dataTable.(header{i}) = t{i};
            end            
        end
    end
    
end