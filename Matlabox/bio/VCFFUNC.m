classdef VCFFUNC < handle
    
    methods (Static)
        function [output, status] = filterByCoverage(vcffn, bedfn, varargin)            
            %both vcf and bed should be sorted
            para = VCFFUNC.processpara(varargin{:});
            if ~isempty(para.returncol)
                columns = sprintf('%d,', para.returncol);
                select = sprintf(' | cut -f%s"', columns(1:end-1));
            else
                select = '"';
            end
            cmd = 'ssh node025 "bedops -e -1 ';
            vcf = sprintf(' <(cat %s | /data/NYGC/Software/python/bin/python2.7 /data/NYGC/Software/bedops/bin/vcf2bed -)', vcffn);
            bed = sprintf(' <( gawk ''\\$4>%d {print \\$0}'' %s )', para.mincov-1, bedfn);
            %vcf = sprintf('-a <(gawk ''{OFS="\\t"; if($1 ~ "^chr"){sub("chr","",$1); if($1 == "M"){$1="MT"}};  print $0}'' %s ) ', vcffn);
            %bed = sprintf('-b <(gawk ''{OFS="\\t"; if($4 >=%d) {if($1 ~ "^chr"){sub("chr","",$1); if($1 == "M"){$1="MT"}};  print $0}}'' %s ) ',para.mincov, covbed);
            fprintf([cmd, vcf, strrep(bed, '\', '\\'), select, '\n']);
            [status, output] = system([cmd, vcf, bed, select]);            
            if status == 0
                output = VCFFUNC.parseOutput(output, para.parseformatstr);
            end
        end
        function [output, status] = extract(vcffn, varargin)
            para = VCFFUNC.processpara(varargin{:});
            inputstring = '';
            if ~isempty(para.format)                
                inputstring = sprintf('%s -v formatname=%s', inputstring, strjoin(para.format, ':'));
            end
            if ~isempty(para.returncol)
                inputstring = sprintf('%s -v col=%s', inputstring, strjoin( num2cellstr( para.returncol), ':'));                
            end
            [status, output] = system(['gawk' inputstring ' -f ' para.awkprog ' ' vcffn]);
            if status == 0
                output = VCFFUNC.parseOutput(output, para.parseformatstr);
            end
        end
        function [output, status] = getHeaderLine(vcffn)
            [status, output] = system(sprintf('egrep ^#CHROM %s',vcffn));
            if status ~= 0
                output = {};
            else
                output = strtrim(output);
                output = textscan(output, '%s');
                output = output{1};
            end
        end
        
        function parsedoutput = parseOutput(output, formatstr)
            output = strtrim(output);
            parsedoutput = textscan(output, formatstr);            
        end
        
        function para = processpara(varargin)
            para.mincov = 1;            
            para.returncol = [1, 2];
            para.format = {};
            para.info = {};
            para.parseformatstr = '%s %f';
            para.awkprog = 'shscripts/vcf.getdata.awk';
            para = assignpara(para, varargin{:});            
            if size(para.returncol,1) > size(para.returncol,2)
                para.returncol = para.returncol';
            end
            if size(para.format,1) > size(para.format,2)
                para.format = para.format';
            end
            if size(para.info,1) > size(para.info,2)
                para.info = para.info';
            end
        end
        
    end
end