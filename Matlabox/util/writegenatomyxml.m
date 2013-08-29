function text = writegenatomyxml(structarray)
    %structarray: structure array containing structure to write as
    %   genatomy filters
    %Fields of each struct:
    %   .name: module name
    %   .objtype: {'simple', 'modnet', 'reg', 'modnet, 'splitreg'}
    %   .Set: cell array; set of rows in each module
    %   .Samples: cell array; set of columns in each module
    %   .regprog
    %       .regulators: cell array
    %       .coefficients: array; #regulators or #regulators + 1 (const)
    %       .regulatorinformation: struct array; each struct has 'name', array
    %           of 'regulator' and 'value'
    %   .split
    %       .type    
    %       .splitdata_regulator
    %       .splitdata_type
    %       .splitdata_point
    %       .left/right: split or regprog
    %   
    %
    %
    
    text = '';
    
    nmodule = length(structarray);
    fds = fieldnames(structarray(1));
    for i = 1:nmodule
        text = sprintf('%s<Module Name="%s">\n', text, structarray(i).name);
        
        printfd = {'Set', 'Samples'};
        for pfi = 1:length(printfd)
            if ~ismember(printfd{pfi}, fds), continue; end
            text = sprintf('%s%s',text, printattr(printfd{pfi}, structarray(i).(printfd{pfi})));            
        end
        
        
        if strcmp(structarray(i).objtype, 'reg')
            %regprog
            if ismember('regprog', fds)
                text = sprintf('%s%s', text, printregprog(structarray(i).regprog));
            end
        elseif strcmp(structarray(i).objtype, 'modnet')
            %modnet
            if ~isempty(structarray(i).split)
                text = sprintf('%s%s', text, printmodnet(structarray(i).split));
            end
        elseif strcmp(structarray(i).objtype, 'splitreg')
            %splitprog
            text = sprintf('%s%s', text, printsplitreg(structarray(i).split));
        end
        
        text = sprintf('%s</Module>\n',text);
    end
end

function text = printattr(name, list)
    %print <name>list(1)\tlist(2)...</name>
    text = '';
    if isempty(list), return; end
    text = sprintf('%s<%s>', text, name);
    if iscellstr(list)
        for ei = 1:length(list)
            text = sprintf('%s%s\t', text, list{ei});
        end
    else
        for ei = 1:length(list)
            text = sprintf('%s%g\t', text, list(ei));
        end
    end
    text = sprintf('%s</%s>', text, name);
end

function text = printregprog(regprog)
    fds = fieldnames(regprog);
    printfd = {'regulators', 'coefficients'};
    
    text = '';
    for pfi = 1:length(printfd)
        if ~ismember(printfd{pfi}, fds), continue; end
        pfd = printfd{pfi}; pfd(1) = upper(pfd(1));
        text = sprintf('%s%s', text, printattr(pfd, regprog.(printfd{pfi})));            
    end
    
    if ismember('regulatorinformation', fds)
        if isempty(regprog.regulatorinformation), return; end
        text = sprintf('%s<Regulatorinformation>\n', text);
        nattr = length(regprog.regulatorinformation);
        for ai = 1:nattr
            text = sprintf('%s<Parameter Name="%s">\n', text, regprog.regulatorinformation(ai).name);
            nr = length(regprog.regulatorinformation(ai).regulator);
            if iscellstr(regprog.regulatorinformation(ai).value)
                for ei = 1:nr                
                    text = sprintf('%s<Regulator Name="%s">%s</Regulator>\n', ...
                        text, regprog.regulatorinformation(ai).name{ei}, ...
                        regprog.regulatorinformation(ai).value{ei});
                end
            else
                for ei = 1:nr                
                    text = sprintf('%s<Regulator Name="%s">%g</Regulator>\n', ...
                        text, regprog.regulatorinformation(ai).name{ei}, ...
                        regprog.regulatorinformation(ai).value(ei));
                end
            end
            text = sprintf('%s</Parameter>\n', text);
        end
        text = sprintf('%s</Regulatorinformation>\n', text);
    end
end

function text = printsplitreg(split)
    %recursive
    %   .split    
    %       .type    
    %       .regprog
    %       .splitdata_regulator
    %       .splitdata_type
    %       .splitdata_point    
    %       .left/right: split 
    text = '';
    text = sprintf('%s<Split Type="%s">\n', text, split.type);    
    if strcmpi(split.type, 'lr')
        text = sprintf('%s%s', text, printregprog(split.regprog));
    elseif strcmpi(split.type, 'regressionwithsplit')
        if ~isempty(split.regprog)
            text = sprintf('%s%s', text, printregprog(split.regprog));
        end
        text = sprintf('%s<SplitData Regulator="%s" SplitPoint="%g" Type="%s"/>\n', ...
            text, split.splitdata_regulator, split.splitdata_point, split.splitdata_type);
        
        printfds = {'left', 'right'};
        for pfi = 1:length(printfds)
            if isempty(split.(printfds{pfi})), continue; end
            pfd = printfds{pfi};
            pfd(1) = upper(pfd(1));
            text = sprintf('%s<%s>\n', text, pfd);
%             if strcmpi(split.(printfds{pfi}).objtype, 'regprog')
%                 text = sprintf('%s%s', text, printregprog(split.(printfds{pfi})));
%             else
            text = sprintf('%s%s', text, printsplitreg(split.(printfds{pfi})));
%             end
            text = sprintf('%s</%s>\n', text, pfd);
        end
    end
    
    
    text = sprintf('%s</Split>\n', text);
end

function text = printmodnet(split)
    %recursive
    %   .split    
    %       .type    
    %       .regprog
    %       .splitdata_regulator
    %       .splitdata_type
    %       .splitdata_point
    %       .splitdata_data
    %       .left/right: split 
    text = '';
    text = sprintf('%s<Split Type="%s">\n', text, split.type);    
    if ~strcmpi(split.type, 'oneregulator')
        error('modnet should have split type OneRegulator');
    end
    
    text = sprintf('%s<SplitData Regulator="%s" SplitPoint="%g" Type="%s"/>\n', ...
        text, split.splitdata_regulator, split.splitdata_point, split.splitdata_type);
    
    printfds = {'left', 'right'};
    for pfi = 1:length(printfds)
        if isempty(split.(printfds{pfi})), continue; end
        pfd = printfds{pfi};
        pfd(1) = upper(pfd(1));
        text = sprintf('%s<%s>\n', text, pfd);
        text = sprintf('%s%s', text, printmodnet(split.(printfds{pfi})));
        text = sprintf('%s</%s>\n', text, pfd);
    end
    
    text = sprintf('%s</Split>\n', text);
end