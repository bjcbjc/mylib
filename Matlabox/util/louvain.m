function [c, Q] = louvain(adjmtx, varargin)

    para.adjfile = 'louvain.tmp.txt';
    para.progpath = '/Users/bjc/Lab/Tools/Programs/Community_latest/';
    
    para = assignpara(para, varargin{:});
    
    [i, j] = find(adjmtx);
    
    dlmwrite(para.adjfile, [i-1 j-1], '\t');
    
    % convert txt file to bin file    
    fprintf('converting to bin\n');
    command = [ para.progpath 'convert -i ' para.adjfile ' -o ' para.adjfile '.bin' ];    
    system( command );
    
    % run louvain algorithm    
    fprintf('running louvain\n');
    command = [ para.progpath 'community ' para.adjfile '.bin -l -1 -v > ' para.adjfile '.graph.tree' ];    
    [~, r] = system( command );
    
    % find each iteration's modularity    
    Q = find_modularity( r )
    
    % find nu. of levels
    command = [ para.progpath 'hierarchy ' para.adjfile '.graph.tree' ];    
    [~, r] = system( command );    
    
    r = strtok( r, 10 );
    r = regexprep( r, 'Number of levels: ', '' );
    nu_levels = str2double( r )-1;
        
    fprintf('max level is %d\n', nu_levels );
    
    c = cell(nu_levels, 1);
    
    % import each of the levels g.t. 0
    for level = 1:nu_levels
        fprintf( 1, 'MATLAB: importing level %d\n', level );
        
        command = [ para.progpath 'hierarchy ' para.adjfile '.graph.tree -l ' num2str( level ) ' > ' para.adjfile '.tmp' ];        
        system( command );
        
        hierarchy_output = load( [ para.adjfile '.tmp' ] );
        
        c{level} = hierarchy_output( :, 2 ) + 1;
    end
end


function Q = find_modularity( r )
    % Q = find_modularity( r )
    %
    % convert the text output into modularity score of each iteration
    signature = '  modularity increased from %f to %f';
    idx = 0;

    while( ~isempty( r ) )
        % read a line and match it to the signature
        [token, r] = strtok( r, char( 10 ) );
        a = sscanf( token, signature );

        if( ~isempty( a ) )
            % signature matched copy new modularity
            idx = idx + 1;
            Q( idx ) = a( 2 );
        end
    end
end






