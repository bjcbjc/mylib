function NCBI = getLatestNCBI(perlscriptDir, varargin)
%update NCBI.mat
% example: newNCBI = getLatestNCBI('/nethome/bjchen/DATA/GenomicInfo/NCBI/','sshCmd','ssh node001', 'savemat', true);
    if perlscriptDir(end) ~= '/'
        perlscriptDir = [perlscriptDir '/'];
    end

    para.basedir = '';
    para.backupdir = '';
    para.cmdGetGeneHistory = '';
    para.sshCmd = '';
    para.savemat = false;
    
    para = assignpara(para, varargin{:});
    
    if isempty(para.basedir), para.basedir = perlscriptDir; end
    if isempty(para.backupdir), para.backupdir = [perlscriptDir, 'backup/']; end
    
    if isempty(para.cmdGetGeneHistory)
        tmpstr = sprintf('cd %s; wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz; gunzip -f gene_history.gz', para.basedir);
        if isempty(para.sshCmd)
            para.cmdGetGeneHistory = tmpstr;
        else
            para.cmdGetGeneHistory = sprintf('%s "%s"', para.sshCmd,tmpstr);
        end
    end
    
    if para.basedir(end) ~= '/'
        para.basedir = [para.basedir '/'];        
    end
    if para.backupdir(end) ~= '/'
        para.backupdir = [para.backupdir '/'];
    end
    
    if ~exist(para.backupdir, 'dir')
        system(sprintf('mkdir %s', para.backupdir));
    end
    
    backupfile = {'gene_history.gz', 'Human_Exons.txt', 'Human_Gene_Info_All.txt', ...
        'Human_Gene_Info.txt', 'Human_nonCoding_Gene_Info.txt'};
    for i = 1:length(backupfile)
        status = system(sprintf('mv %s%s %s',para.basedir, backupfile{i}, para.backupdir) );
        if status ~= 0            
            if i > 1
                recover(para.basedir, para.backupdir, backupfile(1:i-1));
            end
            error('error moving file %s',backupfile{i});            
        end
    end
    
    %call UD's script to download and preprocess data
    if isempty(para.sshCmd)
        status = system([sprintf('cd %s; ', para.basedir), perlscriptDir 'Gene_Info.pl']);
    else
        status = system(sprintf('%s "%s"',para.sshCmd, [sprintf('cd %s; ', para.basedir), perlscriptDir 'Gene_Info.pl']));
    end
    if status ~= 0
        recover(para.basedir, para.backupdir, backupfile);
        error('error calling perlscript');
    end
    %download gene_history.gz
    status = system(para.cmdGetGeneHistory);
    if status ~= 0
        recover(para.basedir, para.backupdir, backupfile);
        error('error downloading gene_history.gz');
    end
    
    NCBI = buildGenomeInfo( [para.basedir, 'Human_Gene_Info_All.txt'], ...
        [para.basedir, 'gene_history'] );    
    
    if para.savemat
        if exist(sprintf('%sNCBI.mat',para.basedir), 'file')
            system(sprintf('mv %sNCBI.mat %s', para.basedir, para.backupdir));
        end
        save(sprintf('%sNCBI.mat',para.basedir), 'NCBI');
    end
end

function recover(basedir, backupdir, files)
    for i = 1:length(files)
        system(sprintf('mv %s%s %s',backupdir, files{i}, basedir));        
    end
end
        