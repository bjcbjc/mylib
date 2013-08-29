function obj = genatomyOBJ(objtype, varargin)
    %generate an structure of genatomy module
    %
    %objtype = {'simple', 'reg', 'modnet', 'splitreg', 'reginfo'}
    %varargin (attribue values, if specified; default otherwise)
    %
    %   .objtype
    %   ----- simple, reg, modnet, splitreg -----
    %   .Set
    %   .Samples
    %   .regprog (reg, splitreg)
    %   .split (modnet, splitreg)   
    %   ------ regprog -----
    %   .regulators
    %   .coefficients
    %   .regulatorinformation: array of reginfo
    %   ------ reginfo -----     
    %   .name
    %   .regulator
    %   .value
    %   ------ split -----    
    %   .type
    %   .splitdata_regulator
    %   .splitdata_point
    %   .splitdata_type
    %   .left/right: split or regprog
    %
    %
    
    para.objtype = 'simple';
    para.Set = {};
    para.Samples = {};
    para.regulators = {};
    para.coefficients = [];
    para.regulatorinformation = [];
    para.name = '';
    para.regulator = '';
    para.value = '';
    para.split = '';
    para.left = '';
    para.right = '';
    para.type = 'OneRegulator';
    para.splitdata_regulator = '';
    para.splitdata_point = 1.0;
    para.splitdata_type = 'less';    
    para.regprog = '';
    
    
    para = assignpara(para, varargin{:});
    
    para.objtype = objtype;
    obj.objtype = lower(para.objtype);
    
    objfield.simple = {'name','Set','Samples'};
    objfield.reg = {'name','Set','Samples', 'regprog'};
    objfield.modnet = {'name','Set','Samples','split'};
    objfield.splitreg = {'name','Set','Samples','split'};
        
    objfield.regprog = {'regulators','coefficients','regulatorinformation'};
    objfield.split = {'type', 'regprog', 'left', 'right', ...
        'splitdata_regulator','splitdata_point','splitdata_type'};   
    
    objfield.reginfo = {'name','regulator','value'};
    
    if ~ismember(obj.objtype, fieldnames(objfield))
        error('unknown objtype %s', objtype);
    end
    
    for i = 1:length(objfield.(obj.objtype))
        obj.(objfield.(obj.objtype){i}) = para.(objfield.(obj.objtype){i});
    end
            
    
    
    