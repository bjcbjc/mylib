function [res] = rtpcr_foldchange(target, targetbase, ref, refbase, teff, reff)
    %calculate fold change for RT-PCR result
    %
    %target: {label Ct std} of target gene, #sample x 3,
    %targetbase: {label Ct std} of target gene of base strain, #sample (or 1) x 3
    %ref: {label Ct std} of reference gene, #sample x 3
    %refbase: {label Ct std} of reference gene of base strain, #sample (or 1) x 3
    %teff: efficiency of primers for target gene
    %reff: efficiency of primers for reference gene
    %
    %
    %
    %Example
    %target: 2 (swap and RM) x 3, Ct for PHO84
    %targetbase: 1 (BY) x 3, Ct for PHO84
    %ref: 2 (swap and RM) x 3, Ct for ERV25
    %refbase: 1 (BY) x 3, Ct for ERV25
    %
   
    
    nsample = size(target,1);
    if nsample ~=size(ref,1)
        error('number of samples is not the same in target and ref');
    end
    
    nbase = size(targetbase,1);
    if nbase ~= size(refbase,1)
        error('number of base strain is not the same in targetbase and refbase');
    end
    
    Et = log2(teff);
    Er = log2(reff);
    
    if nbase ~= nsample
        basemtx_t = repmat(cell2mat(targetbase(:,2:3)), nsample/nbase, 1);
        basemtx_r = repmat(cell2mat(refbase(:,2:3)), nsample/nbase, 1);
    else
        basemtx_t = cell2mat(targetbase(:,2:3));
        basemtx_r = cell2mat(refbase(:,2:3));
    end
    
    target_mtx = cell2mat(target(:,2:3));
    ref_mtx = cell2mat(ref(:,2:3));
    
    delta_t = -(target_mtx(:,1) - basemtx_t(:,1));
    delta_r = -(ref_mtx(:,1) - basemtx_r(:,1));

    %assume Ct of target gene in target strain has no correlation with ref strain!!!!!
    var_deltat = target_mtx(:,2).^2 + basemtx_t(:,2).^2; 
    var_deltar = ref_mtx(:,2).^2 + basemtx_r(:,2).^2;
    
    res.foldchange = delta_t*Et - delta_r*Er;
    res.std = sqrt(Et^2 .* var_deltat + Er^2 .* var_deltar);
    res.target = target(:,1);
    res.ref = ref(:,1);
    res.tbase = targetbase(:,1);
    res.rbase = refbase(:,1);
end