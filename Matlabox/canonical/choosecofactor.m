
%pick cofactor for each marker
%%
fns = {'c54_TPO.mat','c54_IL3.mat','c54_IL7.mat','c54_DMSO.mat',...
    'c54_BCR.mat','c54_Dasatinib_Unstim.mat','c54_GCSF.mat','c54_U0126_Unstim.mat','c54_IFNad.mat'};
%datapath = 'C:\Users\Dylan\Desktop\mr_folder\';
datapath = '/Users/bjc/Desktop/Courses/Spring 2011/Data mining/project/data/';

cf = 1:8;
msethres = -0.1;

%%
chosencf = NaN(31, length(fns));
for fi = 1:length(fns)    
    mse = NaN(31, length(cf));
    for ci = 1:length(cf)
        data = preparedata(datapath, fns{fi}, cf(ci));
        mse(:,ci) = qqdeviation(data.mtx{1}(:,1:31));
        clear data
    end
    msecopy = mse;
    for i = 1:31        
        d = diff(mse(i,:));
        if all(d < 0) 
            j = find(d > msethres, 1, 'first');
            mse(i, j+2:end) = NaN;
        end
        [~, mi] = nanmin(mse(i,:));
        chosencf(i, fi) = cf(mi);        
    end
end

save('chosencf.mat', 'chosencf')
