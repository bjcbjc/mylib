function [pval, Z, O, E, V] = mylogrank(data1, data2)

t = unique([data1(data1(:,2)==0,1); data2(data2(:,2)==0,1)]);
O = 0;
E = 0;
V = 0;
% valid1 = true(size(data1,1),1);
% valid2 = true(size(data2,1),1);
for i = 1:length(t)
    idx1 = data1(:,1) == t(i);
    idx2 = data2(:,1) == t(i);
    observedDeath = [sum(data1(idx1,2)==0), sum(data2(idx2,2)==0)];    
    n = [sum(data1(:,1) >= t(i)), sum(data2(:,1) >= t(i))];    
    observedAlive = n - observedDeath;
    expected = n(1)*sum(observedDeath)/sum(n);
    v = n(1) * n(2) * sum(observedDeath) * sum(observedAlive);
    v = v / (sum(n)^2) / (sum(n)-1); 
    O = O + observedDeath(1);
    E = E + expected(1);
    V = V + v;
%     fprintf('%f, %f, %f\n',observedDeath(1), expected(1), v);
end
Z = (O - E) / sqrt(V);
pval = 2*normcdf( abs(Z), 0, 1, 'upper');
