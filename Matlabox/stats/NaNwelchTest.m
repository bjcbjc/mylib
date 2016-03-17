function [pValue, tValue] = NaNwelchTest(matrix1, matrix2)
n1 = sum( ~isnan(matrix1), 2);
n2= sum( ~isnan(matrix2), 2);

var1= nanvar(matrix1, 0, 2);

var2= nanvar(matrix2, 0, 2);


tValue = (nanmean(matrix1, 2)-nanmean(matrix2, 2))./(sqrt(var1./n1+var2./n2));

u=var2./var1;
temp = (1./n1+u./n2);
v= (temp.*temp)./(1./ (n1.^2.*(n1-1))+   u.*u./(n2.^2.*(n2-1))  );

pValue = tcdf(-abs(tValue),v)*2;
