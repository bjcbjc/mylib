function [pValue, tValue] = NaNwelchTest(matrix1, matrix2)
n1 = sum( (1-isnan(matrix1))');
n2= sum( (1-isnan(matrix2))');
var1= nanvar(matrix1');

var2= nanvar(matrix2');


tValue = (nanmean(matrix1')-nanmean(matrix2'))./(sqrt(var1./n1+var2./n2));

u=var2./var1;
temp = (1./n1+u./n2);
v= (temp.*temp)./(1./ (n1.^2.*(n1-1))+   u.*u./(n2.^2.*(n2-1))  );

pValue = tcdf(-abs(tValue),v)'*2;
