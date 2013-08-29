function [xtraining,xtest,ytraining,ytest] = training_test_datasets (fold,input,x,y)

test_size = input (fold,size(input,2));
training_size = size(input,2) - 1 - test_size;

xtest = x(:,input(fold,1:test_size));
xtraining = x (:,input(fold,test_size + 1 : test_size + training_size));

ytest = y(:,input(fold,1:test_size));
ytraining = y (:,input(fold,test_size + 1 : test_size + training_size));


