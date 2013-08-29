function [training_index, test_index] = training_test_indices (fold,input)

test_size = input (fold,size(input,2));
training_size = size(input,2) - 1 - test_size;

test_index = input (fold,1:test_size);
training_index = input (fold,test_size+1:test_size+training_size);


