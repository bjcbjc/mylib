function [out] = cross_validation (fold,num_samples)

perm = randperm (num_samples);

a = floor(num_samples / fold);
b = num_samples - a * fold;

indx = 0;
out = [];
for i = 1 : fold
	test_size = a + 1 - step(i-1-b);
	out = [out ; perm(indx+1 : indx+test_size) perm(1:indx) perm(indx+test_size+1:num_samples) test_size];
	indx = indx + test_size;
end

