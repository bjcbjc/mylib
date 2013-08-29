function [] = sample_generation (fold,num_samples,filename)

fid = fopen(filename,'w');

perm = randperm (num_samples);

a = floor(num_samples / fold);
b = num_samples - a * fold;

indx = 0;
for i = 1 : fold

   test = a + 1 - step (i-1-b);
   out = perm(indx+1 : indx+test) - 1;
   indx = indx + test;

   fprintf(fid,'%d',out(1));
   for j=2:size(out,2)
      fprintf(fid,'	%d',out(j));
   end
   fprintf(fid,'\n');
end
fclose (fid);

clear out;

