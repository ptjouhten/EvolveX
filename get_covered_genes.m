function [genes1_covered_by_genes2] = get_covered_genes(gene_set1,gene_set2)

genes1_covered_by_genes2 = zeros(length(gene_set1),1);

for i =1:length(gene_set1)
   covered = find(strcmp(gene_set1{i},gene_set2));
   if ~isempty(covered)
       genes1_covered_by_genes2(i) = 1;
   end
end

end

