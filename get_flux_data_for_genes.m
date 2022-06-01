function [gene_data] = get_flux_data_for_genes(gene_set,flux_data,model,bigger)

gene_data = zeros(length(gene_set),1);

for i = 1:length(gene_set)
    rxn_ind_list = find(sum(model.rxnGeneMat(:,find(strcmp(gene_set{i},model.genes))),2));
    for j = 1:length(rxn_ind_list)
        if j == 1
            gene_data(i) = flux_data(rxn_ind_list(j));
        elseif bigger == 1 && flux_data(rxn_ind_list(j)) > gene_data(i) %positive flux
            gene_data(i) = flux_data(rxn_ind_list(j));
        elseif bigger == 0 && flux_data(rxn_ind_list(j)) < gene_data(i) %negative flux
            gene_data(i) = flux_data(rxn_ind_list(j));
        end
    end
end


end

