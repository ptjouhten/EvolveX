function [gene_data] = get_change_for_genes(gene_set,flux_change,model)

gene_data = {};
for i = 1:length(gene_set)
    rxn_ind_list = find(sum(model.rxnGeneMat(:,find(strcmp(gene_set{i},model.genes))),2));
    for j = 1:length(rxn_ind_list)
        if j == 1
            gene_data{i,1} = flux_change{rxn_ind_list(j)};
        elseif ~strcmp(flux_change{rxn_ind_list(j)},gene_data{i,1})
            if strcmp('NO_CHANGE',gene_data{i,1})
               gene_data{i,1} = flux_change(rxn_ind_list(j));
            elseif strcmp('DOWN',flux_change{rxn_ind_list(j)})
               gene_data{i,1} = 'BOTH'; 
            elseif strcmp('UP',flux_change{rxn_ind_list(j)})
               gene_data{i,1} = 'BOTH'; 
            end
        end
    end
end
end

