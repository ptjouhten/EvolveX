function [evolveX_targets, evolveX_target_dirs, gene_annotations] = parse_flux_basis_evolvex_input(model,rxns,dirs,signs)
%This function gets the gene annotations of the input rxns (indeces)
%from the model, and parses the rxns, dirs, and signs into desired trait 
%flux basis (UP and with gene annotations only) input for evolvex 

genes_matrix = model.rxnGeneMat(rxns,:);

gene_annotations = {};
e = 1;
j = 1;
for i = 1:size(genes_matrix,1)   
    genes = find(genes_matrix(i,:));
    for k = 1:length(genes)
        gene_annotations{j,1} = model.genes{genes(k)};
        gene_annotations{j,2} = dirs{i};
        gene_annotations{j,3} = signs(i);
        gene_annotations{j,4} = model.rxns{rxns(i)};
        gene_annotations{j,5} = model.rxnNames{rxns(i)};
        gene_annotations{j,6} = rxns(i);
        %let's collect UP change fluxes with gene annotations
        if strcmp('UP',dirs{i}) && signs(i)==1
           evolveX_targets_all(e)=rxns(i);
           evolveX_target_dirs_all{e}='UP';
           e = e+1;
        elseif strcmp('DOWN',dirs{i}) && signs(i)==-1
           evolveX_targets_all(e)=rxns(i);
           evolveX_target_dirs_all{e}='UP';
           e = e+1;
        end
        j = j+1;
    end
end
[evolveX_targets,ia,ic] = unique(evolveX_targets_all);
evolveX_target_dirs = evolveX_target_dirs_all(ia);

end

