function [up_genes, down_genes] = get_up_down_genes(model,targets,dirs)

up_genes = {};
down_genes = {};
dirs_abs = dirs;

up_targets = targets(strcmp('UP',dirs_abs));
down_targets = targets(strcmp('DOWN',dirs_abs));

up_genes = model.genes(logical(sum(model.rxnGeneMat(up_targets,:),1)));
down_genes = model.genes(logical(sum(model.rxnGeneMat(down_targets,:),1)));

end

