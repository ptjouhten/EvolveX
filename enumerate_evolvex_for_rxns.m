function [strengthDOWN,strengthUP,coverageDOWN,coverageUP] = enumerate_evolvex_for_rxns(model,components,inhibitors,inhibitor_rxn_matrix,reference_minmax,aerobic)
%This functions calls evolvex for each reaction in the model for UP and 
%DOWN selection in the given evolution niche

%for all the reactions on the model
tic
k = 0;
for j=1:length(model.rxns)
    display(j)
    if sum(model.rxnGeneMat(j,:)) ~= 0
        k = k+1;
        [strengthDOWN(j,1),coverageDOWN(j,1)] = evolvex(components,inhibitors,model,j,{'DOWN'},inhibitor_rxn_matrix,reference_minmax,aerobic);
        [strengthUP(j,1),coverageUP(j,1)] = evolvex(components,inhibitors,model,j,{'UP'},inhibitor_rxn_matrix,reference_minmax,aerobic);
    end
end
toc
end

