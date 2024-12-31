function [binRxns] = getPotentialTargets(model, wtminmax, redVector)
% This function is a helper function for SwitchX
%
%INPUT:  model      engineered model containing the production pathway
%        wtminmax   flux ranges in wt phenotype (size: rxns x 2)
%        redVector  binary vector defining rxns being possible targets
%OUTPUT: binRxns    potential target rxns

% 23/08/2018 improved comments pjouhten

%[mets, rxns] = size(model.S);
[wtrxns] = size(wtminmax,1);
if isempty(redVector)
   redVector = ones(wtrxns,1); 
end
%check the endogenous rxns that were constrained by max growth, binary variables
%only for them, binrxns contains the indeces of reactions considered
%this also directly limits the binary variables to the endogenous reactions
%since heterologous reactions are added in the end of the reaction list
s = 1;
for r = 1:wtrxns
    if redVector(r) == 1
    if model.rev(r) == 1
        if ~(wtminmax(r,1) == -1000 && wtminmax(r,2) == 1000)
            binRxns(s) = r;
            s = s+1;
        end
    else
        if ~(wtminmax(r,1) == 0 && wtminmax(r,2) == 1000)
            binRxns(s) = r;
            s = s+1;
        end
    end
    end
end
%wtrxns = length(binRxns);
binRxns = binRxns';
end

