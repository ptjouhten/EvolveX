function [binSolp, fullSolp, dirBinSolp, solutionp] = prune_flux_basis_solutions(wtminmax, binSol, fullSol, binRxns, delta, eps, solution)
%This function prunes the flux basis solutions to remove the non-optimal 
%solutions returned by the solver

dirBinSol = zeros(size(binSol));
binSolp = [];
fullSolp = [];
prune = ones(size(binSol,2),1);

binSol = round(binSol);

%for columns
for j = 1:size(binSol,2)
    if solution{j}.origStat==1
        %for binary reactions
        for i = 1:length(binRxns)
            if binSol(i,j) == 1
                %is it smaller than wt minimum
                minThres = wtminmax(binRxns(i),1) - delta*abs(wtminmax(binRxns(i),1))-eps;
                maxThres = wtminmax(binRxns(i),2) + delta*abs(wtminmax(binRxns(i),2))+eps;
                if maxThres < minThres
                    disp('threshold error') 
                end
                if fullSol(binRxns(i),j) <= minThres
                    dirBinSol(i,j) = -1;
                %is it bigger than wt maximum
                elseif fullSol(binRxns(i),j) >= maxThres
                    dirBinSol(i,j) = 1;
                %other case, e.g. sign change, but no change in capacity
                else
                    disp(i)
                    prune(j) = 0;
                    dirBinSol(i,j) = 1000;
                end
            end
        end
    else
        prune(j) = 0;
    end
end
binSolp = binSol(:,logical(prune));
fullSolp = fullSol(:,logical(prune));
dirBinSolp = dirBinSol(:,logical(prune));
solutionp = solution(logical(prune));

end

