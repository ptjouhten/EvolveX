function [targets, dirs, signs, dirSolp, binRxns, fullSolOpt, solutionp] = determine_trait_flux_basis(model, osense, wtminmax, redVector, alpha, delta, eps, noSols)
%This function calls switchx to enumerate alternative minimum flux sets
%that need to change from one phenotype to another
%INPUT:  model      second phenotype model
%        wtminmax   flux ranges in the current phenotype (size: rxns x 2)
%        alpha      optimality in second phenotype -parameter
%        delta      relative flux change parameter
%        eps        additive flux change parameter
%        noSols     number of solutions to be enumerated
%OUTPUT: targets    rxn indeces of fluxes that need to change
%        dirs       required flux change directions (UP/DOWN)
%        signs      signs of fluxes that need to change
%        dirSolp    directions of flux changes in integers
%        binRxns    fluxes for which the change was evaluated
%        fullSolOpt solution vectors with binary and flux variables
%        solutionp  alternative optimal solutions

limit = 0;
fullSol = [];
binSol = [];
dirSol = [];
binRxns = [];
fullSol = [];
targets = {};
solution = {};

%let's round for the wtminmax
wtminmax = round((10^9)*wtminmax)/10^9;

%let's then enumerate alternative solutions
ultimate_limit = 3*noSols;
while limit < noSols
    %binSol is appended with new solutions
    [solution, binSol] = trait_flux_basis(model, osense, wtminmax, redVector, alpha, binSol, solution, delta, eps);
    cfullSol = size(fullSol,2); %number of previous solutions
    if length(solution) > cfullSol && ~isempty(solution{cfullSol+1}.full)
        fullSol(:,cfullSol+1) = solution{cfullSol+1}.full;
        limit = limit + 1;
    else
        solution = solution(1:cfullSol);
    end
    ultimate_limit = ultimate_limit -1;
    if ultimate_limit == 0
        limit = noSols;
    end
end

%get a logical defining the reactions for which the potential need-to-change-in-flux was
%evaluated, inbuilt also in trait_flux basis
[binRxns] = getPotentialTargets(model, wtminmax, redVector);

if ~isempty(binSol)
    %prune to remove the non-optimal returned by the solver
    [binSolp, fullSolp, dirThresp, solutionp] = prune_flux_basis_solutions(wtminmax, binSol, fullSol, binRxns, delta, eps, solution);
    %identify alternative solutions
    [altoptsols, indeces] = get_alternative_flux_bases(binSolp,1);
    %if solutions found
    if ~isempty(indeces)
        %find the required directions of flux change: -1 DOWN (below min threshold), 
        %1 UP (above max threshold), -1000 DOWN with sign change, 1000 UP with a
        %sign change, 5000 only a sign change
        [targets, dirs, signs, dirSolp] = get_flux_basis_directions(wtminmax, binSol(:,indeces), fullSol(:,indeces), binRxns, dirThresp(:,indeces), delta, eps);
        fullSolOpt = fullSol(:,indeces);
        solutionp = solutionp{indeces};
    else
        targets = [];
        dirs = [];
        signs = [];
        dirSolp = [];
        fullSolOpt = []; 
        solutionp = {};
    end
else
    targets = [];
    dirs = [];
    signs = [];
    dirSolp = [];
    fullSolOpt = []; 
    solutionp = {};
end
end

