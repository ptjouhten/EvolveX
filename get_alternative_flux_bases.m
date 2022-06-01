function [altoptsols, indeces] = get_alternative_flux_bases(solMatrix, notOnlyOptimal)

nSols = size(solMatrix,2);
nTargets = [];
indeces = [];
altoptsols = [];


%calculate the number of targets in each solution
for i = 1:nSols
   nTargets(i) = length(find(solMatrix(:,i)));
   indeces(i) = i;
end

if ~isempty(nTargets)
    if notOnlyOptimal == 0
        minvalue = min(nTargets);

        indeces = indeces(nTargets==minvalue);
        %optimal solutions
        optSolMatrix = solMatrix(:,nTargets==minvalue);
    else
        optSolMatrix = solMatrix;
    end
    %unique optimal solutions
    [altoptsols, ai, ci] = unique(optSolMatrix','rows');
    altoptsols = altoptsols';
    indeces = indeces(ai);

end

end

