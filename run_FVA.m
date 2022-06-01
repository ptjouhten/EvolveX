function [minFlux, maxFlux, minStat, maxStat, time] = run_FVA(model, osense, objLevel, targetRxns)
%this function runs flux variability analysis for a set of reactions

if isempty(osense)
    osense = 'max';
end
if isempty(objLevel)
    objLevel = 100;
end
if isempty(targetRxns)
    targetRxns = model.rxns;
end

delta = 10^6;
minFlux = [];
maxFlux = [];

%run FBA
[solution] = FBApj(model,osense);
objValue = (objLevel/100) * sum(solution.x(logical(model.c)));
disp(objValue)
disp(solution.f)

%run min/max for the targetRxns
A = [full(model.S);model.c'];
if strcmp(osense,'max')
    if solution.f > 0
        %max and pos, relaxing should decrease the value (smaller positive)
        objValue = (1 - (100-objLevel)/100) * sum(solution.x(logical(model.c)));
        b = [full(model.b);floor(delta * objValue')/delta];
    else
        %max and neg, relaxing should decrease the value (bigger negative)
        objValue = (1 + (100-objLevel)/100) * sum(solution.x(logical(model.c)));
        b = [full(model.b);floor(delta * objValue')/delta];
    end
else
    if solution.f > 0
        %min and pos, relaxing should increase the value (bigger positive)
        objValue = (1 + (100-objLevel)/100) * sum(solution.x(logical(model.c)));
        b = [full(model.b);ceil(delta * objValue')/delta];
    else
        %min and neg, relaxing should increase the value (smaller negative)
        objValue = (1 - (100-objLevel)/100) * sum(solution.x(logical(model.c)));
        b = [full(model.b);ceil(delta * objValue')/delta];
    end
end



lb = full(model.lb);
ub = full(model.ub);
csense(1:(length(b)-1),1) = 'S';
%objective constrained to equal or greater to objValue in case of osense
%max
if strcmp(osense,'max')
    csense(length(b),1) = 'L';
else
    csense(length(b),1) = 'U';
end

%params.msglev = 3; % level of verbosity
%params.msglev = 0;
%params.tolbnd = 1e-9; %tolerance
%params.toldj = 1e-6; %tolerance

Aeq = A(csense == 'S',:);
beq = b(csense == 'S',:);
A(csense == 'L',:) = A(csense == 'L',:) *(-1);
b(csense == 'L',:) = b(csense == 'L',:) *(-1);
Aineq = A(csense ~= 'S',:);
bineq = b(csense ~= 'S',:);
    
options = cplexoptimset('cplex');
options.advance = 0;
options.read.scale = -1;
options.Display = 'off';
options.simplex.tolerances.optimality = 10^-9;

tic;
l = length(targetRxns);
for j= 1:l
    disp(j)
    c = zeros(length(model.rxns),1);
    ind = logical(strcmp(targetRxns{j},model.rxns));
    c(ind) = 1;
    
    %max
    f = -1*c;
    [x,fval,exitflag,output] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,[],options);
    
    maxFlux(j) = x(logical(c));
    maxStat(j) = exitflag;
    
    %min
    f = 1*c;
    [x,fval,exitflag,output] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,[],options);
    
    minFlux(j) = x(logical(c));
    minStat(j) = exitflag;

end
    maxFlux = maxFlux';
    minFlux = minFlux';
    maxStat = maxStat';
    minStat = minStat';
    time = toc;
end
