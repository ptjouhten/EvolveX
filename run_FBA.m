function [solution] = FBApj(model,csense)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Aeq = full(model.S);
beq = zeros(length(model.mets),1);

f = model.c;

if isempty(csense)
    f = -1*f;
elseif strcmp(csense,'max')
    f = -1*f;
end

lb = full(model.lb);
ub = full(model.ub);

options = cplexoptimset;
options.Display = 'on';
options.simplex.tolerances.optimality = 10^-9;

[x, fval, exitflag, output] = cplexlp(f, [], [], Aeq, beq, lb, ub, [], options);
if exitflag == 1
solution.f = model.c'*x;
solution.x = x;
solution.stat = exitflag;
else
    solution.f = [];
    solution.x = [];
    solution.stat = exitflag;
end

end

