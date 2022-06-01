function [solutionA, binSol] = trait_flux_basis(model,osense,wtminmax, redVector, alpha, binSol, solutionA, delta, eps)
% minPhenoSwitch algorithm implements an algorithm derived from the 
% regulatory on/off minimization (ROOM) algorithm by Shlomi et al.
% to identify the minimal set of reactions whose fluxes need to change to 
% to switch from one phenotype (wt) to another.
% Possible prior solution is made infeasible to find alternative solutions.
%
%INPUT:  model      the engineered model containing the production pathway
%        wtminmax   flux ranges in the wild type phenotype (size: rxns x 2)
%        alpha      parameter defining the required level optimal production
%        solver     use 'CPLEX'
%        delta      parameter for the threshold for a flux change, [0 1] 
%        eps        parameter for the threshold of a flux change, 0.001
%        noSols     number of solutions to be enumerated
%OUTPUT: solution   solution structure from the solver
%        binSol     binary solutions defining the identified targets

solution = [];
thr = 10^9;
solver = 'cplex';
%% solve optimal product yield
Aeq = full(model.S);
beq = zeros(length(model.mets),1);

%assumes positive coefficient in model
f = model.c;

if isempty(osense)
    f = -1*f;
elseif strcmp(osense,'max')
    f = -1*f;
end

lb = full(model.lb);
ub = full(model.ub);

options = cplexoptimset;
%options.Display = 'on';
%options.Display = 'off'; 
options.simplex.tolerances.optimality = 10^-9;

[x, fval, exitflag, output] = cplexlp(f, [], [], Aeq, beq, lb, ub, [], options);
if exitflag == 1
    objValue = model.c'*x;
else
    objValue = [];
end

if strcmp(osense,'max')
    if objValue > 0
        %max and pos, relaxing should decrease the value (smaller positive)
        vopt = (1 - (1-alpha)) * objValue;
        vopt = floor(thr * (vopt)')/thr;
    else
        %max and neg, relaxing should decrease the value (bigger negative)
        vopt = (1 + (1-alpha)) * objValue;
        vopt = floor(thr * (vopt)')/thr;
    end
else
    if objValue > 0
        %min and pos, relaxing should increase the value (bigger positive)
        vopt = (1 + (1-alpha)) * objValue;
        vopt = ceil(thr * (vopt)')/thr;
    else
        %min and neg, relaxing should increase the value (smaller negative)
        vopt = (1 - (1-alpha)) * objValue;
        vopt = ceil(thr * (vopt)')/thr;
    end
end

%% setting up the problem
[mets, rxns] = size(model.S);
[wtrxns] = size(wtminmax,1);
if isempty(redVector)
   redVector = ones(wtrxns,1); 
end

wtminmax = round((10^9)*wtminmax)/(10^9);

%create binary variables only for the fluxes that were not free under
%optimality, and were not listed to be left out in redVector
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
wtrxns = length(binRxns);

%calculate upper and lower bounds for a substantial change from wt flux
wiu = wtminmax(:,2) + delta*abs(wtminmax(:,2)) + eps;
wil = wtminmax(:,1) - delta*abs(wtminmax(:,1)) - eps;

%minimisation of the sum of binary variables
c = [zeros(1,rxns) ones(1,wtrxns)];

A = [model.S zeros(mets,wtrxns)];
%constraint for production
A(mets+1,:) = [model.c' zeros(1,wtrxns)];
%space for the upper bound constraints for a substantial change from wt flux
A(mets+2:(mets+1+wtrxns),:) = zeros(wtrxns,(rxns+wtrxns));
%space for the lower bound constraints for a substantial change from wt flux
A((mets+1+wtrxns+1):(mets+1+2*wtrxns),:) = zeros(wtrxns,(rxns+wtrxns));

b = zeros(length(model.mets),1);
b(mets+1) = vopt;
b((mets+2):(mets+1+2*wtrxns)) = [wiu(binRxns)' wil(binRxns)'];

%constraints for a substantial change in flux
for i = 1:wtrxns
    A((mets+1)+i,binRxns(i)) = 1;
    A((mets+1)+i,(rxns+i)) = wiu(binRxns(i))-model.ub(binRxns(i));
    A((mets+1+wtrxns)+i,binRxns(i)) = 1;
    A((mets+1+wtrxns)+i,(rxns+i)) = wil(binRxns(i))-model.lb(binRxns(i));
end

csense(1:mets,1) = 'S';
csense(mets+1,1) = 'L';
csense(mets+2:(mets+1+wtrxns),1) = 'U';
csense(mets+wtrxns+2:(mets+1+2*wtrxns),1) = 'L';

%% enumeration of alternative equally optimal solutions
% add a new constraint for each previous binary solution to not to allow
% for the already found solutions
% binSol has each previous binary solution in columns
if ~isempty(binSol)
    nrPrev = size(binSol,2);
    % for each previous solution add a constraint row
    for p = 1:nrPrev
        rA = size(A,1);
        A(rA+1,:) = [zeros(1,rxns) binSol(:,p)'];
        b = [b; sum(binSol(:,p))-1];       
        csense = [csense; 'U'];        
    end
else
    nrPrev = 0;
end


lb = [model.lb; zeros(wtrxns,1)];
ub = [model.ub; ones(wtrxns,1)];

vartype(1:rxns) = 'C';
vartype(rxns+1:(rxns+wtrxns)) = 'B';

osense = 1;
        
%% CPLEX
if strcmp(solver,'cplex')
        
    Aeq = A(csense == 'S',:);
    beq = b(csense == 'S',:);
    A(csense == 'L',:) = A(csense == 'L',:) *(-1);
    b(csense == 'L',:) = b(csense == 'L',:) *(-1);
    Aineq = A(csense ~= 'S',:);
    bineq = b(csense ~= 'S',:);
    f = (osense*c)';
    
   options = cplexoptimset('cplex');
   options.advance = 0;
   options.read.scale = -1;
   options.Display = 'off';
   %options.simplex.tolerances.optimality = 10^-6;
   options.mip.tolerances.integrality = 10^-9;
 
   [x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,vartype,[],options);
   disp(fval)

   solution.obj = fval;
   solution.stat = output;
   solution.origStat = exitflag;
   solution.full = x;
    
   if ~isempty(x)
       solution.flux = x(1:rxns,:);
       solution.targets = binRxns(logical(round(10^4.*x(vartype == 'B')/10^4)))';
       binSol(:,nrPrev+1) = x(vartype == 'B');
       solutionA{nrPrev+1} = solution;
   else
       solution.flux = [];
       solution.targets = [];
   end

end
end

