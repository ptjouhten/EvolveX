function [strength, evolvex_minmax, sol, tarsol] = evolvex_fva(components,inhibitors,model,flux_basis,flux_basis_dirs,inhibitor_rxn_matrix,aerobic)
%This function scores the suitability of an evolution niche consisting of 
%components and inhibitors for adaptively evolving a desired 
%trait defined through its flux_basis and flux_basis_directions
%model should be a version of yeast consensus metabolic model v. 7.6 due to
%some hard coded reaction indeces, inhibitor_rxn_matrix maps inhibitors to 
%reactions, and reference_minmax contains the growth coupled flux ranges in
%reference conditions

%output:
%  strength     worst case response to selection of flux basis 
%  coverage     number of fluxes in flux basis having stronger response to
%               selection in evolution niche than in reference conditions
%  ncomp        number of components in the evolution niche
%  value        suitability score combined of the above
%  bc           
%  sol
%  tarsol

ncomp = length(unique(components));
bc = zeros(length(flux_basis),1);
mtol = 1e-9;
tol = 1e-9;

%response to selection coverage threshold
thres = 0.10;

model.lb = model.lb.*10;
model.ub = model.ub.*10;
model.lb(2985)=10;
model.ub(2985)=10;

%mapping the chemical environment into the bounds of the model
if aerobic == 1
    model.ub(335) = 10000;
    model.ub(780) = 10000;
    model.ub(203) = 10000;
else
    %no respiration, but oxygen for biosynthetic needs
    model.ub(335) = 0;
    model.ub(780) = 0;
    model.ub(203) = 0;
end
ntars = length(flux_basis);

%compounds in the chemical environment, allow uptake
model.lb(components) = -10000;
model.ub(components) = -1; %changed from 0 06/02/21
%create the separate secretion reactions
[rinit,cinit] = size(model.S);
%model.S(:,cinit+1:cinit+length(comp)) = zeros(rinit,length(comp));
model.S(:,cinit+1:cinit+length(components)) = model.S(:,components);
model.ub(cinit+1:cinit+length(components)) = 10000 * ones(length(components),1);
model.lb(cinit+1:cinit+length(components)) = zeros(length(components),1);
model.c(cinit+1:cinit+length(components)) = zeros(length(components),1);

%inhibitors in the chemical environment
if ~isempty(inhibitors)
%find out whether they are annotated to genes or rxns
[a, ninh] = size(inhibitor_rxn_matrix);
if a == length(model.genes)
%find the orfs affected, then identify rxns affected
[inhGenes, c] = find(inhibitor_rxn_matrix(:,inhibitors));
%find the rxns associated with the inhibited gene products using the
%rxnGeneMat
rxnInd = find(any(model.rxnGeneMat(:,inhGenes),2));
if (~isempty(rxnInd))
    x = true(size(model.genes));
    x(inhGenes) = false;
    inhRxn = false(length(rxnInd),1);
    %use the rules to evaluate whether any reaction state is affected by the inhibition?
    for j = 1:length(rxnInd)
        %this evaluation needs the x vector above
        if (~eval(model.rules{rxnInd(j)}))
            inhRxn(j) = true;
        end
    end
    %assume the inhibited rxns are fully inhibited
    if (any(inhRxn))
        model.ub(rxnInd(inhRxn)) = 0;
        model.lb(rxnInd(inhRxn)) = 0;
    end
end
elseif a == length(model.rxns)
    [inhRxns,c] = find(inhibitor_rxn_matrix(:,inhibitors));
    model.ub(inhRxns) = 0;
    model.lb(inhRxns) = 0;
else
    disp('Incorrect inhibitor annotation matrix! Inhibitors not considered')
end
end

model.ub(logical(model.c)) = 10;
model.lb(logical(model.c)) = 10;

growthc = model.c;

model.c(:) = 0;
model.c(components) = 1;

sol = FBApj(model,'max');

if isempty(sol.x)
    %penalize for non-growing solutions
    disp(sol.f)
    strength = -10^8;
    %continue to strength and coverage evaluation only if the chemical
    %environment supports growth
else
   %fluxes as variables
    [nmets, nrxns] = size(model.S);
    %formulate the objectve function to optimise the performance of the
    %target fluxes
    uptars = flux_basis(strcmp('UP',flux_basis_dirs)');
    nuptars = length(uptars);
    %check which target fluxes are reversible, they need special treatment below
    revuptars = flux_basis(strcmp('UP',flux_basis_dirs) & (model.lb(flux_basis) < 0));
    nrevuptars = length(revuptars);
    revdtars = flux_basis(strcmp('DOWN',flux_basis_dirs) & (model.lb(flux_basis) < 0));
    nrevdtars = length(revdtars);
    %reserve space
    A = zeros(2*nrevuptars+4*nrevdtars+2,nrxns+nrevuptars+2*nrevdtars);
    b = zeros(2*nrevuptars+4*nrevdtars+2,1);
    f = zeros(nrxns+nrevuptars+2*nrevdtars,1);
    
    Aeq = zeros(nmets,nrxns+nrevuptars+2*nrevdtars);
    
    %metabolite mass balance constraints
    Aeq(1:nmets,1:nrxns) = full(model.S);
    beq = zeros(nmets,1);
    %make a binary vector for noting the nutritional conditions
    condlog = zeros(size(model.S,2),1);
    condlog(components) = 1;
    
    vartype(1:length(f),1) = 'C';
    
    %lower and upper bounds for cplex
    lbevo = [model.lb; zeros(nrevuptars+2*nrevdtars,1)];
    ubevo = [model.ub; 10000*ones(nrevuptars+2*nrevdtars,1)];
    
    lbevo(1:nrxns,1) = model.lb;
    ubevo(1:nrxns,1) = model.ub;
    %make a vector for target indeces for considering the absolute value
    %variables to be introduced for reversible fluxes in flux basis
    targetsabs = zeros(length(flux_basis),1);
    
    k = 1;
    p = 1;
    for i = 1:ntars
       if strcmp(flux_basis_dirs{i},'UP')
            %the aim is to minimize the absolute values of the
           	%UP-regulation target fluxes given optimal conversion of nutrients to biomass (=worst case scenario)
            %if the flux can be negative, optimisation has to be
            %done on absolute values
            if ismember(flux_basis(i),revuptars)
                targetsabs(i) = nrxns + k; 
                %add a new variable into the objective function
                f(nrxns + k,1) = 1;
                ubevo(nrxns + k,1) = max(abs(model.ub(flux_basis(i))),abs(model.lb(flux_basis(i)))); %200416 pjouhten
                A((2*k-1),flux_basis(i)) = 1;
                A((2*k-1),nrxns+k) = -1;
                A((2*k),flux_basis(i)) = -1;
                A((2*k),nrxns+k) = -1;
                
                k = k+1;
            else
                %if the target flux is constrained to positive direction 
                targetsabs(i) = flux_basis(i);
                f(flux_basis(i),1) = 1;
            end
       else
           %if down-regulation target
           %the aim is to maximize the absolute values of the DOWN-regulation fluxes given optimal conversion of nutrients to biomass (=worst case scenario)
            if ismember(flux_basis(i), revdtars)
                targetsabs(i) = nrxns + nrevuptars + (2*p-1);
                %add the new variables into the objective function
                f(nrxns + nrevuptars + (2*p-1)) = -1;
                ubevo(nrxns + nrevuptars + (2*p-1)) = max(abs(model.ub(flux_basis(i))),abs(model.lb(flux_basis(i))));
                %add new constraints for the new variables
                A(2*nrevuptars + (4*p-3),flux_basis(i)) = 1;
                A(2*nrevuptars + (4*p-3),nrxns + nrevuptars + (2*p-1)) = -1; %v-u <= 0
                
                A(2*nrevuptars + (4*p-2),flux_basis(i)) = -1;
                A(2*nrevuptars + (4*p-2),nrxns + nrevuptars + (2*p-1)) = -1; %-v-u<=0
                
                A(2*nrevuptars + (4*p-1),flux_basis(i)) = -1;
                A(2*nrevuptars + (4*p-1),nrxns + nrevuptars + (2*p-1)) = 1;
                A(2*nrevuptars + (4*p-1),nrxns + nrevuptars + (2*p)) = -20000; %-v-M*B+u<=0
                
                A(2*nrevuptars + (4*p),flux_basis(i)) = 1;
                A(2*nrevuptars + (4*p),nrxns + nrevuptars + (2*p-1)) = 1;
                A(2*nrevuptars + (4*p),nrxns + nrevuptars + (2*p)) = 20000; %v+M*B+u<=20000
                b(2*nrevuptars + (4*p)) = 20000;
                
                vartype(nrxns+nrevuptars+(2*p),1) = 'B';
                ubevo(nrxns+nrevuptars+(2*p),1) = 1;
                p = p+1;
            else
                targetsabs(i) = flux_basis(i);
                f(flux_basis(i),1) = -1;
            end
       end
    end
    
    A((end-1),1:nrxns) = -1 .* growthc';
    b((end-1),1) = -10;
    A(end,1:nrxns) = -1 .* condlog';
    b(end,1) = -floor((model.c'*(sol.x))/tol)*tol;
    
    options = cplexoptimset;
    options.Display = 'off'; 
    options.Diagnostics = 'on';
    options.advance = 0;
    options.read.scale = -1;
    options.simplex.tolerances.optimality = 10^-9;
   
    if nrevdtars == 0
        %no binary variables if no reversible DOWN-regulation fluxes
        [x, fval, exitflag, output] = cplexlp(f, A, b, Aeq, beq, lbevo, ubevo, [], options);
    else
        [x,fval,exitflag,output] = cplexmilp(f,A,b,Aeq,beq,[],[],[],lbevo,ubevo,vartype',[],options);
    end
    
    if exitflag == 1
        tarsol.f = f'*x;
        tarsol.x = x;
        tarsol.stat = exitflag;
        strength = f'*tarsol.x;
    elseif exitflag == 5 && output.cplexstatus == 102
        tarsol.f = f'*x;
        tarsol.x = x;
        tarsol.stat = exitflag;
        strength = f'*tarsol.x;
    else
        tarsol.f = [];
        tarsol.x = [];
        disp(exitflag)
        disp(b(end,1))
        disp(sol.f)
        strength = -10^6;
    end
    
    for j = 1: length(model.rxns)
        f(:) = 0;
        %min
        f(j) = 1;
        [x, fval, exitflag, output] = cplexlp(f, A, b, Aeq, beq, lbevo, ubevo, [], options);
        evolvex_minmax(j,1) = x'*f;
        %max
        f(j) = -1;
        [x, fval, exitflag, output] = cplexlp(f, A, b, Aeq, beq, lbevo, ubevo, [], options);
        evolvex_minmax(j,2) = x'*-f;
    end
    %continue to coverage evaluation only if strength optimization was
    %successful
    %if strength > -10^6
    %    oldf = f;
    
        %create new empty f
    %    f(:)=0;
        %preallocation of space
    %    A(:,nrxns+nrevuptars+2*nrevdtars+1:nrxns+nrevuptars+2*nrevdtars+ntars) = zeros(size(A,1),ntars);
    %    A(2*nrevuptars+4*nrevdtars+3:2*nrevuptars+4*nrevdtars+2+ntars+1,:) = zeros(ntars+1,size(A,2));
    %    Aeq(:,end+1:end+ntars) = zeros(size(Aeq,1),ntars);
    %    b(2*nrevuptars+4*nrevdtars+3:2*nrevuptars+4*nrevdtars+2+ntars+1) = zeros(ntars+1,1);
        %minimization if worst coverage
    %    f(nrxns+nrevuptars+2*nrevdtars+1:nrxns+nrevuptars+2*nrevdtars+ntars,1) = ones(ntars,1);
    %    
    %    vartype(nrxns+nrevuptars+2*nrevdtars+1:nrxns+nrevuptars+2*nrevdtars+ntars,1) = 'B';
    %    ubevo(nrxns+nrevuptars+2*nrevdtars+1:nrxns+nrevuptars+2*nrevdtars+ntars) = ones(ntars,1);
    %    lbevo(nrxns+nrevuptars+2*nrevdtars+1:nrxns+nrevuptars+2*nrevdtars+ntars) = zeros(ntars,1);
   % 
   % 
   %     %optimization for coverage
   %     for i = 1:length(targetsabs)
   %     %if it is an UP flux
   %         if strcmp(flux_basis_dirs{i},'UP')
                %add a new constraint
                %worst (=min) coverage, b can only be zero if the
                %flux is below the threshold of (1+thres)*wt flux max
   %             A(2*nrevuptars+4*nrevdtars+2+i,targetsabs(i)) = 1;
   %             A(2*nrevuptars+4*nrevdtars+2+i,nrxns+nrevuptars+2*nrevdtars+i) = ((1+thres)*max(abs(reference_minmax(flux_basis(i),1)),abs(reference_minmax(flux_basis(i),2)))-ubevo(targetsabs(i)));
   %             b(2*nrevuptars+4*nrevdtars+2+i) = (1+thres)*max(abs(reference_minmax(flux_basis(i),1)),abs(reference_minmax(flux_basis(i),2)));
                
            %if down-regulation flux
   %         elseif strcmp(flux_basis_dirs{i},'DOWN')
                %add a new constraint               
                %worst (=min) coverage, b can only be zero if the
                %flux is above the threshold of (1-thres)*wt flux min
   %             A(2*nrevuptars+4*nrevdtars+2+i,targetsabs(i)) = -1;
   %             A(2*nrevuptars+4*nrevdtars+2+i,nrxns+nrevuptars+2*nrevdtars+i) = (-(1-thres)*min(abs(reference_minmax(flux_basis(i),1)),abs(reference_minmax(flux_basis(i),2)))+(lbevo(targetsabs(i))));
   %             b(2*nrevuptars+4*nrevdtars+2+i) = -(1-thres)*min(abs(reference_minmax(flux_basis(i),1)),abs(reference_minmax(flux_basis(i),2)));
   %         end
   %     end
    
        %under the worst case response to selection
   %     A(end,1:nrxns+nrevuptars+2*nrevdtars) = oldf';
   %     b(end,1) = ceil((oldf'*(tarsol.x))/tol)*tol;
    
   %     [xc,fvalc,exitflagc,outputc] = cplexmilp(f,A,b,Aeq,beq,[],[],[],lbevo,ubevo,vartype',[tarsol.x;zeros(ntars,1)],options);
    
   %     if exitflagc == 1
   %         covsol.f = f'*xc;
   %         covsol.x = xc;
   %         bc = xc(end-(ntars-1):end);
   %         covsol.stat = exitflagc;
   %         %for min coverage
   %         coverage = f'*covsol.x;
   %     elseif exitflagc == 5 && outputc.cplexstatus == 102
   %         covsol.f = f'*xc;
   %         covsol.x = xc;
   %         covsol.stat = exitflagc;
   %         %for min coverage
   %         coverage = f'*covsol.x;
   %     else
   %         disp(exitflagc)
   %         coverage = 0;
   %     end
   % 
   %     %formulation of the score
   %     %for worst-case selection pressure, worst-case coverage
   %     value = 1000 * (ntars-coverage)/ntars + (1000*nuptars - strength)/ntars + ncomp;
   % else
   %     coverage = 0;
   %     value = 10^6;
   % end
end
end
