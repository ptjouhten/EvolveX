function [growth, uptakes] = get_growth_supporting_niche(model,condMap)

%all combinations of 1, 2, or 3 components in evolution niche
C{1}=condMap;
C{2}=nchoosek(condMap,2);
C{3}=nchoosek(condMap,3);

growth{1} = zeros(size(C{1},1),1);
growth{2} = zeros(size(C{2},1),1);
growth{3} = zeros(size(C{3},1),1);

uptakes{1} = cell(size(C{1},1),1);
uptakes{2} = cell(size(C{2},1),1);
uptakes{3} = cell(size(C{3},1),1);

%increase the bounds as below a high arbitrary growth constraint is used
lb = model.lb*10;
ub = model.ub*10;
growthc = logical(model.c);
%growth set to 10
model.ub(growthc) = 10;
model.lb(growthc) = 10;
%ub(growthc) = 10;
%lb(growthc) = 10;

k=0;
%%for i = 1:3
 %for i = 3
 for i = 2
    for j = 1:size(C{i},1)
        k = k+1;
        disp(k)
        model.lb = lb;
        model.ub = ub;

        %compounds in the chemical environment
        model.lb(C{i}(j,:)) = -10000;
        %model.ub(C{i}(j,:)) = 0;
        %26.2.2021 pj to align with evolvex algorithm
        model.ub(C{i}(j,:)) = -1;

        model.c(:) = 0;
        model.c(C{i}(j,:)) = 1;

        sol = run_FBA(model,'max');

        if ~isempty(sol.x)
            u = sol.x(C{i}(j,:));
            display(u) 
            growth{i}(j) = 1;
            uptakes{i}{j} = {u};
        end
    end
end
end

