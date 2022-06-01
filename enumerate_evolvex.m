function [strength,coverage,ncomp,value,comp,bc] = enumerate_evolvex(model,condMap,growth,inh,targets,targetDirs,glcnh3minmax,aer)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nCond = length(condMap);

%all combinations reduced to the ones supporting growth without inhibition
%C3=nchoosek(condMap,3);
%C3 = C3(logical(growth{3}),:);
C2=nchoosek(condMap,2);
C2 = C2(logical(growth{2}),:);
%C1=condMap;
%C1 = C1(logical(growth{1}),:);

strength = [];
coverage = [];
comp = {};

%for i = 1:size(C3,1)
%    %disp(i)
%    [strength(i),coverage(i),ncomp(i),value(i),bc(:,i)] = evolvex(C3(i,:),inh,model,targets,targetDirs,[],glcnh3minmax,aer);
%    comp(i) = {model.rxnNames(C3(i,:))};
%end
%for i = (size(C3,1)+1):(size(C3,1)+size(C2,1))
for i = 1:size(C2,1)
    %disp(i)
    %[strength(i),coverage(i),ncomp(i),value(i),bc(:,i)] = evolvex(C2((i-size(C3,1)),:),inh,model,targets,targetDirs,[],glcnh3minmax,aer);
    [strength(i),coverage(i),ncomp(i),value(i),bc(:,i)] = evolvex(C2(i,:),inh,model,targets,targetDirs,[],glcnh3minmax,aer);
    %comp(i) = {model.rxnNames(C2((i-size(C3,1)),:))};
    comp(i) = {model.rxnNames(C2(i,:))};
end
%for i = (size(C3,1)+size(C2,1)+1):(size(C3,1)+size(C2,1)+size(C1,1))
%    %disp(i)
%    [strength(i),coverage(i),ncomp(i),value(i),bc(:,i)] = evolvex(C1((i-size(C3,1)-size(C2,1)),:),inh,model,targets,targetDirs,[],glcnh3minmax,aer);
%    comp(i) = {model.rxnNames(C1((i-size(C3,1)-size(C2,1)),:))};
%end
end

