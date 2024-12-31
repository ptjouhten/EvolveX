function [env] = parseComp(comp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i=1:size(comp,2)
    c = [];
    for j=1:length(comp{i})
        c = strcat(c,strcat(comp{i}(j),' ')); 
    end
    env(i,1) = c;
end

end

