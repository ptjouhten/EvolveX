function [modelrev] = revise_yeast76(model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
modelrev = model;
modelrev.lb(901:925) = 0; %indeces from yeast76gluc
modelrev.lb(328:332) = 0; %indeces from yeast76gluc

end

