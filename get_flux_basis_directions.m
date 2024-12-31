function [targets, dirs, signs, dirSolp] = get_flux_basis_directions(wtminmax, binSol, fullSol, binRxns, dirThresp, delta, eps)
% This function does post-processing to collect the needed flux change
% directions

dirSolp = zeros(size(dirThresp));

binSol = round(binSol);
fullSol = round((10^9)*fullSol)/10^9;
signChange = ones(size(dirThresp));
signs = cell(size(dirThresp,2));
%for alternative equally optimal solutions
for j = 1:size(binSol,2)
    k = 1;
    %for reactions with binary variables 
    for i = 1:length(binRxns)
        %if the flux has to change
        if binSol(i,j) == 1
            %get the thresholds
            minThres = wtminmax(binRxns(i),1) - delta*abs(wtminmax(binRxns(i),1))-eps;
            maxThres = wtminmax(binRxns(i),2) + delta*abs(wtminmax(binRxns(i),2))+eps;
            if maxThres < minThres
               disp('threshold error') 
            end
            %is there a sign change?
            if (dirThresp(i,j) == -1) && (sign(minThres) * sign(fullSol(binRxns(i),j)) < 0)
               signChange(i,j) = 1000;
               disp(binRxns(i))
               disp(fullSol(binRxns(i),j))
            elseif (dirThresp(i,j) == 1) && (sign(maxThres) * sign(fullSol(binRxns(i),j)) < 0)
               signChange(i,j) = 1000;
               disp(binRxns(i))
               disp(fullSol(binRxns(i),j))
            end
            %is it smaller than min threshold?
            if fullSol(binRxns(i),j) <= minThres
                dirSolp(i,j) = -1 * signChange(i,j);
                dirs{j}{k} = 'DOWN';
                targets{j}(k) = binRxns(i);
                signs{j}(k) = sign(fullSol(binRxns(i),j));
                k = k + 1;
            %is it bigger than max threshold?
            elseif fullSol(binRxns(i),j) >= maxThres
                dirSolp(i,j) = 1 * signChange(i,j);
                dirs{j}{k} = 'UP';
                targets{j}(k) = binRxns(i);
                signs{j}(k) = sign(fullSol(binRxns(i),j));
                k = k + 1;
            %other case, e.g. sign but not capacity change?
            else
                disp(i)
                dirSolp(i,j) = 5000;
            end
        end
    end
end

end

