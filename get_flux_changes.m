function [change] = get_flux_changes(switchx_fullsol,wtminmax,eps)

%minThres = wtminmax(binRxns(i),1) - delta*abs(wtminmax(binRxns(i),1))-eps;
%maxThres = wtminmax(binRxns(i),2) + delta*abs(wtminmax(binRxns(i),2))+eps;
change = {};
for i = 1:size(wtminmax,1)
%compare fullsol flux values to wtminmax
if abs(switchx_fullsol(i)) - min(abs(wtminmax(i,:))) < -eps
   display(i)
   change{i,1} = 'DOWN';
elseif abs(switchx_fullsol(i)) - max(abs(wtminmax(i,:))) > eps
   display(i)
   change{i,1} = 'UP';
else
   change{i,1} = 'NO_CHANGE';
end

end

end

