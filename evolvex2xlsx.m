function [compcomb] = evolvex2xlsx(strength,coverage,ncomp,value,comp,bc)
%post-processing of evolveX results into an .xlsx file

[compcomb] = parseComp(comp);
xlswrite('evolveXnew.xlsx',compcomb,'Sheet1','A1')
xlswrite('evolveXnew.xlsx',strength','Sheet1','B1')
xlswrite('evolveXnew.xlsx',coverage','Sheet1','C1')
xlswrite('evolveXnew.xlsx',value','Sheet1','E1')
xlswrite('evolveXnew.xlsx',ncomp','Sheet1','D1')
xlswrite('evolveXnew.xlsx',bc','Sheet1','F1')

end

