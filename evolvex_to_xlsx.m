function [compcomb] = evolvex_to_xlsx(strength,coverage,ncomp,value,comp,bc,filename)
%post-processing of evolveX results into an .xlsx file

[compcomb] = parseComp(comp);
xlswrite(filename,compcomb,'Sheet1','A1')
xlswrite(filename,strength','Sheet1','B1')
xlswrite(filename,coverage','Sheet1','C1')
xlswrite(filename,value','Sheet1','E1')
xlswrite(filename,ncomp','Sheet1','D1')
xlswrite(filename,bc','Sheet1','F1')

end

