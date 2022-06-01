%workflow of EvolveX simulations as in submission (18.4.2021)
%% model revision
%let's introduce the missing parts of Ehrlich pathway in the S. cerevisiae
%consensus model v. 7.6
model = add_Ehrlich_yeast76(yeast76);
%let's set growth as an objective
model.c(:)=0;
model.c(2985)=1;
%let's revise the flux bounds of FA activation reactions
[model] = revise_yeast76(model);
%let's introduce a glucose minimal medium (ammonium as nitrogen source)
model.lb(2665) = -10;
model.lb(2612) = -1000;

%mimic anaerobiosis
model.ub(335) = 0;
model.ub(780) = 0;
model.lb(780) = 0;
model.ub(203) = 0;
model.ub(526) = 0; %isocitrate lyase

%let's close threonine aldolase assumed to have minor role during growth 
%in minimal medium, close to equilibrium -7.2 +- 4.4 [kJ/mol] 
%(http://equilibrator.weizmann.ac.il)
model.ub(798) = 0;
%let's close FMN reductases with relevance in apoptosis
%if these are on, anaerobic conditions generate in simulations succinate and no glycerol
%this does not match with the real biological state
model.ub(338) = 0;
model.ub(339) = 0;
%let's close glycerol utilization pathway as glucose is used as c-source
model.ub(381) = 0;
model.ub(382) = 0;
model.lb(381) = 0;
model.lb(382) = 0;
%let's close higher alcohol acetate ester -esterases
model.ub(521) = 0;
model.ub(522) = 0;
%let's close phenylacetaldehyde exchange assuming redox status preferring
%further conversion to phenyl ethanol synthesis
model.ub(2902) = 0;
%let's constrain isocitrate dehydrogenase in NADPH generating direction
model.lb(524) = 0;
%glycine synthesis from serine when glucose is the carbon source
model.lb(391) = 0;
model.lb(392) = 0;
%C1 metabolism not the main source of NADPH
model.ub(589) = 0; %5,10-methylenetetrahydrofolate + NADP+ = 5,10-methenyltetrahydrofolate + NADPH + H+, deltaG ~10
model.ub(590) = 0;

%% EvolveX 
%let's prepare for the EvolveX runs
%let's prepare growth-constrained models
model_evolveX = model;
model_evolveX.ub(2549)=1000;
model_evolveX.ub(2549)=1000;
model_evolveX.ub(2985) = 10;
model_evolveX.lb(2985) = 10;
%let's close glucose and ammonium uptakes for the evolveX runs
model_evolveX.lb(2665) = 0;
model_evolveX.lb(2612) = 0;

%let's create growth constrained wt model for the evolvex runs 
model_evolveX_wt = model_evolveX;
model_evolveX_wt.lb(2665) = -1000;
model_evolveX_wt.lb(2612) = -1000;
model_evolveX_wt.c(:)=0;
%let's assign glucose and ammonium uptakes as objectives
model_evolveX_wt.c(2665)=1;
model_evolveX_wt.c(2612)=1;

%mimic anaerobiosis off, let's comapre against aerobic glc + nh4
model_evolveX_wt.ub(335) = 1000; %COX
model_evolveX_wt.ub(780) = 1000; %Succinate dehydrogenase
model_evolveX_wt.lb(780) = -1000; %Succinate dehydrogenase
model_evolveX_wt.ub(203) = 1000; %ATP synthase
model_evolveX_wt.ub(526) = 1000; %isocitrate lyase

%let's run FVA for the growth-constrained model minmax
[model_evolveX_wt_wtminmax(:,1), model_evolveX_wt_wtminmax(:,2), minStat, maxStat, time] = run_FVA(model_evolveX_wt, 'max', 100, model_evolveX_wt.rxns);

%let's allow for glycerol utilization as possible evolution niche
%component
model_evolveX.ub(382) = 1000;
model_evolveX.lb(382) = 0;
%let's allow the use of glutamate synthase
model_evolveX.ub(367) = 1000;
%mimic anaerobiosis
model_evolveX.ub(335) = 1000;
model_evolveX.ub(780) = 1000;
model_evolveX.lb(780) = -1000;
model_evolveX.ub(203) = 1000;
model_evolveX.ub(526) = 1000; %isocitrate lyase
%serine synthesis from glycine possible
model_evolveX.lb(391) = -1000;
model_evolveX.lb(392) = -1000;
model_evolveX.lb(524) = -1000;
%C1 metabolism not the main source of NADPH
model_evolveX.ub(589) = 1000; %5,10-methylenetetrahydrofolate + NADP+ = 5,10-methenyltetrahydrofolate + NADPH + H+, deltaG ~10
model_evolveX.ub(590) = 1000;

%let's pre-process by identifying those niche that support growth
growth_supporting_niche = get_growth_supporting_niche(model_evolveX,CD_extended);

for rxn=1:856
    %let's run evolvex for all the growth supporting niche
    [strength_rxns{rxn},coverage_rxns{rxn},ncomp_rxns{rxn},value_rxns{rxn},comp_rxns{rxn},bc_rxns{rxn}] = enumerate_evolvex(model_evolveX,CD_extended,growth_supporting_niche,[],rxn,{'UP'},model_evolveX_wt_wtminmax,1);
end

%[compcombBCAA1245] = evolvex_to_xlsx(strengthBCAA1245,coverageBCAA1245,ncompBCAA1245,valueBCAA1245,compBCAA1245,bcBCAA1245,'evolvex_enumerated_BCAA1245_anaerobic_160421.xlsx');
%[compcombPHE] = evolvex_to_xlsx(strengthPHE,coveragePHE,ncompPHE,valuePHE,compPHE,bcPHE,'evolvex_enumerated_PHE_anaerobic_160421.xlsx');
