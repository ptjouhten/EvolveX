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
%% trait's flux basis determination
%let's create a separate model version for that
model_flux_basis_wt = model;

%let's run flux variability analysis to determine wt flux ranges
[model_flux_basis_wt_wtminmax(:,1), model_flux_basis_wt_wtminmax(:,2), minStat, maxStat, time] = run_FVA(model_flux_basis_wt, 'max', 100, model_flux_basis_wt.rxns);
model_flux_basis_wt_wtminmax_unrounded = model_flux_basis_wt_wtminmax;
model_flux_basis_wt_wtminmax = round(model_flux_basis_wt_wtminmax.*10^9)/10^9;

%let's create separate models for the two desired traits 
%model for L-Phenylalanine derived aromas
model_flux_basis_PHE = model;
model_flux_basis_PHE.c(:)=0;
model_flux_basis_PHE.c(2901)=1;
model_flux_basis_PHE.c(2557)=1;
%let's block the L-serine cycling as a non-growing fermenting state 
%will be simulated, technical change as minimum sets of fluxes that have to
%change from the optimally growing state are determined
model_flux_basis_PHE.ub(688) = 0;

%ester ratios constraints
model_flux_basis_PHE.b(end+1) = 0;
model_flux_basis_PHE.mets(end+1) = {'flux_ratio_constraint1'};
model_flux_basis_PHE.S(end+1,2557) = 2; %2-Phenyl ethanol exchange
model_flux_basis_PHE.S(end,2901) = -1; %Phenylethyl acetate exchange

%extension starts here
%two models for branch chain amino acids' derived aromas to gather the flux
%bases of aromas of different stoichiometric costs 
model_flux_basis_BCAA1 = model;
%let's block the L-serine cycling as a non-growing fermenting state 
%will be simulated, technical change as minimum sets of fluxes that have to
%change from the optimally growing state are determined
model_flux_basis_BCAA1.ub(688) = 0;

model_flux_basis_BCAA1.mets(end+1) = {'flux_ratio_constraint2'};
model_flux_basis_BCAA1.b(end+1) = 0;
model_flux_basis_BCAA1.S(end+1,2549) = 2; %2-Methylbutanol exchange
model_flux_basis_BCAA1.S(end,2550) = -1; %2-Methylbutyl acetate exchange

model_flux_basis_BCAA2 = model;
%let's block the L-serine cycling as a non-growing fermenting state 
%will be simulated, technical change as minimum sets of fluxes that have to
%change from the optimally growing state are determined
model_flux_basis_BCAA2.ub(688) = 0;

model_flux_basis_BCAA2.b(end+1) = 0;
model_flux_basis_BCAA2.mets(end+1) = {'flux_ratio_constraint3'};
model_flux_basis_BCAA2.S(end+1,2792) = 2; %Isoamyl alcohol exchange
model_flux_basis_BCAA2.S(end,2789) = -1; %Isoamyl alcohol acetate exchange

model_flux_basis_BCAA1.c(:)=0;
model_flux_basis_BCAA2.c(:)=0;
model_flux_basis_BCAA1.c(2550)=1;
model_flux_basis_BCAA1.c(2549)=1;
model_flux_basis_BCAA2.c(2792)=1;
model_flux_basis_BCAA2.c(2789)=1;

%let's determine the minimum sets of fluxes that have to change for aroma 
%production (FBA solutions)
%the simulation needs three parameters, alpha, delta, and epsilon
%trait flux bases for evolvex will be determined from these simulations
%with alpha 1, delta 1, and epsilon 0.00001
[targetsPHE_110, dirsPHE_110, signsPHE_110, dirSolpPHE_110, binRxnsPHE_110, fullSolOptPHE_110, solutionpPHE_110] = determine_trait_flux_basis(model_flux_basis_PHE, 'max', model_flux_basis_wt_wtminmax, [], 1, 1, 0.00001, 5);
model_flux_basis_BCAA1.ub(2549)=1000;
model_flux_basis_BCAA2.ub(2549)=1000;
[targetsBCAA_1101, dirsBCAA_1101, signsBCAA_1101, dirSolpBCAA_1101, binRxnsBCAA_1101, fullSolOptBCAA_1101, solutionpBCAA_1101] = determine_trait_flux_basis(model_flux_basis_BCAA1, 'max', model_flux_basis_wt_wtminmax, [], 1, 1, 0.00001, 5);
[targetsBCAA_1102, dirsBCAA_1102, signsBCAA_1102, dirSolpBCAA_1102, binRxnsBCAA_1102, fullSolOptBCAA_1102, solutionpBCAA_1102] = determine_trait_flux_basis(model_flux_basis_BCAA2, 'max', model_flux_basis_wt_wtminmax, [], 1, 1, 0.00001, 5);

%let's parse the flux bases with absolute flux increase required here into 
%input for evolvex, this takes into account the flux signs
[evolveX_targets_PHE, evolveX_dirs_PHE, gene_annotations_PHE_110] = parse_flux_basis_evolvex_input(model_flux_basis_PHE,targetsPHE_110{1,1}',dirsPHE_110{1,1}',signsPHE_110{1,1}');
[evolveX_targets_BCAA1, evolveX_dirs_BCAA1, gene_annotations_BCAA_1101] = parse_flux_basis_evolvex_input(model_flux_basis_BCAA1,targetsBCAA_1101{1,1}',dirsBCAA_1101{1,1}',signsBCAA_1101{1,1}');
[evolveX_targets_BCAA2, evolveX_dirs_BCAA2, gene_annotations_BCAA_1102] = parse_flux_basis_evolvex_input(model_flux_basis_BCAA2,targetsBCAA_1102{1,1}',dirsBCAA_1102{1,1}',signsBCAA_1102{1,1}');

%let's combine the flux bases and get unique reactions
evolveX_targets_BCAA_1245 = [evolveX_targets_BCAA1'; evolveX_targets_BCAA2'];
evolveX_dirs_BCAA_1245 = [evolveX_dirs_BCAA1'; evolveX_dirs_BCAA2'];
[evolveX_targets_BCAA_1245u,evolveX_targets_BCAA_1245uind,evolveX_targets_BCAA_1245ucind]=unique(evolveX_targets_BCAA_1245);
evolveX_dirs_BCAA_1245u = evolveX_dirs_BCAA_1245(evolveX_targets_BCAA_1245uind);
evolveX_targets_PHE = evolveX_targets_PHE';
evolveX_dirs_PHE = evolveX_dirs_PHE';
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

%condMap & growth allowing combinations of 1, 2, and 3 
%condMap_ = condMapGminYeast76red(1:(end-1));
%C contains the reduced condMap with rxnNames

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

%let's set arbitrary ester to corresponding alcohol ratio constraints,
%without an excretion of the stoichiometrically cheaper, alcohol, would 
%always be preferred in the simulations
[mets, rxns] = size(model_evolveX.S);
model_evolveX.S(mets+1:mets+3,:) = zeros(3,rxns);
model_evolveX.b(mets+1:mets+3) = zeros(3,1);
model_evolveX.mets(mets+1) = {'flux_ratio_constraint1'};
model_evolveX.S(mets+1,2557) = 2; %2-Phenyl ethanol exchange
model_evolveX.S(mets+1,2901) = -1; %Phenylethyl acetate exchange
model_evolveX.mets(mets+2) = {'flux_ratio_constraint2'};
model_evolveX.S(mets+2,2549) = 2; %2-Methylbutanol exchange
model_evolveX.S(mets+2,2550) = -1; %2-Methylbutyl acetate exchange
model_evolveX.mets(mets+3) = {'flux_ratio_constraint3'};
model_evolveX.S(mets+3,2792) = 2; %Isoamyl alcohol exchange
model_evolveX.S(mets+3,2789) = -1; %Isoamyl alcohol acetate exchange

%let's pre-process by identifying those niche that support growth
growth_supporting_niche = get_growth_supporting_niche(model_evolveX,CD_limited);

%let's run evolvex for all the growth supporting niche
[strengthBCAA1245,coverageBCAA1245,ncompBCAA1245,valueBCAA1245,compBCAA1245,bcBCAA1245] = enumerate_evolvex(model_evolveX,CD_limited,growth_supporting_niche,[],evolveX_targets_BCAA_1245u,evolveX_dirs_BCAA_1245u,model_evolveX_wt_wtminmax,1);
[strengthPHE,coveragePHE,ncompPHE,valuePHE,compPHE,bcPHE] = enumerate_evolvex(model_evolveX,CD_limited,growth_supporting_niche,[],evolveX_targets_PHE,evolveX_dirs_PHE,model_evolveX_wt_wtminmax,1);

%[compcombBCAA1245] = evolvex_to_xlsx(strengthBCAA1245,coverageBCAA1245,ncompBCAA1245,valueBCAA1245,compBCAA1245,bcBCAA1245,'evolvex_enumerated_BCAA1245_anaerobic.xlsx');
%[compcombPHE] = evolvex_to_xlsx(strengthPHE,coveragePHE,ncompPHE,valuePHE,compPHE,bcPHE,'evolvex_enumerated_PHE_anaerobic.xlsx');
[compcombBCAA1245] = evolvex_to_xlsx(strengthBCAA1245,coverageBCAA1245,ncompBCAA1245,valueBCAA1245,compBCAA1245,bcBCAA1245,'evolvex_enumerated_BCAA1245_anaerobic_160421.xlsx');
[compcombPHE] = evolvex_to_xlsx(strengthPHE,coveragePHE,ncompPHE,valuePHE,compPHE,bcPHE,'evolvex_enumerated_PHE_anaerobic_160421.xlsx');
%% After evolution

%let's predict response to selection for all the reactions in the model
%in the selected evolution niche to compare against the RNA-seq and
%proteomics data on evolved clones in evolution niche
%Note: the worst case selection pressure constraint is different when 
%evolve is run for each reactions individually than in case of a set of 
%fluxes in a flux basis
[strengthDOWN_PHE_G,strengthUP_PHE_G,coverageDOWN_PHE_G,coverageUP_PHE_G] = enumerate_evolvex_for_rxns(model_evolveX,[2737,2828,2836],[],[],model_evolveX_wt_wtminmax,1);
[strengthDOWN_BCAA_E,strengthUP_BCAA_E,coverageDOWN_BCAA_E,coverageUP_BCAA_E] = enumerate_evolvex_for_rxns(model_evolveX,[2707,2806,2739],[],[],model_evolveX_wt_wtminmax,1);

%let's get the genes annotated to the reactions having stronger response 
%to selection in the evolution niche than in reference conditions
%let's first prune the reaction binary vectors
coverageDOWN_PHE_G_pruned = prune_binary_vector(coverageDOWN_PHE_G);
coverageUP_PHE_G_pruned = prune_binary_vector(coverageUP_PHE_G);
coverageDOWN_BCAA_E_pruned = prune_binary_vector(coverageDOWN_BCAA_E);
coverageUP_BCAA_E_pruned = prune_binary_vector(coverageUP_BCAA_E);
%let's then get the genes annotated to the reactions
coverageDOWN_genes_PHE_G = yeast76.genes(find(sum(yeast76.rxnGeneMat(find(coverageDOWN_PHE_G_pruned),:))));
coverageUP_genes_PHE_G = yeast76.genes(find(sum(yeast76.rxnGeneMat(find(coverageUP_PHE_G_pruned),:))));
coverageDOWN_genes_BCAA_E = yeast76.genes(find(sum(yeast76.rxnGeneMat(find(coverageDOWN_BCAA_E_pruned),:))));
coverageUP_genes_BCAA_E = yeast76.genes(find(sum(yeast76.rxnGeneMat(find(coverageUP_BCAA_E_pruned),:))));


%E31E3283
[metabolic] = get_covered_genes(proteomics_E31E3283_evolution_medium(:,1),yeast76.genes);
proteomics_metabolic_E31E3283_evolution_medium = proteomics_E31E3283_evolution_medium(logical(metabolic),:);
[proteomics_genes_E31E3283_covered_evolution_medium(:,1)] = get_covered_genes(proteomics_metabolic_E31E3283_evolution_medium(:,1),coverageUP_genes_BCAA_E);
[proteomics_genes_E31E3283_covered_evolution_medium(:,2)] = get_covered_genes(proteomics_metabolic_E31E3283_evolution_medium(:,1),coverageDOWN_genes_BCAA_E);

[metabolic] = get_covered_genes(rnaseq_E31E3283_evolution_medium(:,1),yeast76.genes);
rnaseq_metabolic_E31E3283_evolution_medium = rnaseq_E31E3283_evolution_medium(logical(metabolic),:);
[rnaseq_genes_E31E3283_covered_evolution_medium(:,1)] = get_covered_genes(rnaseq_metabolic_E31E3283_evolution_medium(:,1),coverageUP_genes_BCAA_E);
[rnaseq_genes_E31E3283_covered_evolution_medium(:,2)] = get_covered_genes(rnaseq_metabolic_E31E3283_evolution_medium(:,1),coverageDOWN_genes_BCAA_E);

%E31G321
[metabolic] = get_covered_genes(proteomics_E31G321_evolution_medium(:,1),yeast76.genes);
proteomics_metabolic_E31G321_evolution_medium = proteomics_E31G321_evolution_medium(logical(metabolic),:);
[proteomics_genes_E31G321_covered_evolution_medium(:,1)] = get_covered_genes(proteomics_metabolic_E31G321_evolution_medium(:,1),coverageUP_genes_PHE_G);
[proteomics_genes_E31G321_covered_evolution_medium(:,2)] = get_covered_genes(proteomics_metabolic_E31G321_evolution_medium(:,1),coverageDOWN_genes_PHE_G);

[metabolic] = get_covered_genes(rnaseq_E31G321_evolution_medium(:,1),yeast76.genes);
rnaseq_metabolic_E31G321_evolution_medium = rnaseq_E31G321_evolution_medium(logical(metabolic),:);
[rnaseq_genes_E31G321_covered_evolution_medium(:,1)] = get_covered_genes(rnaseq_metabolic_E31G321_evolution_medium(:,1),coverageUP_genes_PHE_G);
[rnaseq_genes_E31G321_covered_evolution_medium(:,2)] = get_covered_genes(rnaseq_metabolic_E31G321_evolution_medium(:,1),coverageDOWN_genes_PHE_G);

%E32E3246
[metabolic] = get_covered_genes(proteomics_E32E3246_evolution_medium(:,1),yeast76.genes);
proteomics_metabolic_E32E3246_evolution_medium = proteomics_E32E3246_evolution_medium(logical(metabolic),:);
[proteomics_genes_E32E3246_covered_evolution_medium(:,1)] = get_covered_genes(proteomics_metabolic_E32E3246_evolution_medium(:,1),coverageUP_genes_BCAA_E);
[proteomics_genes_E32E3246_covered_evolution_medium(:,2)] = get_covered_genes(proteomics_metabolic_E32E3246_evolution_medium(:,1),coverageDOWN_genes_BCAA_E);

[metabolic] = get_covered_genes(rnaseq_E32E3246_evolution_medium(:,1),yeast76.genes);
rnaseq_metabolic_E32E3246_evolution_medium = rnaseq_E32E3246_evolution_medium(logical(metabolic),:);
[rnaseq_genes_E32E3246_covered_evolution_medium(:,1)] = get_covered_genes(rnaseq_metabolic_E32E3246_evolution_medium(:,1),coverageUP_genes_BCAA_E);
[rnaseq_genes_E32E3246_covered_evolution_medium(:,2)] = get_covered_genes(rnaseq_metabolic_E32E3246_evolution_medium(:,1),coverageDOWN_genes_BCAA_E);

%E33G3213
[metabolic] = get_covered_genes(proteomics_E33G3213_evolution_medium(:,1),yeast76.genes);
proteomics_metabolic_E33G3213_evolution_medium = proteomics_E33G3213_evolution_medium(logical(metabolic),:);
[proteomics_genes_E33G3213_covered_evolution_medium(:,1)] = get_covered_genes(proteomics_metabolic_E33G3213_evolution_medium(:,1),coverageUP_genes_PHE_G);
[proteomics_genes_E33G3213_covered_evolution_medium(:,2)] = get_covered_genes(proteomics_metabolic_E33G3213_evolution_medium(:,1),coverageDOWN_genes_PHE_G);

[metabolic] = get_covered_genes(rnaseq_E33G3213_evolution_medium(:,1),yeast76.genes);
rnaseq_metabolic_E33G3213_evolution_medium = rnaseq_E33G3213_evolution_medium(logical(metabolic),:);
[rnaseq_genes_E33G3213_covered_evolution_medium(:,1)] = get_covered_genes(rnaseq_metabolic_E33G3213_evolution_medium(:,1),coverageUP_genes_PHE_G);
[rnaseq_genes_E33G3213_covered_evolution_medium(:,2)] = get_covered_genes(rnaseq_metabolic_E33G3213_evolution_medium(:,1),coverageDOWN_genes_PHE_G);


%let's vary alpha and delta and keep epsilon as 0.00001, in addition to 
%the alpha 1, delta 1 predictions used as input for EvolveX for comparison
%against proteomics and RNA-seq data on evolved clones in wine must
[targetsPHE_195, dirsPHE_195, signsPHE_195, dirSolpPHE_195, binRxnsPHE_195, fullSolOptPHE_195, solutionpPHE_195] = determine_trait_flux_basis(model_flux_basis_PHE, 'max', model_flux_basis_wt_wtminmax, [], 1, 0.95, 0.00001, 5);
[targetsBCAA_1951, dirsBCAA_1951, signsBCAA_1951, dirSolpBCAA_1951, binRxnsBCAA_1951, fullSolOptBCAA_1951, solutionpBCAA_1951] = determine_trait_flux_basis(model_flux_basis_BCAA1, 'max', model_flux_basis_wt_wtminmax, [], 1, 0.95, 0.00001, 5);
[targetsBCAA_1952, dirsBCAA_1952, signsBCAA_1952, dirSolpBCAA_1952, binRxnsBCAA_1952, fullSolOptBCAA_1952, solutionpBCAA_1952] = determine_trait_flux_basis(model_flux_basis_BCAA2, 'max', model_flux_basis_wt_wtminmax, [], 1, 0.95, 0.00001, 5);

[targetsPHE_120, dirsPHE_120, signsPHE_120, dirSolpPHE_120, binRxnsPHE_120, fullSolOptPHE_120, solutionpPHE_120] = determine_trait_flux_basis(model_flux_basis_PHE, 'max', model_flux_basis_wt_wtminmax, [], 1, 0.2, 0.00001, 5);
[targetsBCAA_1201, dirsBCAA_1201, signsBCAA_1201, dirSolpBCAA_1201, binRxnsBCAA_1201, fullSolOptBCAA_1201, solutionpBCAA_1201] = determine_trait_flux_basis(model_flux_basis_BCAA1, 'max', model_flux_basis_wt_wtminmax, [], 1, 0.2, 0.00001, 8);
[targetsBCAA_1202, dirsBCAA_1202, signsBCAA_1202, dirSolpBCAA_1202, binRxnsBCAA_1202, fullSolOptBCAA_1202, solutionpBCAA_1202] = determine_trait_flux_basis(model_flux_basis_BCAA2, 'max', model_flux_basis_wt_wtminmax, [], 1, 0.2, 0.00001, 5);

[targetsPHE_0210, dirsPHE_0210, signsPHE_0210, dirSolpPHE_0210, binRxnsPHE_0210, fullSolOptPHE_0210, solutionpPHE_0210] = determine_trait_flux_basis(model_flux_basis_PHE, 'max', model_flux_basis_wt_wtminmax, [], 0.2, 1, 0.00001, 5);
[targetsBCAA_02101, dirsBCAA_02101, signsBCAA_02101, dirSolpBCAA_02101, binRxnsBCAA_02101, fullSolOptBCAA_02101, solutionpBCAA_02101] = determine_trait_flux_basis(model_flux_basis_BCAA1, 'max', model_flux_basis_wt_wtminmax, [], 0.2, 1, 0.00001, 5);
[targetsBCAA_02102, dirsBCAA_02102, signsBCAA_02102, dirSolpBCAA_02102, binRxnsBCAA_02102, fullSolOptBCAA_02102, solutionpBCAA_02102] = determine_trait_flux_basis(model_flux_basis_BCAA2, 'max', model_flux_basis_wt_wtminmax, [], 0.2, 1, 0.00001, 5);

[targetsPHE_0110, dirsPHE_0110, signsPHE_0110, dirSolpPHE_0110, binRxnsPHE_0110, fullSolOptPHE_0110, solutionpPHE_0110] = determine_trait_flux_basis(model_flux_basis_PHE, 'max', model_flux_basis_wt_wtminmax, [], 0.1, 1, 0.00001, 5);
[targetsBCAA_01101, dirsBCAA_01101, signsBCAA_01101, dirSolpBCAA_01101, binRxnsBCAA_01101, fullSolOptBCAA_01101, solutionpBCAA_01101] = determine_trait_flux_basis(model_flux_basis_BCAA1, 'max', model_flux_basis_wt_wtminmax, [], 0.1, 1, 0.00001, 5);
[targetsBCAA_01102, dirsBCAA_01102, signsBCAA_01102, dirSolpBCAA_01102, binRxnsBCAA_01102, fullSolOptBCAA_01102, solutionpBCAA_01102] = determine_trait_flux_basis(model_flux_basis_BCAA2, 'max', model_flux_basis_wt_wtminmax, [], 0.1, 1, 0.00001, 5);

%let's get the gene annotations of the desired traits' flux bases to be 
%compared against proteomics and rna-seq differentially abundant sets on 
%evolved clones in wine must
[up_genesPHE_110, down_genesPHE_110] = get_up_down_genes(model_flux_basis_PHE,targetsPHE_110{1,1},dirsPHE_110{1,1});
[up_genesPHE_0210, down_genesPHE_0210] = get_up_down_genes(model_flux_basis_PHE,targetsPHE_0210{1,1},dirsPHE_0210{1,1});
[up_genesPHE_0110, down_genesPHE_0110] = get_up_down_genes(model_flux_basis_PHE,targetsPHE_0110{1,1},dirsPHE_0110{1,1});
[up_genesPHE_120, down_genesPHE_120] = get_up_down_genes(model_flux_basis_PHE,targetsPHE_120{1,1},dirsPHE_120{1,1});
[up_genesPHE_195, down_genesPHE_195] = get_up_down_genes(model_flux_basis_PHE,targetsPHE_195{1,1},dirsPHE_195{1,1});

[up_genesBCAA_1101, down_genesBCAA_1101] = get_up_down_genes(model_flux_basis_BCAA1,targetsBCAA_1101{1,1},dirsBCAA_1101{1,1});
[up_genesBCAA_02101, down_genesBCAA_02101] = get_up_down_genes(model_flux_basis_BCAA1,targetsBCAA_02101{1,1},dirsBCAA_02101{1,1});
[up_genesBCAA_01101, down_genesBCAA_01101] = get_up_down_genes(model_flux_basis_BCAA1,targetsBCAA_01101{1,1},dirsBCAA_01101{1,1});
[up_genesBCAA_1201, down_genesBCAA_1201] = get_up_down_genes(model_flux_basis_BCAA1,targetsBCAA_1201{1,1},dirsBCAA_1201{1,1});
[up_genesBCAA_1951, down_genesBCAA_1951] = get_up_down_genes(model_flux_basis_BCAA1,targetsBCAA_1951{1,1},dirsBCAA_1951{1,1});

[up_genesBCAA_1102, down_genesBCAA_1102] = get_up_down_genes(model_flux_basis_BCAA2,targetsBCAA_1102{1,1},dirsBCAA_1102{1,1});
[up_genesBCAA_02102, down_genesBCAA_02102] = get_up_down_genes(model_flux_basis_BCAA2,targetsBCAA_02102{1,1},dirsBCAA_02102{1,1});
[up_genesBCAA_01102, down_genesBCAA_01102] = get_up_down_genes(model_flux_basis_BCAA2,targetsBCAA_01102{1,1},dirsBCAA_01102{1,1});
[up_genesBCAA_1202, down_genesBCAA_1202] = get_up_down_genes(model_flux_basis_BCAA2,targetsBCAA_1202{1,1},dirsBCAA_1202{1,1});
[up_genesBCAA_1952, down_genesBCAA_1952] = get_up_down_genes(model_flux_basis_BCAA2,targetsBCAA_1952{1,1},dirsBCAA_1952{1,1});

%let's analyse the overlap of the predictions and the RNA-seq and
%proteomics data on evolved clones

%clone E2-2 (E31E3283)
%wine must
%let's get metabolic enzymes from the set of differentially abundant
%proteins
[metabolic] = get_covered_genes(proteomics_E31E3283_wine_must(:,1),yeast76.genes);
proteomics_metabolic_E31E3283_wine_must = proteomics_E31E3283_wine_must(logical(metabolic),:);

[proteomics_genes_E31E3283_covered1(:,1)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1101);
[proteomics_genes_E31E3283_covered1(:,2)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1101);
[proteomics_genes_E31E3283_covered1(:,3)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1201);
[proteomics_genes_E31E3283_covered1(:,4)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1201);
[proteomics_genes_E31E3283_covered1(:,5)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_02101);
[proteomics_genes_E31E3283_covered1(:,6)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_02101);
[proteomics_genes_E31E3283_covered1(:,7)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_01101);
[proteomics_genes_E31E3283_covered1(:,8)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_01101);
[proteomics_genes_E31E3283_covered1(:,9)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1951);
[proteomics_genes_E31E3283_covered1(:,10)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1951);

[proteomics_genes_E31E3283_covered2(:,1)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1102);
[proteomics_genes_E31E3283_covered2(:,2)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1102);
[proteomics_genes_E31E3283_covered2(:,3)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1202);
[proteomics_genes_E31E3283_covered2(:,4)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1202);
[proteomics_genes_E31E3283_covered2(:,5)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_02102);
[proteomics_genes_E31E3283_covered2(:,6)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_02102);
[proteomics_genes_E31E3283_covered2(:,7)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_01102);
[proteomics_genes_E31E3283_covered2(:,8)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_01102);
[proteomics_genes_E31E3283_covered2(:,9)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1952);
[proteomics_genes_E31E3283_covered2(:,10)] = get_covered_genes(proteomics_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1952);

%let's get metabolic genes from the set of differentially abundant
%mRNAs
[metabolic] = get_covered_genes(rnaseq_E31E3283_wine_must(:,1),yeast76.genes);
rnaseq_metabolic_E31E3283_wine_must = rnaseq_E31E3283_wine_must(logical(metabolic),:);

[rnaseq_genes_E31E3283_covered1(:,1)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1101);
[rnaseq_genes_E31E3283_covered1(:,2)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1101);
[rnaseq_genes_E31E3283_covered1(:,3)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1201);
[rnaseq_genes_E31E3283_covered1(:,4)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1201);
[rnaseq_genes_E31E3283_covered1(:,5)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_02101);
[rnaseq_genes_E31E3283_covered1(:,6)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_02101);
[rnaseq_genes_E31E3283_covered1(:,7)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_01101);
[rnaseq_genes_E31E3283_covered1(:,8)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_01101);
[rnaseq_genes_E31E3283_covered1(:,9)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1951);
[rnaseq_genes_E31E3283_covered1(:,10)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1951);

[rnaseq_genes_E31E3283_covered2(:,1)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1102);
[rnaseq_genes_E31E3283_covered2(:,2)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1102);
[rnaseq_genes_E31E3283_covered2(:,3)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1202);
[rnaseq_genes_E31E3283_covered2(:,4)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1202);
[rnaseq_genes_E31E3283_covered2(:,5)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_02102);
[rnaseq_genes_E31E3283_covered2(:,6)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_02102);
[rnaseq_genes_E31E3283_covered2(:,7)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_01102);
[rnaseq_genes_E31E3283_covered2(:,8)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_01102);
[rnaseq_genes_E31E3283_covered2(:,9)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),up_genesBCAA_1952);
[rnaseq_genes_E31E3283_covered2(:,10)] = get_covered_genes(rnaseq_metabolic_E31E3283_wine_must(:,1),down_genesBCAA_1952);

%clone G2-2 (E31G321)
%wine must
%let's get metabolic enzymes from the set of differentially abundant
%proteins
[metabolic] = get_covered_genes(proteomics_E31G321_wine_must(:,1),yeast76.genes);
proteomics_metabolic_E31G321_wine_must = proteomics_E31G321_wine_must(logical(metabolic),:);

[proteomics_genes_E31G321_covered(:,1)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),up_genesPHE_110);
[proteomics_genes_E31G321_covered(:,2)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),down_genesPHE_110);
[proteomics_genes_E31G321_covered(:,3)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),up_genesPHE_120);
[proteomics_genes_E31G321_covered(:,4)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),down_genesPHE_120);
[proteomics_genes_E31G321_covered(:,5)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),up_genesPHE_0210);
[proteomics_genes_E31G321_covered(:,6)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),down_genesPHE_0210);
[proteomics_genes_E31G321_covered(:,7)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),up_genesPHE_0110);
[proteomics_genes_E31G321_covered(:,8)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),down_genesPHE_0110);
[proteomics_genes_E31G321_covered(:,9)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),up_genesPHE_195);
[proteomics_genes_E31G321_covered(:,10)] = get_covered_genes(proteomics_metabolic_E31G321_wine_must(:,1),down_genesPHE_195);

%let's get metabolic genes from the set of differentially abundant
%mRNAs
[metabolic] = get_covered_genes(rnaseq_E31G321_wine_must(:,1),yeast76.genes);
rnaseq_metabolic_E31G321_wine_must = rnaseq_E31G321_wine_must(logical(metabolic),:);

[rnaseq_genes_E31G321_covered(:,1)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),up_genesPHE_110);
[rnaseq_genes_E31G321_covered(:,2)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),down_genesPHE_110);
[rnaseq_genes_E31G321_covered(:,3)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),up_genesPHE_120);
[rnaseq_genes_E31G321_covered(:,4)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),down_genesPHE_120);
[rnaseq_genes_E31G321_covered(:,5)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),up_genesPHE_0210);
[rnaseq_genes_E31G321_covered(:,6)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),down_genesPHE_0210);
[rnaseq_genes_E31G321_covered(:,7)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),up_genesPHE_0110);
[rnaseq_genes_E31G321_covered(:,8)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),down_genesPHE_0110);
[rnaseq_genes_E31G321_covered(:,9)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),up_genesPHE_195);
[rnaseq_genes_E31G321_covered(:,10)] = get_covered_genes(rnaseq_metabolic_E31G321_wine_must(:,1),down_genesPHE_195);

%clone E2-1 (E32E3246)
%wine must
%let's get metabolic enzymes from the set of differentially abundant
%proteins
[metabolic] = get_covered_genes(proteomics_E32E3246_wine_must(:,1),yeast76.genes);
proteomics_metabolic_E32E3246_wine_must = proteomics_E32E3246_wine_must(logical(metabolic),:);

[proteomics_genes_E32E3246_covered1(:,1)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1101);
[proteomics_genes_E32E3246_covered1(:,2)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1101);
[proteomics_genes_E32E3246_covered1(:,3)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1201);
[proteomics_genes_E32E3246_covered1(:,4)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1201);
[proteomics_genes_E32E3246_covered1(:,5)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_02101);
[proteomics_genes_E32E3246_covered1(:,6)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_02101);
[proteomics_genes_E32E3246_covered1(:,7)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_01101);
[proteomics_genes_E32E3246_covered1(:,8)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_01101);
[proteomics_genes_E32E3246_covered1(:,9)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1951);
[proteomics_genes_E32E3246_covered1(:,10)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1951);

[proteomics_genes_E32E3246_covered2(:,1)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1102);
[proteomics_genes_E32E3246_covered2(:,2)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1102);
[proteomics_genes_E32E3246_covered2(:,3)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1202);
[proteomics_genes_E32E3246_covered2(:,4)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1202);
[proteomics_genes_E32E3246_covered2(:,5)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_02102);
[proteomics_genes_E32E3246_covered2(:,6)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_02102);
[proteomics_genes_E32E3246_covered2(:,7)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_01102);
[proteomics_genes_E32E3246_covered2(:,8)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_01102);
[proteomics_genes_E32E3246_covered2(:,9)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1952);
[proteomics_genes_E32E3246_covered2(:,10)] = get_covered_genes(proteomics_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1952);

%let's get metabolic genes from the set of differentially abundant
%mRNAs
[metabolic] = get_covered_genes(rnaseq_E32E3246_wine_must(:,1),yeast76.genes);
rnaseq_metabolic_E32E3246_wine_must = rnaseq_E32E3246_wine_must(logical(metabolic),:);

[rnaseq_genes_E32E3246_covered1(:,1)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1101);
[rnaseq_genes_E32E3246_covered1(:,2)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1101);
[rnaseq_genes_E32E3246_covered1(:,3)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1201);
[rnaseq_genes_E32E3246_covered1(:,4)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1201);
[rnaseq_genes_E32E3246_covered1(:,5)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_02101);
[rnaseq_genes_E32E3246_covered1(:,6)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_02101);
[rnaseq_genes_E32E3246_covered1(:,7)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_01101);
[rnaseq_genes_E32E3246_covered1(:,8)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_01101);
[rnaseq_genes_E32E3246_covered1(:,9)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1951);
[rnaseq_genes_E32E3246_covered1(:,10)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1951);

[rnaseq_genes_E32E3246_covered2(:,1)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1102);
[rnaseq_genes_E32E3246_covered2(:,2)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1102);
[rnaseq_genes_E32E3246_covered2(:,3)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1202);
[rnaseq_genes_E32E3246_covered2(:,4)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1202);
[rnaseq_genes_E32E3246_covered2(:,5)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_02102);
[rnaseq_genes_E32E3246_covered2(:,6)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_02102);
[rnaseq_genes_E32E3246_covered2(:,7)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_01102);
[rnaseq_genes_E32E3246_covered2(:,8)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_01102);
[rnaseq_genes_E32E3246_covered2(:,9)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),up_genesBCAA_1952);
[rnaseq_genes_E32E3246_covered2(:,10)] = get_covered_genes(rnaseq_metabolic_E32E3246_wine_must(:,1),down_genesBCAA_1952);


%clone G2-1 (E33G3213)
%wine must
%let's get metabolic enzymes from the set of differentially abundant
%proteins
[metabolic] = get_covered_genes(proteomics_E33G3213_wine_must(:,1),yeast76.genes);
proteomics_metabolic_E33G3213_wine_must = proteomics_E33G3213_wine_must(logical(metabolic),:);

[proteomics_genes_E33G3213_covered(:,1)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),up_genesPHE_110);
[proteomics_genes_E33G3213_covered(:,2)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),down_genesPHE_110);
[proteomics_genes_E33G3213_covered(:,3)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),up_genesPHE_120);
[proteomics_genes_E33G3213_covered(:,4)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),down_genesPHE_120);
[proteomics_genes_E33G3213_covered(:,5)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),up_genesPHE_0210);
[proteomics_genes_E33G3213_covered(:,6)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),down_genesPHE_0210);
[proteomics_genes_E33G3213_covered(:,7)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),up_genesPHE_0110);
[proteomics_genes_E33G3213_covered(:,8)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),down_genesPHE_0110);
[proteomics_genes_E33G3213_covered(:,9)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),up_genesPHE_195);
[proteomics_genes_E33G3213_covered(:,10)] = get_covered_genes(proteomics_metabolic_E33G3213_wine_must(:,1),down_genesPHE_195);

%let's get metabolic genes from the set of differentially abundant
%mRNAs
[metabolic] = get_covered_genes(rnaseq_E33G3213_wine_must(:,1),yeast76.genes);
rnaseq_metabolic_E33G3213_wine_must = rnaseq_E33G3213_wine_must(logical(metabolic),:);

[rnaseq_genes_E33G3213_covered(:,1)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),up_genesPHE_110);
[rnaseq_genes_E33G3213_covered(:,2)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),down_genesPHE_110);
[rnaseq_genes_E33G3213_covered(:,3)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),up_genesPHE_120);
[rnaseq_genes_E33G3213_covered(:,4)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),down_genesPHE_120);
[rnaseq_genes_E33G3213_covered(:,5)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),up_genesPHE_0210);
[rnaseq_genes_E33G3213_covered(:,6)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),down_genesPHE_0210);
[rnaseq_genes_E33G3213_covered(:,7)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),up_genesPHE_0110);
[rnaseq_genes_E33G3213_covered(:,8)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),down_genesPHE_0110);
[rnaseq_genes_E33G3213_covered(:,9)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),up_genesPHE_195);
[rnaseq_genes_E33G3213_covered(:,10)] = get_covered_genes(rnaseq_metabolic_E33G3213_wine_must(:,1),down_genesPHE_195);


