////////////////////////////////////////////////////////////////////////////////
// SENSITIVITY ANALYSIS OF BLT FROM BREATHING AND ANXIETY CLINICAL STUDY ///////
////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------- //
// Author: Sophie Morris													  //
// Created: 2/12/2025														  //
// -------------------------------------------------------------------------- //
// Description: This script runs the sensitivity analysis of BLT behavioural  //
// data for Sophie M's Masters' project. The data was collected as part of    //
// the breathing and anxiety study conducted collaboratively with the         //
// Department of Psychology, School of Sport, Physical Education and Sport    //
// Science, and the School of Pharmacy at the University of Otago. Data from  //
// both clinical cohorts are used within this analysis. A timestamped         //
// analysis plan is available on (https://github.com/IMAGEotago). 			  //
// -------------------------------------------------------------------------- //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Import & Prepare Data ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Graph settings
graph set window fontface "Helvetica"
graph set ps fontface "Helvetica"
graph set svg fontface "Helvetica"
set graphics off

// Set working directory
cd /Users/sophiemorris/Desktop/study/analysis

// Combine BLT + questionnaire datasets with group labels
tempfile data_beh_all

import excel using "matlab/blt_beh_data_all.xlsx", sheet("Subjects") firstrow clear
gen group = "treatment"
replace group = "control" if substr(PPID, 1, 2) == "35"
save `data_beh_all'

import excel using "stata/data/results-clinQuest.xlsx", sheet("Sums_Data") firstrow clear
gen group = "treatment"
append using `data_beh_all'
save `data_beh_all', replace

import excel using "stata/data/results-clinControlQuest.xlsx", sheet("Sums_Data") firstrow clear
gen group = "control"
append using `data_beh_all'
save `data_beh_all', replace

// Generate session variable and clean PPIDs
gen session = "pre"
replace session = "post" if strpos(PPID, "-2")
replace PPID = substr(PPID, 1, strpos(PPID, "-2") - 1) if session == "post"

// Keep only variables of interest & collapse across duplicates
keep PPID session group STAI_T GAD7 muTheta muPhi AvgCertainty
collapse (mean) STAI_T GAD7 muTheta muPhi AvgCertainty, by(PPID group session)

// Keep only complete sessions
egen n_miss = rowmiss(STAI_T GAD7 muTheta muPhi AvgCertainty)
keep if n_miss == 0
drop n_miss

// Keep only participants with both pre and post
bysort PPID: egen n_sessions = count(session)
keep if n_sessions == 2
drop n_sessions

save "stata/data/data_beh_all.dta", replace

// Variable list
global myvars muTheta muPhi AvgCertainty GAD7 STAI_T



////////////////////////////////////////////////////////////////////////////////
// Summary Measures ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

use "stata/data/data_beh_all.dta", clear

// Descriptive summary
describe
sort group session
by group session: summarize

// For tabstat
gen str SessionGroup = session + "_" + group

// Summaries by session × group to Excel
tabstat $myvars, by(SessionGroup) stats(mean p50 sd n) ///
    columns(statistics) long save
	
return list

putexcel set "stata/results/sens/sensSummaryStat.xlsx", replace

putexcel set "stata/results/sens/sensSummaryStat.xlsx", sheet("post_control") modify
putexcel A1 = matrix(r(Stat1)), names

putexcel set "stata/results/sens/sensSummaryStat.xlsx", sheet("post_treatment") modify
putexcel A1 = matrix(r(Stat2)), names

putexcel set "stata/results/sens/sensSummaryStat.xlsx", sheet("pre_control") modify
putexcel A1 = matrix(r(Stat3)), names

putexcel set "stata/results/sens/sensSummaryStat.xlsx", sheet("pre_treatment") modify
putexcel A1 = matrix(r(Stat4)), names



////////////////////////////////////////////////////////////////////////////////
// Repeated Measures Observations //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Histograms by group and session
foreach var of global myvars {
    histogram `var', ///
		percent bin(20) ///
        by(session group, col(2) note("") ///
        title("Distribution of `var' by Group & Session (Sensitivity)")) ///
        ytitle("Percent") xtitle("Value") ///
        name(hist_`var', replace)
    graph export "stata/plots/sens/sens_hist_`var'.pdf", replace
}

// Clean histograms for combined figure
foreach var of global myvars {
    histogram `var', ///
		percent bin(10) ///
        by(session group, col(2) note("") ///
		title("`var'")) ///
        ytitle("") ///
		xtitle("") ///
        name(sens_hist_clean_`var', replace)
}

local graphlist ""

foreach var of global myvars {
    local graphlist "`graphlist' sens_hist_clean_`var'"
}

graph combine `graphlist', ///
    title("Distributions of All Variables by Session and Group (Sensitivity)") ///
    b1title("Variable Value") ///
    l1title("Percent") ///
    name(sens_combined_hists, replace)
	
graph export "stata/plots/sens/sens_combined_hists.pdf", replace

// Reshape to wide for pre vs post plots and correlations
drop SessionGroup
reshape wide $myvars, i(PPID group) j(session) string
save "stata/data/data_beh_all_wide.dta", replace

// Pre vs post scatterplots by group
foreach var of global myvars {
    // Overall correlation
    pwcorr `var'pre `var'post, obs

    // By group
    foreach g in treatment control {
        di as text "Correlation for `var' in group: `g'"
        pwcorr `var'pre `var'post if group == "`g'", obs
    }

    // Individual scatterplots with labels
    twoway (scatter `var'post `var'pre if group == "treatment", mcolor(blue)) ///
           (scatter `var'post `var'pre if group == "control",  mcolor(red)), ///
           title("Pre vs Post for `var' (Sensitivity)") ///
           legend(order(1 "Treatment" 2 "Control")) ///
           xtitle("`var' Pre") ///
           ytitle("`var' Post") ///
           name(sens_scatter_`var', replace)
    graph export "stata/plots/sens/sens_scatter_`var'.pdf", replace

    // Clean version (for combined figure)
    twoway (scatter `var'post `var'pre if group == "treatment", mcolor(blue)) ///
           (scatter `var'post `var'pre if group == "control",  mcolor(red)), ///
           title("`var'") ///
		   legend(off) ///
		   xtitle("") ///
		   ytitle("") ///
           name(sens_scatter_clean_`var', replace)
}

// Combined scatterplots
local graphlist ""

foreach var of global myvars {
    local graphlist "`graphlist' sens_scatter_clean_`var'"
}

graph combine `graphlist', ///
    title("Pre vs Post Correlations - All Variables (Sensitivity)") ///
    b1title("Pre-treatment Value") ///
    l1title("Post-treatment Value") ///
    commonscheme ///
    note("Blue = Treatment, Red = Control") ///
    name(sens_combined_scatter, replace)	
	
graph export "stata/plots/sens/sens_combined_scatterplots.pdf", replace



////////////////////////////////////////////////////////////////////////////////
// Individual Trajectories /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

use "stata/data/data_beh_all.dta", clear

// Numeric session (0 = pre, 1 = post)
gen session_num = .
replace session_num = 0 if session == "pre"
replace session_num = 1 if session == "post"
label define session_lbl 0 "Pre" 1 "Post"
label values session_num session_lbl

// Numeric group (0 = control, 1 = treatment)
gen g = cond(group=="control", 0, 1)
label define g_lbl 0 "control" 1 "treatment"
label values g g_lbl

// Numeric participant ID
encode PPID, gen(ppid_num)
xtset ppid_num session_num

// Individual trajectories (spaghetti plots)
foreach var of global myvars {
    twoway (scatter `var' session_num if group=="treatment", ///
           s(i) connect(L) msymbol(i) mcolor(blue%50) lcolor(blue%50)) ///
       (scatter `var' session_num if group=="control", ///
           s(i) connect(L) msymbol(i) mcolor(red%50) lcolor(red%50)), ///
       legend(order(1 "Treatment" 2 "Control")) ///
       xtitle("Session") ///
       xlabel(0 "Pre" 1 "Post", valuelabel) ///
       xscale(range(-0.2 1.2)) ///
       ytitle("`var'") ///
       title("Individual Trajectories: `var' (Sensitivity)") ///
       name(sens_spaghetti_`var', replace) ///
       aspectratio(1)
    graph export "stata/plots/sens/sens_spaghetti_`var'.pdf", replace
}

// Combined trajectories
local graphlist ""

foreach var of global myvars {
    local graphlist "`graphlist' sens_spaghetti_`var'"
}

graph combine `graphlist', ///
    title("Individual Trajectories of All Variables (Sensitivity)") ///
    name(sens_combined_spaghetti, replace)
	
graph export "stata/plots/sens/sens_combined_spaghetti.pdf", replace



////////////////////////////////////////////////////////////////////////////////
// Linear Mixed Effects Models /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Outcomes: STAI_T, GAD7, muTheta, muPhi, AvgCertainty
// Fixed effects: session, group, interaction
// Random effects: intercept

use "stata/data/data_beh_all.dta", clear

// Recreate numeric variables
capture drop session_num group_num ppid_num
gen session_num = cond(session == "pre", 0, 1)
encode group, gen(group_num)
tab group group_num // coding check
encode PPID, gen(ppid_num)
xtset ppid_num session_num

// Define outcomes
local outcomes STAI_T GAD7 muTheta muPhi AvgCertainty

// Excel setup
putexcel set "stata/results/sens/sens_mixed_models_both_results.xlsx", replace

putexcel A1 = "Variable"
putexcel B1 = "LRT_ChiSq"
putexcel C1 = "LRT_pvalue"
putexcel D1 = "Treatment_Intercept"
putexcel E1 = "Treatment_Int_pvalue"
putexcel F1 = "Control_Intercept"
putexcel G1 = "Control_Int_pvalue"
putexcel H1 = "Treatment_Slope"
putexcel I1 = "Treatment_Slope_pvalue"
putexcel J1 = "Control_Slope"
putexcel K1 = "Control_Slope_pvalue"
putexcel L1 = "Group_Slope_Diff"
putexcel M1 = "Group_Slope_pvalue"
putexcel N1 = "Group_Slope_LL"
putexcel O1 = "Group_Slope_UL"
putexcel P1 = "Std_Group_Slope_Diff"
putexcel Q1 = "Std_Group_Slope_LL"
putexcel R1 = "Std_Group_Slope_UL"
putexcel S1 = "Within_Control_PostMinusPre"
putexcel T1 = "Within_Control_p"
putexcel U1 = "Within_Control_LL"
putexcel V1 = "Within_Control_UL"
putexcel W1 = "Within_Treatment_PostMinusPre"
putexcel X1 = "Within_Treatment_p"
putexcel Y1 = "Within_Treatment_LL"
putexcel Z1 = "Within_Treatment_UL"
putexcel AA1 = "Between_Pre_TreatMinusControl"
putexcel AB1 = "Between_Pre_p"
putexcel AC1 = "Between_Pre_LL"
putexcel AD1 = "Between_Pre_UL"
putexcel AE1 = "Between_Post_TreatMinusControl"
putexcel AF1 = "Between_Post_p"
putexcel AG1 = "Between_Post_LL"
putexcel AH1 = "Between_Post_UL"

local row = 2

// Loop through each outcome
foreach var of local outcomes {
    di "{res}{bf}Working on variable: `var'"

////// 1. IDENTIFY PREFERRED MODEL /////////////////////////////////////////////
    
	// Fit random intercept only model
    quietly mixed `var' c.session_num##i.group_num || ppid_num:, reml
    estimates store int_only
	
	// Fit random intercept and slope model
	quietly mixed `var' c.session_num##i.group_num || ppid_num: session_num, ///
	covariance(unstructured) reml
	
    estimates store int_slope
	
	// Compare models using LRT
	di "{res}{bf}Likelihood-Ratio Test for: `var'"
    lrtest int_only int_slope
	local lrt_chi2 = r(chi2)
	local lrt_p = r(p)
	
	// Identify preferred model
	if `lrt_p' < 0.05 {
		local final_model "int_slope"
		di "{res}{bf} Using random intercept and slope model for `var'"
	}
	
	else {
		local final_model "int_only"
		di "{res}{bf} Using random intercept only model for `var'"
	}
	
	estimates restore `final_model'
	
////// 2. DIAGNOSTICS OF PREFERRED MODEL ///////////////////////////////////////
	
	// Standardised residuals and fitted values
	predict double r_`var', rstandard
	predict double xb_`var', xb
	
	// Residuals: Q-Q plots
	qnorm r_`var', ///
        name(sens_qq_resid_`var', replace) ///
        title("Q–Q Plot of Residuals: `var' (Sensitivity)")
		
    graph export "stata/plots/sens/sens_qq_resid_`var'.pdf", replace
	
	// Residuals: Histograms
	histogram r_`var', normal ///
        title("Residual Distribution: `var' (Sensitivity)") ///
        name(sens_hist_resid_`var', replace)
		
    graph export "stata/plots/sens/sens_resid_hist_`var'.pdf", replace
	
	// Residuals: Homoscedasticity (residuals vs fitted)
	scatter r_`var' xb_`var', ///
        yline(0) ///
        title("Residuals vs Fitted: `var' (Sensitivity)") ///
        xtitle("Fitted values") ///
        ytitle("Standardised residuals") ///
        name(sens_resfit_`var', replace)
		
    graph export "stata/plots/sens/sens_resfit_`var'.pdf", replace
	
	// Residuals: Variance over time (residuals vs session)
	scatter r_`var' session_num, ///
		yline(0) ///
		title("Residuals vs Session: `var' (Sensitivity)") ///
		xtitle("Session number") ///
		ytitle("Standardised residuals") ///
		name(sens_restime_`var', replace)
		
	graph export "stata/plots/sens/sens_restime_`var'.pdf", replace
	
	// Residuals: Outliers
	egen resid_z_`var' = std(r_`var')
    count if abs(resid_z_`var') > 3
    local n_outliers = r(N)
    di "{res}{bf}N(|standardised residual| > 3) for `var': " `n_outliers'
	
    if `n_outliers' > 0 {
        di "{res}{bf}Residual outliers for `var':"
        list PPID session group `var' r_`var' resid_z_`var' ///
            if abs(resid_z_`var') > 3, noobs
        preserve
            keep PPID session group ppid_num session_num ///
                 `var' r_`var' resid_z_`var'
            keep if abs(resid_z_`var') > 3
            export delimited using "stata/outliers/sens/sens_outliers_`var'.csv", replace
        restore
    }
	
    drop resid_z_`var'
	
	// ICC (correlation between repeated measures)
    estat icc
    local icc = r(icc2)
    
////// 3. RANDOM-EFFECTS DIAGNOSTICS ///////////////////////////////////////////
    
	di "{res}{bf}Random-effects diagnostics for: `var'"
    
    // Random-effects covariance (intercept/slope, if present)
    estat recov

    // Predict random effects for each participant
    if "`final_model'" == "int_slope" {
		
        // Random intercept and slope
        capture drop re_int_`var' re_slope_`var'
        predict double re_int_`var' re_slope_`var', reffects
        
        // Q–Q plot and histogram for random intercepts
        qnorm re_int_`var', ///
            name(sens_qq_reint_`var', replace) ///
            title("Q–Q Plot: Random Intercepts (`var') (Sensitivity)")
			
        graph export "stata/plots/sens/sens_qq_reint_`var'.pdf", replace
        
        histogram re_int_`var', normal ///
            title("Random Intercept Distribution: `var' (Sensitivity)") ///
            name(sens_hist_reint_`var', replace)
			
        graph export "stata/plots/sens/sens_hist_reint_`var'.pdf", replace
        
        // Q–Q plot and histogram for random slopes
        qnorm re_slope_`var', ///
            name(sens_qq_reslope_`var', replace) ///
            title("Q–Q Plot: Random Slopes (`var') (Sensitivity)")
			
        graph export "stata/plots/sens/sens_qq_reslope_`var'.pdf", replace
        
        histogram re_slope_`var', normal ///
            title("Random Slope Distribution: `var' (Sensitivity)") ///
            name(sens_hist_reslope_`var', replace)
			
        graph export "stata/plots/sens/sens_hist_reslope_`var'.pdf", replace
        
        // Scatterplot of random slopes vs random intercepts
        scatter re_slope_`var' re_int_`var', ///
            yline(0) xline(0) ///
            xtitle("Random Intercept") ///
            ytitle("Random Slope") ///
            title("Random Intercept vs Slope: `var' (Sensitivity)") ///
            name(sens_scatter_re_`var', replace)
			
        graph export "stata/plots/sens/sens_scatter_re_`var'.pdf", replace
    }
    
	else {
		
        // Random intercept only
        capture drop re_int_`var'
        predict double re_int_`var', reffects
        
        // Q–Q plot and histogram for random intercepts
        qnorm re_int_`var', ///
            name(sens_qq_reint_`var', replace) ///
            title("Q–Q Plot: Random Intercepts (`var') (Sensitivity)")
			
        graph export "stata/plots/sens/sens_qq_reint_`var'.pdf", replace
        
        histogram re_int_`var', normal ///
            title("Random Intercept Distribution: `var' (Sensitivity)") ///
            name(sens_hist_reint_`var', replace)
			
        graph export "stata/plots/sens/sens_hist_reint_`var'.pdf", replace
    }
	
////// 4. FIXED-EFFECTS SUMMARIES OF PREFERRED MODEL ///////////////////////////
	
	// Control intercept 
	di "{res}{bf}Intercepts (session = 0/pre):"
	lincom _cons // control
	local control_int = r(estimate)
	local control_int_p = r(p)
	
	// Treatment intercept
	lincom _cons + 2.group_num // treatment (ref + group effect)
	local treat_int = r(estimate)
	local treat_int_p = r(p)
	
	// Baseline group difference via group term
	lincom 2.group_num
	local base_diff = r(estimate)
	local base_diff_p  = r(p)
	local base_diff_ll = r(lb)
	local base_diff_ul = r(ub)
	
	// Control slope
	di "{res}{bf}Slopes (change per session_num):"
	lincom session_num 
	local control_slope = r(estimate)
	local control_slope_se = r(se)
	local control_slope_p = r(p)
	local control_slope_ll = r(lb)
	local control_slope_ul = r(ub)
	
	// Treatment slope
	lincom session_num + 2.group_num#c.session_num 
	local treat_slope    = r(estimate)
	local treat_slope_se = r(se)
	local treat_slope_p  = r(p)
	local treat_slope_ll = r(lb)
	local treat_slope_ul = r(ub)
	
	// Group x time interaction (difference in slopes)
	di "{res}{bf}Slope comparison (group × session) for: `var'"
	lincom 2.group_num#c.session_num
	local group_slope_diff    = r(estimate)
	local group_slope_diff_se = r(se)
	local group_slope_p       = r(p)
	local group_slope_diff_ll = r(lb)
	local group_slope_diff_ul = r(ub)
	
////// 5. LINEAR CONTRASTS /////////////////////////////////////////////////////

	// CONTROL: Within-group (group_num==1): change = session_num
	local within_control = `control_slope'
	local within_control_p = `control_slope_p'
	local within_control_ll = `control_slope_ll'
	local within_control_ul = `control_slope_ul'

	// TREATMENT: Within-group(group_num==2): change = session_num + interaction
	local within_treat = `treat_slope'
	local within_treat_p = `treat_slope_p'
	local within_treat_ll = `treat_slope_ll'
	local within_treat_ul = `treat_slope_ul'

	// PRE: Between-group (session_num==0): difference = group main effect
	local between_pre = `base_diff'
	local between_pre_p = `base_diff_p'
	local between_pre_ll = `base_diff_ll'
	local between_pre_ul = `base_diff_ul'

	// POST: Between-group (session_num==1): difference = group + interaction
	lincom 2.group_num + 2.group_num#c.session_num
	local between_post = r(estimate)
	local between_post_p = r(p)
	local between_post_ll = r(lb)
	local between_post_ul = r(ub)
	
////// 6. MARGINS & PLOTS //////////////////////////////////////////////////////
	
	margins group_num, at(session_num=(0 1))
	
	marginsplot, ///
		plot1opts(lcolor(navy)) ///
		plot2opts(lcolor(maroon)) ///
		title("Predicted `var' by Group and Session (Sensitivity)") ///
        name(sens_mplot_`var', replace) ///
		ytitle("Predicted `var'") xtitle("Session")
		
    graph export "stata/plots/sens/sens_margins_`var'.pdf", replace
	
////// 7. EFFECT SIZES /////////////////////////////////////////////////////////
	
	// SD of outcome across all observations
	summ `var'
	local sd_y = r(sd)

	// SD of predictor session_num (0/1 coding; SD = 0.5)
	local sd_x = sqrt(.25)
	
	// Scaling factor for standardisation
	local k = `sd_x' / `sd_y'
	
	// Standardised group x time interaction + CIs
	local group_slope_diff_std    = `group_slope_diff' * `k'
	local group_slope_diff_std_ll = `group_slope_diff_ll' * `k'
	local group_slope_diff_std_ul = `group_slope_diff_ul' * `k'

////// 8. EXPORT RESULTS ///////////////////////////////////////////////////////
    
	putexcel A`row' = "`var'" ///
        B`row' = `lrt_chi2' ///
        C`row' = `lrt_p' ///
        D`row' = `treat_int' ///
        E`row' = `treat_int_p' ///
        F`row' = `control_int' ///
        G`row' = `control_int_p' ///
        H`row' = `treat_slope' ///
        I`row' = `treat_slope_p' ///
        J`row' = `control_slope' ///
        K`row' = `control_slope_p' ///
        L`row' = `group_slope_diff' ///
        M`row' = `group_slope_p' ///
		N`row' = `group_slope_diff_ll' ///
		O`row' = `group_slope_diff_ul' ///
        P`row' = `group_slope_diff_std' ///
        Q`row' = `group_slope_diff_std_ll' ///
        R`row' = `group_slope_diff_std_ul' ///
		S`row' = `within_control' ///
		T`row' = `within_control_p' ///
		U`row' = `within_control_ll' ///
		V`row' = `within_control_ul' ///
		W`row' = `within_treat' ///
		X`row' = `within_treat_p' ///
		Y`row' = `within_treat_ll' ///
		Z`row' = `within_treat_ul' ///
		AA`row' = `between_pre' ///
		AB`row' = `between_pre_p' ///
		AC`row' = `between_pre_ll' ///
		AD`row' = `between_pre_ul' ///
		AE`row' = `between_post' ///
		AF`row' = `between_post_p' ///
		AG`row' = `between_post_ll' ///
		AH`row' = `between_post_ul'

    // Drop stored models
    estimates drop int_only
    estimates drop int_slope
    
    local row = `row' + 1
}

di "{res}{bf}Results exported to sens_mixed_models_both_results.xlsx"

// Combine margin plots 
local graphlist ""

foreach var of local outcomes {
    local graphlist "`graphlist' sens_mplot_`var'"
}

graph combine `graphlist', ///
    title("All Predicted Variable Values (Sensitivity)") ///
    name(sens_combined_margins, replace)
	
graph export "stata/plots/sens/sens_combined_margins.pdf", replace

save "stata/data/blt_beh_all.dta", replace



////////////////////////////////////////////////////////////////////////////////
// Correlations: Changes in Learning Rates, Questionnaires, Certainties ////////
////////////////////////////////////////////////////////////////////////////////

// Ensure data is reshaped to wide format, specify only treatment group
use "stata/data/data_beh_all_wide.dta", clear
keep if group == "treatment"

// Calculate change scores
gen gad7_change = GAD7post - GAD7pre
gen stait_change = STAI_Tpost - STAI_Tpre
gen certainty_change = AvgCertaintypost - AvgCertaintypre
gen muTheta_change = muThetapost - muThetapre

// Set up matrix 
tempname corr_results
matrix `corr_results' = J(4,3,.)

matrix rownames `corr_results' = gad7_change__muTheta ///
	stait_change__muTheta ///
	gad7_change__certainty ///
	stait_change__certainty
	
matrix colnames `corr_results' = r p N

local row = 1

// 1. CALCULATE CORRELATIONS ///////////////////////////////////////////////////

// GAD-7 change vs muTheta change
quietly pwcorr gad7_change muTheta_change, sig obs
matrix `corr_results'[`row',1] = r(C)[1,2]
matrix `corr_results'[`row',2] = r(sig)[1,2]
matrix `corr_results'[`row',3] = r(N)
local row = `row' + 1

// STAI-T change vs muTheta change
quietly pwcorr stait_change muTheta_change, sig obs
matrix `corr_results'[`row',1] = r(C)[1,2]
matrix `corr_results'[`row',2] = r(sig)[1,2]
matrix `corr_results'[`row',3] = r(N)
local row = `row' + 1

// GAD-7 change vs certainty change
quietly pwcorr gad7_change certainty_change, sig obs
matrix `corr_results'[`row',1] = r(C)[1,2]
matrix `corr_results'[`row',2] = r(sig)[1,2]
matrix `corr_results'[`row',3] = r(N)
local row = `row' + 1

// STAI-T change vs certainty change
quietly pwcorr stait_change certainty_change, sig obs
matrix `corr_results'[`row',1] = r(C)[1,2]
matrix `corr_results'[`row',2] = r(sig)[1,2]
matrix `corr_results'[`row',3] = r(N)

matrix list `corr_results'

// Bonferroni for questionnaires vs muTheta
local p_gad7_muTheta = `corr_results'[1,2]
local p_stait_muTheta = `corr_results'[2,2]
local p_gad7_muTheta_bonf = `p_gad7_muTheta'  * 2
local p_stait_muTheta_bonf = `p_stait_muTheta' * 2

// Bonferroni for questionnaires vs certainty
local p_certainty_gad7 = `corr_results'[3,2]
local p_certainty_stait = `corr_results'[4,2]
local p_certainty_gad7_bonf = `p_certainty_gad7'  * 2
local p_certainty_stait_bonf = `p_certainty_stait' * 2

// 2. FISHER TESTS  ////////////////////////////////////////////////////////////

// Learning rate
// Extract correlations
local r_mutheta_gad7 = `corr_results'[1,1]
local r_mutheta_stait = `corr_results'[2,1] 

// Extract N
local N_mutheta = `corr_results'[1,3]

// Fisher r-to-z test (independent)
local z_mutheta_gad7 = 0.5 * ln((1 + `r_mutheta_gad7') / (1 - `r_mutheta_gad7'))
local z_mutheta_stait = 0.5 * ln((1 + `r_mutheta_stait') / (1 - `r_mutheta_stait'))
local se_ind_mutheta = sqrt(2/(`N_mutheta' - 3))
local z_ind_mutheta = (`z_mutheta_stait' - `z_mutheta_gad7') / `se_ind_mutheta'
local p_ind_mutheta = 2 * (1 - normal(abs(`z_ind_mutheta')))

// Steiger test (dependent)
quietly pwcorr gad7_change stait_change, sig obs
local r_gad7_stait = r(C)[1,2]

local numerator_mutheta = (`r_mutheta_gad7' - `r_mutheta_stait') * ///
    sqrt((`N_mutheta' - 1) * (1 + `r_gad7_stait'))

local denominator_mutheta = sqrt( ///
    2*(`N_mutheta' - 1)/(`N_mutheta' - 3) * (1 - `r_gad7_stait') * ///
    (1 + `r_gad7_stait') + ///
    ((`r_mutheta_gad7' + `r_mutheta_stait')^2/4) * (1 - `r_gad7_stait')^3 )

local z_steiger_mutheta = `numerator_mutheta' / `denominator_mutheta'
local p_steiger_mutheta = 2 * (1 - normal(abs(`z_steiger_mutheta')))

// Certainty
// Extract correlations
local r_certainty_gad7  = `corr_results'[3,1]
local r_certainty_stait = `corr_results'[4,1]

// Extract N
local N_certainty = `corr_results'[3,3]

// Fisher r-to-z test (independent)
local z_certainty_gad7  = 0.5 * ln((1 + `r_certainty_gad7')  / (1 - `r_certainty_gad7'))
local z_certainty_stait = 0.5 * ln((1 + `r_certainty_stait') / (1 - `r_certainty_stait'))
local se_ind_certainty  = sqrt(2/(`N_certainty' - 3))
local z_ind_certainty   = (`z_certainty_stait' - `z_certainty_gad7') / `se_ind_certainty'
local p_ind_certainty   = 2 * (1 - normal(abs(`z_ind_certainty')))

// Steiger test (dependent)
quietly pwcorr gad7_change stait_change, sig obs
local r_gad7_stait_certainty = r(C)[1,2]

local numerator_certainty   = (`r_certainty_gad7' - `r_certainty_stait') * ///
    sqrt((`N_certainty' - 1) * (1 + `r_gad7_stait_certainty'))

local denominator_certainty = sqrt( ///
    2*(`N_certainty' - 1)/(`N_certainty' - 3) * (1 - `r_gad7_stait_certainty') * ///
    (1 + `r_gad7_stait_certainty') + ///
    ((`r_certainty_gad7' + `r_certainty_stait')^2/4) * (1 - `r_gad7_stait_certainty')^3 )

local z_steiger_certainty = `numerator_certainty' / `denominator_certainty'
local p_steiger_certainty = 2 * (1 - normal(abs(`z_steiger_certainty')))

// 3. SCATTERPLOTS  ////////////////////////////////////////////////////////////

// GAD-7, muTheta
twoway ///
    (scatter gad7_change muTheta_change) ///
    (lfit    gad7_change muTheta_change), ///
    title("ΔGAD-7 vs ΔmuTheta (Sensitivity)") ///
    xtitle("ΔmuTheta (post - pre)") ///
    ytitle("ΔGAD-7 (post - pre)") ///
    name(sens_sc_muTheta_gad7, replace)
	
graph export "stata/plots/sens/sens_scatter_muTheta_gad7_change.pdf", replace

// STAI-T, muTheta
twoway ///
    (scatter stait_change muTheta_change) ///
    (lfit    stait_change muTheta_change), ///
    title("ΔSTAI-T vs ΔmuTheta (Sensitivity)") ///
    xtitle("ΔmuTheta (post - pre)") ///
    ytitle("ΔSTAI-T (post - pre)") ///
    name(sens_sc_muTheta_stait, replace)
	
graph export "stata/plots/sens/sens_scatter_muTheta_stait_change.pdf", replace

// Combine muTheta plots
graph combine sens_sc_muTheta_gad7 sens_sc_muTheta_stait, ///
    cols(2) ///
    name(sens_combined_muTheta_scatter, replace)
	
graph export "stata/plots/sens/sens_combined_muTheta_scatterplots.pdf", replace

// GAD-7, certainty
twoway ///
    (scatter gad7_change certainty_change) ///
    (lfit gad7_change certainty_change), ///
    title("ΔGAD-7 vs ΔCertainty (Sensitivity)") ///
    xtitle("ΔCertainty (post - pre)") ///
    ytitle("ΔGAD-7 (post - pre)") ///
    name(sens_sc_cert_gad7, replace)
	
graph export "stata/plots/sens/sens_scatter_certainty_gad7_change.pdf", replace

// STAI-T, certainty
twoway ///
    (scatter stait_change certainty_change) ///
    (lfit stait_change certainty_change), ///
    title("ΔSTAI-T vs ΔCertainty (Sensitivity)") ///
    xtitle("ΔCertainty (post - pre)") ///
    ytitle("ΔSTAI-T (post - pre)") ///
    name(sens_sc_cert_stait, replace)
	
graph export "stata/plots/sens/sens_scatter_certainty_stait_change.pdf", replace

// Combine certainty plots
graph combine sens_sc_cert_gad7 sens_sc_cert_stait, ///
    cols(2) ///
    name(sens_combined_certainty_scatter, replace)
	
graph export "stata/plots/sens/sens_combined_certainty_scatterplots.pdf", replace

// 4. EXPORT RESULTS ///////////////////////////////////////////////////////////

// Pull matrix values into locals for putexcel
local r_gad_mu = el(`corr_results', 1, 1)
local p_gad_mu = el(`corr_results', 1, 2)
local n_gad_mu = el(`corr_results', 1, 3)

local r_stai_mu = el(`corr_results', 2, 1)
local p_stai_mu = el(`corr_results', 2, 2)
local n_stai_mu = el(`corr_results', 2, 3)

local r_gad_cert = el(`corr_results', 3, 1)
local p_gad_cert = el(`corr_results', 3, 2)
local n_gad_cert = el(`corr_results', 3, 3)

local r_stai_cert = el(`corr_results', 4, 1)
local p_stai_cert = el(`corr_results', 4, 2)
local n_stai_cert = el(`corr_results', 4, 3)

// Export
putexcel set "stata/results/sens/sens_correlation_results.xlsx", replace

putexcel A1 = "Correlations using change scores (post - pre)"

putexcel A3 = "Association"
putexcel B3 = "r"
putexcel C3 = "p"
putexcel D3 = "N"
putexcel E3 = "p (Bonferroni, questionnaires with muTheta only)"

putexcel A4 = "GAD-7 vs muTheta" ///
    B4 = `r_gad_mu' ///
    C4 = `p_gad_mu' ///
    D4 = `n_gad_mu' ///
    E4 = `p_gad7_muTheta_bonf'

putexcel A5 = "STAI-T vs muTheta" ///
    B5 = `r_stai_mu' ///
    C5 = `p_stai_mu' ///
    D5 = `n_stai_mu' ///
    E5 = `p_stait_muTheta_bonf'

putexcel A6 = "GAD-7 vs Certainty" ///
    B6 = `r_gad_cert' ///
    C6 = `p_gad_cert' ///
    D6 = `n_gad_cert'

putexcel A7 = "STAI-T vs Certainty" ///
    B7 = `r_stai_cert' ///
    C7 = `p_stai_cert' ///
    D7 = `n_stai_cert'

putexcel A9  = "Comparison of r(muTheta, GAD-7) vs r(muTheta, STAI-T)"

putexcel A11 = "Fisher r-to-z test (independent correlations)"
putexcel A12 = "z"  B12 = `z_ind_mutheta'
putexcel A13 = "p"  B13 = `p_ind_mutheta'

putexcel A15 = "Steiger test (dependent correlations)"
putexcel A16 = "z"  B16 = `z_steiger_mutheta'
putexcel A17 = "p"  B17 = `p_steiger_mutheta'

putexcel A19 = "r(GAD-7, STAI-T) used in Steiger test"
putexcel B19 = `r_gad7_stait'
putexcel A20 = "N used in Fisher/Steiger"
putexcel B20 = `N_mutheta'

putexcel A22 = "Comparison of r(Certainty, GAD-7) vs r(Certainty, STAI-T)"

putexcel A24 = "Fisher r-to-z test (independent correlations)"
putexcel A25 = "z"  B25 = `z_ind_certainty'
putexcel A26 = "p"  B26 = `p_ind_certainty'

putexcel A28 = "Steiger test (dependent correlations)"
putexcel A29 = "z"  B29 = `z_steiger_certainty'
putexcel A30 = "p"  B30 = `p_steiger_certainty'

putexcel A32 = "r(GAD-7, STAI-T) used in Steiger test"
putexcel B32 = `r_gad7_stait_certainty'

putexcel A33 = "N used in Fisher/Steiger"
putexcel B33 = `N_certainty'

di as text "{bf}Correlation results exported to sens_correlation_results.xlsx"

set graphics on
