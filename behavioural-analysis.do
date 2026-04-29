////////////////////////////////////////////////////////////////////////////////
// DATA ANALYSIS OF BLT FROM BREATHING AND ANXIETY CLINICAL STUDY //////////////
////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------- //
// Author: Sophie Morris										     		  //
// Created: 23/07/2025														  //
// -------------------------------------------------------------------------- //
// Description: This script runs the analysis of BLT behavioural data for     //
// Sophie M's Masters' project. The data was collected as part of the         //
// breathing and anxiety study conducted collaboratively with the Department  //
// of Psychology, School of Sport, Physical Education and Sport Science, and  //
// the School of Pharmacy at the University of Otago. Data from both clinical //
// cohorts are used within this analysis. A timestamped analysis plan is      //
// available on	(https://github.com/IMAGEotago). 							  //
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
tempfile data_beh_included

import excel using "matlab/blt_beh_data_included.xlsx", sheet("Subjects") firstrow clear
gen group = "treatment"
replace group = "control" if substr(PPID, 1, 2) == "35"
save `data_beh_included'

import excel using "stata/data/results-clinQuest.xlsx", sheet("Sums_Data") firstrow clear
gen group = "treatment"
append using `data_beh_included'
save `data_beh_included', replace

import excel using "stata/data/results-clinControlQuest.xlsx", sheet("Sums_Data") firstrow clear
gen group = "control"
append using `data_beh_included'
save `data_beh_included', replace

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

save "stata/data/data_beh_included.dta", replace

// Variable list
global myvars muTheta muPhi AvgCertainty GAD7 STAI_T



////////////////////////////////////////////////////////////////////////////////
// Summary Measures ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

use "stata/data/data_beh_included.dta", clear

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

putexcel set "stata/results/summaryStat.xlsx", replace

putexcel set "stata/results/summaryStat.xlsx", sheet("post_control") modify
putexcel A1 = matrix(r(Stat1)), names

putexcel set "stata/results/summaryStat.xlsx", sheet("post_treatment") modify
putexcel A1 = matrix(r(Stat2)), names

putexcel set "stata/results/summaryStat.xlsx", sheet("pre_control") modify
putexcel A1 = matrix(r(Stat3)), names

putexcel set "stata/results/summaryStat.xlsx", sheet("pre_treatment") modify
putexcel A1 = matrix(r(Stat4)), names



////////////////////////////////////////////////////////////////////////////////
// Repeated Measures Observations //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Histograms by group and session
foreach var of global myvars {
    histogram `var', ///
		percent bin(20) ///
        by(session group, col(2) note("") ///
        title("Distribution of `var' by Group & Session")) ///
        ytitle("Percent") xtitle("Value") ///
        name(hist_`var', replace)
    graph export "stata/plots/hist_`var'.pdf", replace
}

// Clean histograms for combined figure
foreach var of global myvars {
    histogram `var', ///
		percent bin(10) ///
        by(session group, col(2) note("") ///
		title("`var'")) ///
        ytitle("") ///
		xtitle("") ///
        name(hist_clean_`var', replace)
}

local graphlist ""

foreach var of global myvars {
    local graphlist "`graphlist' hist_clean_`var'"
}

graph combine `graphlist', ///
    title("Distributions of All Variables by Session and Group") ///
    b1title("Variable Value") ///
    l1title("Percent") ///
    name(combined_hists, replace)
	
graph export "stata/plots/combined_hists.pdf", replace

// Reshape to wide for pre vs post plots and correlations
drop SessionGroup
reshape wide $myvars, i(PPID group) j(session) string
save "stata/data/data_beh_included_wide.dta", replace

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
           title("Pre vs Post for `var'") ///
           legend(order(1 "Treatment" 2 "Control")) ///
           xtitle("`var' Pre") ///
           ytitle("`var' Post") ///
           name(scatter_`var', replace)
		   
    graph export "stata/plots/scatter_`var'.pdf", replace

    // Clean version (for combined figure)
    twoway (scatter `var'post `var'pre if group == "treatment", mcolor(blue)) ///
           (scatter `var'post `var'pre if group == "control",  mcolor(red)), ///
           title("`var'") ///
		   legend(off) ///
		   xtitle("") ///
		   ytitle("") ///
           name(scatter_clean_`var', replace)
}

// Combined scatterplots
local graphlist ""

foreach var of global myvars {
    local graphlist "`graphlist' scatter_clean_`var'"
}

graph combine `graphlist', ///
    title("Pre vs Post Correlations - All Variables") ///
    b1title("Pre-treatment Value") ///
    l1title("Post-treatment Value") ///
    commonscheme ///
    note("Blue = Treatment, Red = Control") ///
    name(combined_scatter, replace)
	
graph export "stata/plots/combined_scatterplots.pdf", replace



////////////////////////////////////////////////////////////////////////////////
// Baseline Comparisons ////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

use "stata/data/data_beh_included_wide.dta", clear

local outcomes STAI_T GAD7 muTheta muPhi AvgCertainty

// Prepare Excel output
putexcel set "stata/results/baseline_group_differences.xlsx", replace

putexcel A1 = "Variable" ///
	B1 = "Mean_Control" ///
	C1 = "SD_Control" ///
	D1 = "Mean_Treatment" ///
	E1 = "SD_Treatment" ///
	F1 = "Difference (Treat - Control)" ///
	G1 = "t-statistic" ///
	H1 = "p-value" ///
	I1 = "Cohen_d" ///
	J1 = "Hedges_g" ///
	K1 = "Diff_LL" ///
	L1 = "Diff_UL"

local row = 2

// Run t-test comparisons
foreach var of local outcomes {
    di "{res}{bf}Baseline comparison for `var':"
    
	// T-test
	quietly ttest `var'pre, by(group)
	
	// Extract group stats
    local mean_control = r(mu_1)
    local sd_control = r(sd_1)
	local n_control = r(N_1)
    local mean_treat = r(mu_2)
    local sd_treat = r(sd_2)
	local n_treat = r(N_2)

    // Differences and t-test results
    local diff = `mean_treat' - `mean_control'
    local tstat = r(t)
    local pval = r(p)
	local diff_ll = -r(lb)
	local diff_ul = -r(ub)
	
	// Pooled SD
	local pooled_sd = sqrt(((`n_control' - 1)*`sd_control'^2 + ///
		(`n_treat' - 1)*`sd_treat'^2) / ///
		(`n_control' + `n_treat' - 2))

    // Cohen's d
    if `pooled_sd' != 0 {
        local d = (`diff') / `pooled_sd'
    }
	
    else local d = .

    // Hedges' g (small-sample corrected d)
    local J = 1 - 3/(4*(`n_control' + `n_treat' - 2) - 1)
    local g = `d' * `J'
	
	// Export to Excel
    putexcel A`row' = "`var' (Pre)" ///
		B`row' = `mean_control' ///
		C`row' = `sd_control' ///
		D`row' = `mean_treat' ///
		E`row' = `sd_treat' ///
		F`row' = `diff' ///
		G`row' = `tstat' ///
		H`row' = `pval' ///
		I`row' = `d' ///
		J`row' = `g' ///
		K`row' = `diff_ll' ///
        L`row' = `diff_ul'

    local row = `row' + 1
}

di "{res}{bf}Baseline comparison results exported to baseline_group_differences.xlsx"



////////////////////////////////////////////////////////////////////////////////
// Individual Trajectories /////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

use "stata/data/data_beh_included.dta", clear

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
       title("Individual Trajectories: `var'") ///
       name(spaghetti_`var', replace) ///
       aspectratio(1)
    graph export "stata/plots/spaghetti_`var'.pdf", replace
}

// Combined trajectories
local graphlist ""

foreach var of global myvars {
    local graphlist "`graphlist' spaghetti_`var'"
}

graph combine `graphlist', ///
    title("Individual Trajectories of All Variables") ///
    name(combined_spaghetti, replace)
	
graph export "stata/plots/combined_spaghetti.pdf", replace



////////////////////////////////////////////////////////////////////////////////
// Linear Mixed Effects Models /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// Outcomes: STAI_T, GAD7, muTheta, muPhi, AvgCertainty
// Fixed effects: session, group, interaction
// Random effects: intercept, slope

use "stata/data/data_beh_included.dta", clear

// Recreate numeric variables
capture drop session_num group_num ppid_num
gen session_num = cond(session == "pre", 0, 1)
encode group, gen(group_num)
tab group group_num   // coding check
encode PPID, gen(ppid_num)
xtset ppid_num session_num

// Define outcomes
local outcomes STAI_T GAD7 muTheta muPhi AvgCertainty

// Excel setup
putexcel set "stata/results/mixed_models_both_results.xlsx", replace

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
        name(qq_resid_`var', replace) ///
        title("Q–Q Plot of Residuals: `var'")
		
    graph export "stata/plots/qq_resid_`var'.pdf", replace
	
	// Residuals: Histograms
	histogram r_`var', normal ///
        title("Residual Distribution: `var'") ///
        name(hist_resid_`var', replace)
		
    graph export "stata/plots/resid_hist_`var'.pdf", replace
	
	// Residuals: Homoscedasticity (residuals vs fitted)
	scatter r_`var' xb_`var', ///
        yline(0) ///
        title("Residuals vs Fitted: `var'") ///
        xtitle("Fitted values") ///
        ytitle("Standardised residuals") ///
        name(resfit_`var', replace)
		
    graph export "stata/plots/resfit_`var'.pdf", replace
	
	// Residuals: Variance over time (residuals vs session)
	scatter r_`var' session_num, ///
		yline(0) ///
		title("Residuals vs Session: `var'") ///
		xtitle("Session number") ///
		ytitle("Standardised residuals") ///
		name(restime_`var', replace)
		
	graph export "stata/plots/restime_`var'.pdf", replace
	
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
            export delimited using "stata/outliers/outliers_`var'.csv", replace
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
            name(qq_reint_`var', replace) ///
            title("Q–Q Plot: Random Intercepts (`var')")
			
        graph export "stata/plots/qq_reint_`var'.pdf", replace
        
        histogram re_int_`var', normal ///
            title("Random Intercept Distribution: `var'") ///
            name(hist_reint_`var', replace)
			
        graph export "stata/plots/hist_reint_`var'.pdf", replace
        
        // Q–Q plot and histogram for random slopes
        qnorm re_slope_`var', ///
            name(qq_reslope_`var', replace) ///
            title("Q–Q Plot: Random Slopes (`var')")
			
        graph export "stata/plots/qq_reslope_`var'.pdf", replace
        
        histogram re_slope_`var', normal ///
            title("Random Slope Distribution: `var'") ///
            name(hist_reslope_`var', replace)
			
        graph export "stata/plots/hist_reslope_`var'.pdf", replace
        
        // Scatterplot of random slopes vs random intercepts
        scatter re_slope_`var' re_int_`var', ///
            yline(0) xline(0) ///
            xtitle("Random Intercept") ///
            ytitle("Random Slope") ///
            title("Random Intercept vs Slope: `var'") ///
            name(scatter_re_`var', replace)
			
        graph export "stata/plots/scatter_re_`var'.pdf", replace
    }
    
	else {
        
		// Random intercept only
        capture drop re_int_`var'
        predict double re_int_`var', reffects
        
        // Q–Q plot and histogram for random intercepts
        qnorm re_int_`var', ///
            name(qq_reint_`var', replace) ///
            title("Q–Q Plot: Random Intercepts (`var')")
			
        graph export "stata/plots/qq_reint_`var'.pdf", replace
        
        histogram re_int_`var', normal ///
            title("Random Intercept Distribution: `var'") ///
            name(hist_reint_`var', replace)
			
        graph export "stata/plots/hist_reint_`var'.pdf", replace
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
		title("Predicted `var' by Group and Session") ///
        name(mplot_`var', replace) ///
		ytitle("Predicted `var'") xtitle("Session")
		
    graph export "stata/plots/margins_`var'.pdf", replace
	
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

di "{res}{bf}Results exported to mixed_models_both_results.xlsx"

// Combine margin plots 
local graphlist ""

foreach var of local outcomes {
    local graphlist "`graphlist' mplot_`var'"
}

graph combine `graphlist', ///
    title("All Predicted Variable Values") ///
    name(combined_margins, replace)
	
graph export "stata/plots/combined_margins.pdf", replace

save "stata/data/blt_beh.dta", replace



////////////////////////////////////////////////////////////////////////////////
// Outlier Analysis: Residual-Based Sensitivity  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

use "stata/data/data_beh_included.dta", clear

// Recreate numeric session, group, and ID
capture drop session_num group_num ppid_num
gen session_num = cond(session == "pre", 0, 1)
encode group, gen(group_num) // 1 = control, 2 = treatment
encode PPID, gen(ppid_num)

// Threshold for outliers
local threshold = 3

// Set up excel file
putexcel set "stata/results/residual_outlier_sensitivity_results.xlsx", replace

putexcel A1 = "Variable" ///
	B1 = "n_outliers" ///
	C1 = "main_coef" ///
	D1 = "main_p" ///
	E1 = "sens_coef" ///
	F1 = "sens_p" ///
	G1 = "percent_change" ///
	H1 = "sig_change_flag"
	
local row = 2

foreach var of local outcomes {
	
	capture drop r_`var' outlier_`var'
	capture scalar drop r_`var'
	
    // Fit main model on all data
	quietly mixed `var' c.session_num##i.group_num || ppid_num:, reml
	
	// Interaction term: 2.group_num#c.session_num (treatment vs control × time)
    local main_int_coef = _b[2.group_num#c.session_num]
    local main_int_se = _se[2.group_num#c.session_num]
    local main_int_p = 2 * ttail(e(df_r), abs(`main_int_coef'/`main_int_se'))
	
	// Standardized residuals from main model
    quietly predict double r_`var', rstandard

    // Flag residual outliers
    quietly gen byte outlier_`var' = abs(r_`var') > `threshold'

    // Count outliers
    quietly count if outlier_`var' == 1
    local n_outliers = r(N)

    di "{res}{bf}Total residual outliers identified for `var' (|rstandard| > `threshold'): " `n_outliers'
	
	// Defaults in case of no outliers
	local sens_int_coef = `main_int_coef'
	local sens_int_p    = `main_int_p'
    local coef_change   = 0
    local sig_flag      "No outliers"

    if `n_outliers' > 0 {
        di "{res}{bf}Outlier cases for `var':"
        list PPID session group `var' r_`var' if outlier_`var' == 1, ///
            sepby(PPID) compress noobs

		di "{res}{bf}Main model (all data):"
		di "{res}{bf}Interaction coef (session_num × g): " %8.4f `main_int_coef'
		di "{res}{bf}p-value: " %8.4f `main_int_p'

		// Sensitivity model: exclude residual outliers 
		di "{res}{bf}Sensitivity model (excluding residual outliers):"

		quietly mixed `var' c.session_num##i.group_num || ppid_num: ///
            if outlier_`var' == 0, reml

		local sens_int_coef = _b[2.group_num#c.session_num]
        local sens_int_se = _se[2.group_num#c.session_num]
        local sens_int_p = 2 * ttail(e(df_r), abs(`sens_int_coef'/`sens_int_se'))

		di "{res}{bf}Interaction coef (no outliers): " %8.4f `sens_int_coef'
		di "{res}{bf}p-value (no outliers): " %8.4f `sens_int_p'

		// Summary of change
		di "{res}{bf}Sensitivity summary for `var': "

		if `main_int_coef' != 0 {
            local coef_change = (`sens_int_coef' - `main_int_coef') / ///
				abs(`main_int_coef') * 100
            di "{res}{bf}% change in interaction coef: " %6.2f `coef_change' "%"
        }
		
    else {
        di "{res}{bf}Main interaction coef is 0; cannot compute % change."
        local coef_change = .
		}

    // Significance change flag
	local sig_flag "No change"
    if (`main_int_p' < 0.05 & `sens_int_p' >= 0.05) | ///
       (`main_int_p' >= 0.05 & `sens_int_p' < 0.05) {
        di "{err}{bf}CRITICAL: Significance changes when outliers removed"
		local sig_flag "Significance changed"
    }

    if `coef_change' != . {
        if abs(`coef_change') < 10 {
            di "{res}{bf}Interpretation: Results robust to outliers (< 10% change)."
        }
        else if abs(`coef_change') < 25 {
            di "{res}{bf}Interpretation: Moderate sensitivity to outliers (10–25% change)."
        }
        else {
            di "{err}{bf}Interpretation: High sensitivity to outliers (> 25% change)."
            di "{err}{bf}Consider reporting both sets of results."
            }
        }
    }
	
    else {
        di "{res}{bf}No residual outliers detected for `var' (|rstandard| > `threshold')."
    }
	
	// Export summary for this variable
	putexcel A`row' = "`var'" ///
             B`row' = `n_outliers' ///
             C`row' = `main_int_coef' ///
             D`row' = `main_int_p' ///
             E`row' = `sens_int_coef' ///
             F`row' = `sens_int_p' ///
             G`row' = `coef_change' ///
             H`row' = "`sig_flag'"
	
	local row = `row' + 1
	
	// Clean up variables
    drop r_`var' outlier_`var'
}

di "{res}{bf}Residual outlier sensitivity results exported to residual_outlier_sensitivity_results.xlsx"



////////////////////////////////////////////////////////////////////////////////
// Correlations: Changes in Learning Rates, Questionnaires, Certainties ////////
////////////////////////////////////////////////////////////////////////////////

// Ensure data is reshaped to wide format, specify only treatment group
use "stata/data/data_beh_included_wide.dta", clear
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
local p_gad7_muTheta_bonf  = `p_gad7_muTheta'  * 2
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
local r_certainty_gad7 = `corr_results'[3,1]
local r_certainty_stait = `corr_results'[4,1]

// Extract N
local N_certainty = `corr_results'[3,3]

// Fisher r-to-z test (independent)
local z_certainty_gad7 = 0.5 * ln((1 + `r_certainty_gad7') / (1 - `r_certainty_gad7'))
local z_certainty_stait = 0.5 * ln((1 + `r_certainty_stait') / (1 - `r_certainty_stait'))
local se_ind_certainty = sqrt(2/(`N_certainty' - 3))
local z_ind_certainty = (`z_certainty_stait' - `z_certainty_gad7') / `se_ind_certainty'
local p_ind_certainty = 2 * (1 - normal(abs(`z_ind_certainty')))

// Steiger test (dependent)
quietly pwcorr gad7_change stait_change, sig obs
local r_gad7_stait_certainty = r(C)[1,2]

local numerator_certainty = (`r_certainty_gad7' - `r_certainty_stait') * ///
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
    title("ΔGAD-7 vs ΔmuTheta") ///
    xtitle("ΔmuTheta (post - pre)") ///
    ytitle("ΔGAD-7 (post - pre)") ///
    name(sc_muTheta_gad7, replace)
	
graph export "stata/plots/scatter_muTheta_gad7_change.pdf", replace

// STAI-T, muTheta
twoway ///
    (scatter stait_change muTheta_change) ///
    (lfit    stait_change muTheta_change), ///
    title("ΔSTAI-T vs ΔmuTheta") ///
    xtitle("ΔmuTheta (post - pre)") ///
    ytitle("ΔSTAI-T (post - pre)") ///
    name(sc_muTheta_stait, replace)
	
graph export "stata/plots/scatter_muTheta_stait_change.pdf", replace

// Combine muTheta plots
graph combine sc_muTheta_gad7 sc_muTheta_stait, ///
    cols(2) ///
    name(combined_muTheta_scatter, replace)
	
graph export "stata/plots/combined_muTheta_scatterplots.pdf", replace

// GAD-7, certainty
twoway ///
    (scatter gad7_change certainty_change) ///
    (lfit gad7_change certainty_change), ///
    title("ΔGAD-7 vs ΔCertainty") ///
    xtitle("ΔCertainty (post - pre)") ///
    ytitle("ΔGAD-7 (post - pre)") ///
    name(sc_cert_gad7, replace)
	
graph export "stata/plots/scatter_certainty_gad7_change.pdf", replace

// STAI-T, certainty
twoway ///
    (scatter stait_change certainty_change) ///
    (lfit stait_change certainty_change), ///
    title("ΔSTAI-T vs ΔCertainty") ///
    xtitle("ΔCertainty (post - pre)") ///
    ytitle("ΔSTAI-T (post - pre)") ///
    name(sc_cert_stait, replace)
	
graph export "stata/plots/scatter_certainty_stait_change.pdf", replace

// Combine certainty plots
graph combine sc_cert_gad7 sc_cert_stait, ///
    cols(2) ///
    name(combined_certainty_scatter, replace)
	
graph export "stata/plots/combined_certainty_scatterplots.pdf", replace

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
putexcel set "stata/results/correlation_results.xlsx", replace

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

di as text "{bf}Correlation results exported to correlation_results.xlsx"

set graphics on
