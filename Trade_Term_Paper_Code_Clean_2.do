********************************************************************************
* Brazil VAR Analysis with Uncertainty Shocks
********************************************************************************

clear all
set more off
capture log close
set scheme s2color

* Install required packages
capture ssc install coefplot, replace
if _rc {
    di as error "Warning: Could not install coefplot package. Will use alternative plotting method."
    global use_alt_plots = 1
}
else {
    global use_alt_plots = 0
}

* Start log file
log using "var_analysis_detailed_log.txt", replace

* Load the data
use "C:\Users\azzam\OneDrive\Documents\PhD semesters\4 Spring 2025\Trade\Final Paper\Brazil_with_4Uncertainty.dta", clear

* Define time series
gen date_str = time
gen year = substr(time, 1, 4)
gen month = substr(time, 6, 2)
destring year month, replace
gen time_num = ym(year, month)
format time_num %tm
tsset time_num

* Create BRICS dummy (Brazil joined in December 2010)
gen brics_dummy = (time_num >= ym(2009, 6))
label var brics_dummy "BRICS membership dummy (=1 after June 2009)"

********************************************************************************
* 1. Variable Definitions and Standardization
********************************************************************************

* Define macro variables
global macro_vars "cpi ipi reer tradebalance stockmarket"
global all_vars "vix gepu_ppp tpu opu_index $macro_vars"



* Standardize all variables (mean 0, sd 1)
foreach var in $all_vars {
    egen `var'_std = std(`var')
    label var `var'_std "Standardized `var'"
}

* Create first differences of standardized variables
foreach var in $all_vars {
    gen d_`var' = D.`var'_std
    label var d_`var' "First Difference of Standardized `var'"
}

* Create interaction terms with BRICS dummy
foreach var in vix gepu_ppp tpu opu_index {
    gen brics_`var' = brics_dummy * d_`var'
    label var brics_`var' "Interaction: BRICS membership * First Diff of `var'"
}

* Check stationarity of differenced variables
foreach var in $all_vars {
    dfuller d_`var', lags(12)
}


* Graphing gepu and tpu
tsline gepu_ppp tpu, ///
    lcolor(red blue) ///
    legend(order(1 "GEPU" 2 "TPU")) ///
    ytitle("Index") ///
    xtitle("Date")

* cpi already inflation
	
********************************************************************************
* 2. GEPU Analysis with BRICS Interaction
********************************************************************************

* Determine optimal lag length
*varsoc d_gepu_ppp brics_gepu_ppp d_cpi d_ipi d_reer d_tradebalance d_stockmarket, maxlag(12)
varsoc gepu d_ipi d_reer d_cpi d_tradebalance d_stockmarket, maxlag(12)

scalar optimal_lags_gepu = e(ic_bic)
if (optimal_lags_gepu < 1) | (optimal_lags_gepu > 12) {
    scalar optimal_lags_gepu = 4
}
*var d_gepu_ppp d_cpi d_ipi d_reer d_tradebalance d_stockmarket, lags(1/4)
var gepu d_ipi d_reer d_cpi d_tradebalance d_stockmarket, lags(1/4)


* Generate IRFs with bootstrap confidence intervals - Pre-BRICS period
irf create gepu_pre, set(myirf2, replace) step(36) bs reps(50)
irf set myirf2
irf graph oirf, impulse(gepu) response(gepu d_ipi d_reer d_cpi d_tradebalance d_stockmarket) ///
    ci level(68) lstep(0) ustep(36) ///
    title("Responses to GEPU Shock", size(medium)) ///
    note("68% Bootstrap CI, Horizon: 36 months", size(small)) ///
    xlabel(0(6)36, labsize(small)) ylabel(,grid labsize(small)) ///
    graphregion(color(white)) bgcolor(white) ///
    byopts(rows(3) yrescale xrescale compact) ///
    name(combined_gepu_pre, replace) ysize(35) xsize(40)
graph export "irf_gepu_pre_brics.pdf", replace as(pdf)

* Forecast Error Variance Decomposition (%)
irf table fevd, impulse(gepu) response(d_ipi d_reer d_cpi d_tradebalance d_stockmarket) 


* Estimate VAR with GEPU and BRICS interaction
gen gepu_brics = gepu_ppp * brics_dummy

*varsoc d_gepu_ppp brics_gepu_ppp d_cpi d_ipi d_reer d_tradebalance d_stockmarket, maxlag(12)
*var d_gepu_ppp brics_gepu_ppp d_cpi d_ipi d_reer d_tradebalance d_stockmarket, lags(1/4)

varsoc gepu_ppp gepu_brics d_ipi d_reer d_cpi d_tradebalance d_stockmarket, maxlag(12)

var d_ln_gepu gepu_brics d_ipi d_reer d_cpi d_tradebalance d_stockmarket, lags(1/4)

estimates store var_gepu_brics

* Generate IRFs for BRICS interaction effect
irf create gepu_post, set(myirf2b, replace) step(36) bs reps(50)
irf set myirf2b
irf graph oirf, impulse(gepu_brics) response(gepu d_ipi d_reer d_cpi d_tradebalance d_stockmarket) ///
    ci level(68) lstep(0) ustep(36) ///
    title("Add. Effect of GEPU Shock post-BRICS Period", size(medium)) ///
    note("68% Bootstrap CI, Horizon: 36 months", size(small)) ///
    xlabel(0(6)36, labsize(small)) ylabel(,grid labsize(small)) ///
    graphregion(color(white)) bgcolor(white) ///
    byopts(rows(3) yrescale xrescale compact) ///
    name(combined_gepu_brics, replace) ysize(35) xsize(40)
graph export "irf_gepu_brics_effect.pdf", replace as(pdf)



********************************************************************************
* 3. TPU Analysis with BRICS Interaction
********************************************************************************
* Determine optimal lag length
* varsoc d_tpu brics_tpu d_cpi d_ipi d_reer d_tradebalance d_stockmarket, maxlag(12)
varsoc tpu d_ipi d_reer d_cpi d_tradebalance d_stockmarket, maxlag(12)

var tpu d_ipi d_reer d_cpi d_tradebalance d_stockmarket, lags(1/4)

scalar optimal_lags_tpu = e(ic_bic)
if (optimal_lags_tpu < 1) | (optimal_lags_tpu > 12) {
    scalar optimal_lags_tpu = 4
}


* Generate IRFs with bootstrap confidence intervals - Pre-BRICS period
irf create tpu_pre, set(myirf3, replace) step(36) bs reps(50)
irf set myirf3
irf graph oirf, impulse(tpu) response(tpu d_ipi d_reer d_cpi d_tradebalance d_stockmarket) ///
    ci level(68) lstep(0) ustep(36) ///
    title("Responses to TPU Shock", size(medium)) ///
    note("68% Bootstrap CI, Horizon: 36 months", size(small)) ///
    xlabel(0(6)36, labsize(small)) ylabel(,grid labsize(small)) ///
    graphregion(color(white)) bgcolor(white) ///
    byopts(rows(3) yrescale xrescale compact) ///
    name(combined_tpu_pre, replace) ysize(35) xsize(40)
graph export "irf_tpu_pre_brics.pdf", replace as(pdf)


* Forecast Error Variance Decomposition (%)
irf table fevd, impulse(tpu) response(tpu d_ipi d_reer d_cpi d_tradebalance d_stockmarket)


* Estimate VAR with TPU and BRICS interaction
gen tpu_brics = tpu * brics_dummy
varsoc tpu tpu_brics d_ipi d_reer d_cpi d_tradebalance d_stockmarket, maxlag(12)

var tpu tpu_brics d_ipi d_reer d_cpi d_tradebalance d_stockmarket, lags(1/4)
estimates store var_tpu_brics

* Generate IRFs for BRICS interaction effect
irf create tpu_post, set(myirf3b, replace) step(36) bs reps(50)
irf set myirf3b
irf graph oirf, impulse(tpu_brics) response(tpu d_ipi d_reer d_cpi d_tradebalance d_stockmarket) ///
    ci level(68) lstep(0) ustep(36) ///
    title("Ad. Effect of TPU Shock Post-BRICS Period", size(medium)) ///
    note("68% Bootstrap CI, Horizon: 36 months", size(small)) ///
    xlabel(0(6)36, labsize(small)) ylabel(,grid labsize(small)) ///
    graphregion(color(white)) bgcolor(white) ///
    byopts(rows(3) yrescale xrescale compact) ///
    name(combined_tpu_brics, replace) ysize(35) xsize(40)
graph export "irf_tpu_brics_effect.pdf", replace as(pdf)
