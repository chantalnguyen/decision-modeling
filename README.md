mle_Ind50: determine best-fit parameters of power-law decision model for Ind50 trials with maximum likelihood estimation, and solve master equation for each trial

	requires:
	llfun: calculate log-likelihood function
	mastereq: solve master equation
	
mle_Ind50_evac: does the same as mle_Ind50 but fits to evacuations instead of decisions

cv_Ind50: perform leave-one-out cross-validation for Ind50 trials

plot_trials: plot data, model prediction, cross-validation prediction, and P_hit trajectories

	requires:
	shadedErrorBar: creates shaded regions representing confidence intervals
	tightfig: reduces whitespace in saved figure
	
group_scores: calculates and plots scores in group games
	requires:
	calc_scores: calculates scores from raw evacuation data
