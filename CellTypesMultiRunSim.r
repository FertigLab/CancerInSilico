dist_inc <- 5 # percent increment of cell type A distribution

sim_runs_list = list() # list to hold CellModel output from each sim

for (a_dist in seq(0, 100, dist_inc)) {
	
	percent_a_dist = a_dist/100
	single_cell_run <- runCancerSim(initialNum = 1, runTime = 84, cellTypeInitFreq = c(percent_a_dist, 1 - percent_a_dist))
	sim_runs_list[[length(sim_runs_list) + 1]] <- single_cell_run
	
}

# analyze data