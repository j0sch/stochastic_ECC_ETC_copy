# stochastic_ECC_ETC

To generate simulations of 3 dyads stimulated 3 times:

Run Coupl_rep_uncoupl_Mult_clusters_multiple_runs_save_commented.m
	Dependencies:
		- parsave.m
		- dec2str.m
		- RK4_Sobie_Cao_2_diffusion.m
			- Sobie_Cao_fluxes_combined_2.m
		- sobie_etal_2002_combined_parameters_3.m
		- Cao_et_al_2014.m
		- update_ip3r_states.m
	Running on Huia, this takes approx. 8hrs per simulation. (On 25 cores, 8hrs/25 sims)

Run plot_rep_stimulated_uncoupled_data.m to condense simulations into one file.
	Dependencies:
		- make sure it can find sim files



To run simulations of 3 dyads unstimulated:
Run Coupl_rep_uncoupl_Mult_clusters_multiple_runs_save_0_st.m
	Dependencies:
		- parsave.m
		- dec2str.m
		- RK4_Sobie_Cao_2_diffusion.m
			- Sobie_Cao_fluxes_combined_2.m
		- sobie_etal_2002_combined_parameters_3.m
		- Cao_et_al_2014.m
		- update_ip3r_states.m
	Running on Huia, this takes approx. 8hrs per simulation. (On 25 cores, 8hrs/25 sims)

Run plot_rep_stimulated_uncoupled_data_no_st.m to condense simulations into one file.
	Dependencies:
		- make sure it can find sim files


To run 1 dyad unstimulated:
Run one_dyad_Mult_clusters_multiple_runs_save_tst_2_only_ryr.m (one dyad with only RyRs) and one_dyad_Mult_clusters_multiple_runs_save_tst_2.m (one dyad with both RyRs and IP3Rs)
	Dependencies:
		- parsave.m
		- dec2str.m
		- RK4_Sobie_Cao_2_diffusion.m
			- Sobie_Cao_fluxes_combined_2.m
		- sobie_etal_2002_combined_parameters_3.m
		- Cat_et_al_2014.m
		- update_ip3r_states.m
	Running on Huia, this takes approx. 5hrs per simulation. (On 25 cores, 5hrs/25 sims)

Run plot_rep_stimulated_uncoupled_data_one_d.m to condense simulations into one file.
	Dependencies:
		- make sure it can find sim files

To generate plots for paper:
Run the first section of plot_numbers_distances_and_all_that_jazz.m. The code for the plots used in Llew's paper should now work. It is located at the bottom of the file. The rest of the code is for other plots.
