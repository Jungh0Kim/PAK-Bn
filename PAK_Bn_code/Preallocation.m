% Pre-allocation for data

N_samp_R_cell = cell(max_round,1);
CP_true = zeros(max_round,n_dim);
iter_cell = cell(max_round,1);
x_cell = cell(max_round,1);
G_cell = cell(max_round,1);
Bn_chg_cell = cell(max_round,1);
CP_eval_cell = cell(max_round,1);
CP_neval_cell = cell(max_round,1);
P_f_cell = cell(max_round,1);
CP_hat_cell = cell(max_round,1);