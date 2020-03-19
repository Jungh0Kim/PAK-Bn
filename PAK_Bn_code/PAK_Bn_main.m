% Jungho Kim, Junho Song
% Probabilisy-Adaptive Kriging in n-ball (PAK-Bn) for Reliability Analysis,
% Structural Safety (2020)
% https://github.com/Jungh0Kim/PAK-Bn

clear; close all; tic;

%% Set Initial Parameters

n_dim = 2;             % dimension of random variables
N_sample = 5e4;        % number of samples
N_iniDOE = 10;         % initial number of DoE
max_iter = 100;
max_round = 10;

LHS_bound = 6.0;       % LHS bound
R1 = 5.0;              % initial n-ball radius
del_R = 0.20;
Beta_tol1 = 1e-4;      % tolerance1
Beta_tol2 = 1e-5;      % tolerance2
P_f_MCS = 2.21e-3;     % reference for numerical example

% Initial DoE
rng(13) % For reproducibility
x_doe = DoE_byLHS(n_dim,N_iniDOE,LHS_bound);
G_doe = G_function(x_doe(:,1), x_doe(:,2));

%% Specify the parameters of gp

meanfunc = [];
covfunc = {@covMaterniso, 5}; hyp.cov = [1/2 1/2];   % Matern covariance function
likfunc = @likGauss; sn = 1e-3; hyp.lik = log(sn);   % Gaussian likelihood
% hyp_MLE = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x_training, G_training);

%% Active-Learning

Preallocation;
R_max = R1+del_R*max_round;
p = haltonset(n_dim,'Skip',1e3,'Leap',1e2); % quasi-random samples
p = scramble(p,'RR2');
total_sample = -R_max + 2*R_max*p(1:N_sample,:); % total sample (memory)

k = 1; t_j = 0;
Ad_radi = [0 R1:del_R:R_max];
while (k <= max_round)

    % n-ball sampling
    I_t = (nB_domain(total_sample, Ad_radi(k)) >= 0 & nB_domain(total_sample, Ad_radi(k+1)) < 0);
    [sp_pt,~] = find(I_t==1);
    Bn_samp = total_sample(sp_pt,:);
    Bn_samp_chg = Bn_samp;
    N_samp_R = size(Bn_samp,1);
    
    Vol_m = pi^(n_dim/2)/gamma(n_dim/2+1)*Ad_radi(k+1)^n_dim; % Volume of n-ball
    index_s = cumprod(normpdf(Bn_samp_chg,0,1),2);
    index = index_s(:,end);
    N_samp_R_cell{k,1} = N_samp_R;
    if k>1
        N_samp_R_cell{k,1} = N_samp_R + N_samp_R_cell{k-1,1};
    end
    clear P_f CP_hat
    
    %% k = 1    
    if k==1 % 1st iteration
        iter = 1;
        % divide eval and non-eval points
        index_eval = []; index_neval =  index(:);
        
        while (iter <= max_iter)
            % Kriging estimation        
            [mu, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x_doe, G_doe, Bn_samp_chg);
                   
            if iter==1
                [CP_hat_s,find_1_id] = Find_CP1(mu,index_neval,Bn_samp);
                CP_hat(iter,:) = CP_hat_s;
                P_f(iter,1) = sum(index(find_1_id)) * Vol_m / N_samp_R;

                % Penalty function
                if ~isempty(find_1_id)
                    penalty_factor = abs(Sqrt_ss(Bn_samp_chg) - norm(CP_hat(iter,:),2)) / Ad_radi(k+1);
                else
                    penalty_factor = 1;
                end
            else % iter > 1
                [CP_hat_s,find_1_id,find_2_id] = Find_CP2(mu,x_doe,G_doe,N_iniDOE,index_neval,index_eval,Bn_samp_chg);
                CP_hat(iter,:) = CP_hat_s;
                P_f(iter,1) = (sum(index_neval(find_1_id)) + sum(index_eval(find_2_id))) * Vol_m / N_samp_R;
                
                % Penalty function
                if (~isempty(find_1_id) || ~isempty(find_2_id))
                    penalty_factor = abs(Sqrt_ss(Bn_samp_chg) - norm(CP_hat(iter,:),2)) / Ad_radi(k+1);
                    if abs( norm(CP_hat(iter,:)-CP_hat(iter-1,:),2) ) < 0.01
                        penalty_factor = 1;
                    end 
                else % not find CP
                    penalty_factor = 1;                
                end
                
                % Convergence condition
                Beta_g_1 = -norminv(P_f(iter-1,1));
                Beta_g_2 = -norminv(P_f(iter,1));
                Beta_var = abs(Beta_g_2 - Beta_g_1) / Beta_g_1;
                if Beta_var < Beta_tol1 || isnan(Beta_var)==1
                    break
                end  
            end  % if end
            
            % Minimize the learning function
            [LF_min,i_min] = Learning_func(mu, sqrt(s2),penalty_factor);  
            
            % Add the new evaluation point to the DoE
            x_doe = vertcat(x_doe, Bn_samp_chg(i_min,:));
            G_doe = vertcat(G_doe, G_function(Bn_samp_chg(i_min,1),Bn_samp_chg(i_min,2)));  % +1 G evaluation    
            index_eval = vertcat(index_eval, index_neval(i_min));  
            Bn_samp_chg = [Bn_samp_chg(1:i_min-1,:); Bn_samp_chg(i_min+1:end,:)];
            index_neval = [index_neval(1:i_min-1); index_neval(i_min+1:end)];
            iter = iter + 1;
        end % while (iter) end
        
        x_cell{k,1} = x_doe;
        G_cell{k,1} = G_doe;
        Bn_chg_cell{k,1} = Bn_samp_chg;
        CP_eval_cell{k,1} = index_eval;
        CP_neval_cell{k,1} = index_neval;

        % Estimate failure probability of the last model
        [mu, ~] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x_doe, G_doe, Bn_samp_chg);
        [CP_hat_s,find_1_id,find_2_id] = Find_CP2(mu,x_doe,G_doe,N_iniDOE,index_neval,index_eval,Bn_samp_chg);
        CP_hat(iter,:) = CP_hat_s;
        CP_hat_cell{k,1} = CP_hat;
        
        P_f(iter,1) = (sum(index_neval(find_1_id)) + sum(index_eval(find_2_id))) * Vol_m / N_samp_R;
        P_f_cell{k,1} = P_f;
        iter = iter - 1; 
        iter_cell{k,1} = iter;
        
    %% k > 1    
    else  % from 2nd round to end (k > 1)
        iter = 1;
        x_doe = x_cell{k-1,1};
        G_doe = G_cell{k-1,1};
        Bn_samp_chg = vertcat(Bn_chg_cell{k-1,1},Bn_samp_chg);
        index_eval = CP_eval_cell{k-1,1};
        index_neval = vertcat(CP_neval_cell{k-1,1}, index);
        
        while (iter <= max_iter)
            % Kriging estimation
            [mu, s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x_doe, G_doe, Bn_samp_chg);
            [CP_hat_s,find_1_id,find_2_id] = Find_CP2(mu,x_doe,G_doe,N_iniDOE,index_neval,index_eval,Bn_samp_chg);
            CP_hat(iter,:) = CP_hat_s;
        
            P_f(iter,1) = (sum(index_neval(find_1_id)) + sum(index_eval(find_2_id))) * Vol_m / N_samp_R_cell{k,1};
            
            % Penalty function
            if (~isempty(find_1_id) || ~isempty(find_2_id))
                penalty_factor = abs(Sqrt_ss(Bn_samp_chg) - norm(CP_hat(iter,:),2)) / Ad_radi(k+1); 
                if iter > 1
                    if abs( norm(CP_hat(iter,:)-CP_hat(iter-1,:),2) ) < 0.01
                        penalty_factor = 1; 
                    end
                end
                if iter == 1 && isempty(CP_hat_cell{k-1,1})==0
                    if abs( norm(CP_hat(iter,:)-CP_hat_cell{k-1,1}(end,:),2) ) < 0.01
                        penalty_factor = 1;  
                    end
                end   
            else % not find CP
                penalty_factor = 1;                
            end 

            Beta_g_j1 = -norminv(P_f_cell{k-1,1}(end));
            Beta_g_j2 = -norminv(P_f(1,1));
            Beta_var_j = abs(Beta_g_j2 - Beta_g_j1) / Beta_g_j1;
            if  Beta_var_j < Beta_tol2
                t_j = 1;
                break % break iteration (whole while loop)
            end
            
            if iter > 1
                Beta_g_1 = -norminv(P_f(iter-1,1));
                Beta_g_2 = -norminv(P_f(iter,1));
                Beta_var = abs(Beta_g_2 - Beta_g_1) / Beta_g_1;
                
                if Beta_var < Beta_tol1 || isnan(Beta_var)==1
                    break % break iteration (while iter loop)
                end
            end  % if end

            % Minimize the learning function
            [LF_min,i_min] = Learning_func(mu, sqrt(s2),penalty_factor);

            % Add the new evaluation point to the DoE
            x_doe = vertcat(x_doe, Bn_samp_chg(i_min,:));
            G_doe = vertcat(G_doe, G_function(Bn_samp_chg(i_min,1),Bn_samp_chg(i_min,2))); % +1 G evaluation 
            index_eval = vertcat(index_eval, index_neval(i_min));      
            Bn_samp_chg = [Bn_samp_chg(1:i_min-1,:); Bn_samp_chg(i_min+1:end,:)];
            index_neval = [index_neval(1:i_min-1); index_neval(i_min+1:end)];
            iter = iter + 1;
        end  % while (iter) end
        
        if  t_j==1
            break % stop whole process
        end
        x_cell{k,1} = x_doe;
        G_cell{k,1} = G_doe;
        Bn_chg_cell{k,1} = Bn_samp_chg;
        CP_eval_cell{k,1} = index_eval;
        CP_neval_cell{k,1} = index_neval;

        % Estimate failure probability of the last model
        [mu, ~] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, x_doe, G_doe, Bn_samp_chg);
        [CP_hat_s,find_1_id,find_2_id] = Find_CP2(mu,x_doe,G_doe,N_iniDOE,index_neval,index_eval,Bn_samp_chg);
        CP_hat(iter,:) = CP_hat_s;
        CP_hat_cell{k,1} = CP_hat;

        P_f(iter,1) = (sum(index_neval(find_1_id)) + sum(index_eval(find_2_id))) * Vol_m / N_samp_R_cell{k,1};
        P_f_cell{k,1} = P_f;
        iter = iter - 1; 
        iter_cell{k,1} = iter;
    end  % Active-learning end
    k = k + 1; % expanding the n-ball
end  % k end 
k = k - 1; toc;

%% Results

Plot_Design;
Plot_Converg;

% Display the added DoE
fprintf ('added DOE : x =\n');
disp(x_doe(N_iniDOE+1:end,:)');

% Total number of evaluations
fprintf('Number of evaluations: %d + %d = %d.\n', N_iniDOE, sum([iter_cell{:}]), N_iniDOE + sum([iter_cell{:}]));
