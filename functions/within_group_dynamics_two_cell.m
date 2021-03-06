function [n_output, fe_count_output, fe_b_output, fe_d_output] = within_group_dynamics_two_cell(verbose, cooperation_experiment, replicates, t_0, t_by, t_end, rho, K, s_c, phi, mu, nu, num_traits, num_types, Z)

tic;

%% Parameters

two_cell_experiment = 1;

%% Within-group calculated parameters

numu = nu * mu;

Q = [
    [(1 - numu) * (1 - numu) * (1 - numu), (1 - numu) * (1 - numu) * numu, (1 - numu) * numu * (1 - numu), (1 - numu) * numu * numu, numu * (1 - numu) * (1 - numu), numu * (1 - numu) * numu, numu * numu * (1 - numu), numu * numu * numu];
    [(1 - numu) * (1 - numu) * mu, (1 - numu) * (1 - numu) * (1 - mu), (1 - numu) * numu * mu, (1 - numu) * numu * (1 - mu), numu * (1 - numu) * mu, numu * (1 - numu) * (1 - mu), numu * numu * mu, numu * numu * (1 - mu)];
    [(1 - numu) * mu * (1 - numu), (1 - numu) * mu * numu, (1 - numu) * (1 - mu) * (1 - numu), (1 - numu) * (1 - mu) * numu, numu * mu * (1 - numu), numu * mu * numu, numu * (1 - mu) * (1 - numu), numu * (1 - mu) * numu];
    [(1 - numu) * mu * mu, (1 - numu) * mu * (1 - mu), (1 - numu) * (1 - mu) * mu, (1 - numu) * (1 - mu) * (1 - mu), numu * mu * mu, numu * mu * (1 - mu), numu * (1 - mu) * mu, numu * (1 - mu) * (1 - mu)];
    [mu * (1 - numu) * (1 - numu), mu * (1 - numu) * numu, mu * numu * (1 - numu), mu * numu * numu, (1 - mu) * (1 - numu) * (1 - numu), (1 - mu) * (1 - numu) * numu, (1 - mu) * numu * (1 - numu), (1 - mu) * numu * numu];
    [mu * (1 - numu) * mu, mu * (1 - numu) * (1 - mu), mu * numu * mu, mu * numu * (1 - mu), (1 - mu) * (1 - numu) * mu, (1 - mu) * (1 - numu) * (1 - mu), (1 - mu) * numu * mu, (1 - mu) * numu * (1 - mu)];
    [mu * mu * (1 - numu), mu * mu * numu, mu * (1 - mu) * (1 - numu), mu * (1 - mu) * numu, (1 - mu) * mu * (1 - numu), (1 - mu) * mu * numu, (1 - mu) * (1 - mu) * (1 - numu), (1 - mu) * (1 - mu) * numu];
    [mu * mu * mu, mu * mu * (1 - mu), mu * (1 - mu) * mu, mu * (1 - mu) * (1 - mu), (1 - mu) * mu * mu, (1 - mu) * mu  * (1 - mu), (1 - mu) * (1 - mu) * mu, (1 - mu) * (1 - mu) * (1 - mu)];
    ];

P = [
    [(1 - numu) * (1 - numu) * (1 - numu), (1 - numu) * (1 - numu) * numu, (1 - numu) * numu * (1 - numu), (1 - numu) * numu * numu, numu * (1 - numu) * (1 - numu), numu * (1 - numu) * numu, numu * numu * (1 - numu), numu * numu * numu];
    [(1 - numu) * (1 - numu) * mu, (1 - numu) * (1 - numu) * (1 - mu), (1 - numu) * numu * mu, (1 - numu) * numu * (1 - mu), numu * (1 - numu) * mu, numu * (1 - numu) * (1 - mu), numu * numu * mu, numu * numu * (1 - mu)];
    [(1 - numu) * mu * (1 - numu), (1 - numu) * mu * numu, (1 - numu) * (1 - mu) * (1 - numu), (1 - numu) * (1 - mu) * numu, numu * mu * (1 - numu), numu * mu * numu, numu * (1 - mu) * (1 - numu), numu * (1 - mu) * numu];
    [(1 - numu) * mu * mu, (1 - numu) * mu * (1 - mu), (1 - numu) * (1 - mu) * mu, (1 - numu) * (1 - mu) * (1 - mu), numu * mu * mu, numu * mu * (1 - mu), numu * (1 - mu) * mu, numu * (1 - mu) * (1 - mu)];
    [mu * (1 - numu) * (1 - numu), mu * (1 - numu) * numu, mu * numu * (1 - numu), mu * numu * numu, (1 - mu) * (1 - numu) * (1 - numu), (1 - mu) * (1 - numu) * numu, (1 - mu) * numu * (1 - numu), (1 - mu) * numu * numu];
    [mu * (1 - numu) * mu, mu * (1 - numu) * (1 - mu), mu * numu * mu, mu * numu * (1 - mu), (1 - mu) * (1 - numu) * mu, (1 - mu) * (1 - numu) * (1 - mu), (1 - mu) * numu * mu, (1 - mu) * numu * (1 - mu)];
    [mu * mu * (1 - numu), mu * mu * numu, mu * (1 - mu) * (1 - numu), mu * (1 - mu) * numu, (1 - mu) * mu * (1 - numu), (1 - mu) * mu * numu, (1 - mu) * (1 - mu) * (1 - numu), (1 - mu) * (1 - mu) * numu];
    [mu + mu - mu * mu, 0, 0, 0, 0, 0, (1 - mu) * (1 - mu) * mu, 1 - (mu + mu - mu * mu) - (1 - mu) * (1 - mu) * mu];
    ];

H = (1 - phi) * Q + phi * P;

%% Output

id = string(java.util.UUID.randomUUID.toString);

filename_within_output = sprintf('cell_results/two%d_coop%d_rep%d_tend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', two_cell_experiment, cooperation_experiment, replicates, t_end, rho, K, s_c, phi, mu, nu);
filename_within_output = strrep(filename_within_output, '.', '-');
filename_within_output = strcat(filename_within_output, '_id_', id);
filename_within_output = strcat(filename_within_output, '.mat');

n_output = zeros(t_end / t_by, num_types, replicates, num_types, num_types);

fe_count_output = zeros(t_end / t_by, num_types, num_types, replicates, num_types, num_types);
fe_b_output = zeros(t_end / t_by, num_types, num_types, replicates, num_types, num_types);
fe_d_output = zeros(t_end / t_by, num_types, num_types, replicates, num_types, num_types);

%% Simulate

for j = 1:num_types
    
    for i = 1:num_types
        
        focal_type_a = i;
        focal_type_b = j;
        
        if i >= j
            
            replicate = 1;
            
            while replicate <= replicates
                
                disp(strcat('ID:', {' '}, id, {' '}, ' Type: ', {' ('}, num2str(focal_type_a), {', '}, num2str(focal_type_b), {') '}, 'Replicate: ', {' '}, num2str(replicate)))
                
                %% Specify initial conditions
                
                t_step = 1;
                
                n_t0 = zeros(1, num_types)';
                n_t0(focal_type_a) = n_t0(focal_type_a) + 1;
                n_t0(focal_type_b) = n_t0(focal_type_b) + 1;
                
                Z_t0 = Z;
                H_t0 = H;
                
                %% Set initial conditions
                
                t = t_0;
                n_t = n_t0;
                Z_t = Z_t0;
                H_t = H_t0;
                
                %% Initial output
                
                n_output(1, :, replicate, focal_type_a, focal_type_b) = n_t0;
                
                while t <= t_end
                    
                    N_t = sum(n_t);
                    x_t = n_t / N_t;
                    Z_t = Z_t0;
                    z_bar_t = Z_t' * x_t;
                    
                    b_t = rho .* ((1 - s_c + s_c * (Z_t(:, 1))) ./ (1 - s_c + s_c * (z_bar_t(1)))) * 1;
                    if cooperation_experiment
                        b_t = rho .* ((1 - s_c + s_c * (1 - Z_t(:, 1))) ./ (1 - s_c + s_c * (1 - z_bar_t(1)))) * 1;
                    end
                    d_t = rho .* ((1 - s_c + s_c * (1 - Z_t(:, 2))) ./ (1 - s_c + s_c * (1 - z_bar_t(2)))) * (N_t - 1) / K;
                    p_t = H_t' * diag(n_t) * b_t;
                    
                    b_bar_t = b_t' * x_t;
                    d_bar_t = d_t' * x_t;
                    p_bar_t = sum(p_t) / N_t;
                    if N_t == 0
                        p_bar_t = 0;
                    end
                    
                    B_t = b_t' * n_t;
                    D_t = d_t' * n_t;
                    P_t = sum(p_t);
                    E_t = P_t + D_t;
                    
                    b_p = (b_t .* n_t) / B_t;
                    if B_t == 0
                        b_p = zeros(1, num_types)';
                    end
                    
                    d_p = (d_t .* n_t) / D_t;
                    if D_t == 0
                        d_p = zeros(1, num_types)';
                    end
                    
                    p_p = p_t / P_t;
                    if P_t == 0
                        p_p = zeros(1, num_types)';
                    end
                    
                    %     sum(p_p) == sum(b_p)
                    
                    random_event = rand() * E_t;
                    
                    v_t = zeros(1, num_types)';
                    if random_event <= P_t % birth event
                        idx_birth = randsample(1:num_types, 1, true, b_p);
                        idx_produced = randsample(1:num_types, 1, true, H_t(idx_birth, :));
                        v_t(idx_produced) = 1;
                        
                        if verbose == 1
                            
                            % fitness effect
                            fe_n_t = n_t + v_t;
                            fe_N_t = sum(fe_n_t);
                            fe_x_t = fe_n_t / fe_N_t;
                            fe_Z_t = Z_t0;
                            fe_z_bar_t = fe_Z_t' * fe_x_t;
                            
                            fe_b_t = rho .* ((1 - s_c + s_c * (fe_Z_t(:, 1))) ./ (1 - s_c + s_c * (fe_z_bar_t(1)))) * 1;
                            if cooperation_experiment
                                fe_b_t = rho .* ((1 - s_c + s_c * (1 - fe_Z_t(:, 1))) ./ (1 - s_c + s_c * (1 - fe_z_bar_t(1)))) * 1;
                            end
                            fe_d_t = rho .* ((1 - s_c + s_c * (1 - fe_Z_t(:, 2))) ./ (1 - s_c + s_c * (1 - fe_z_bar_t(2)))) * (fe_N_t - 1) / K;
                            
                            fe_count_output(t_step, focal_type_a, focal_type_b, replicate, idx_birth, idx_produced) = fe_count_output(t_step, focal_type_a, focal_type_b, replicate, idx_birth, idx_produced) + 1;
                            fe_b_output(t_step, focal_type_a, focal_type_b, replicate, idx_birth, idx_produced) = fe_b_output(t_step, focal_type_a, focal_type_b, replicate, idx_birth, idx_produced) + fe_b_t(idx_produced) / fe_b_t(idx_birth);
                            fe_d_output(t_step, focal_type_a, focal_type_b, replicate, idx_birth, idx_produced) = fe_d_output(t_step, focal_type_a, focal_type_b, replicate, idx_birth, idx_produced) + fe_d_t(idx_produced) / fe_d_t(idx_birth);
                            
                        end
                        
                    else % death event
                        idx_death = randsample(1:num_types, 1, true, d_p);
                        v_t(idx_death) = -1;
                    end
                    
                    n_t = n_t + v_t;
                    t = t - log(rand()) / E_t;
                    
                    %% Record output
                    
                    if t >= t_step
                        t_step = t_step + t_by;
                        n_output(t_step, :, replicate, focal_type_a, focal_type_b) = n_t;
                    end
                    
                end
                
                replicate = replicate + 1;
                
            end
            
        end
        
    end
    
end

save(filename_within_output, 'n_output', 'fe_count_output', 'fe_b_output', 'fe_d_output', '-v7.3')

toc

end

