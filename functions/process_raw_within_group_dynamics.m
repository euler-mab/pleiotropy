function [] = process_raw_within_group_dynamics(verbose, two_cell_experiment, cooperation_experiment, y_end, rho, K, s, phi, mu, nu, num_traits, num_types, num_group_types, Z)

desired_replicates = 10000;
filename_within_output = sprintf('processed_cell_results/two%d_coop%d_rep%d_tend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', two_cell_experiment, cooperation_experiment, desired_replicates, y_end, rho, K, s, phi, mu, nu);
filename_within_output = strrep(filename_within_output, '.', '-');
filename_within_output = strcat(filename_within_output, '.mat');

if isfile(filename_within_output)
    
    % Do nothing.
    
else
    
    % Load the abundances - verbose.
    [n_output, fe_b_output, fe_count_output, fe_d_output, replicates] = load_raw_within_group_dynamics(verbose, two_cell_experiment, cooperation_experiment, y_end, rho, K, s, phi, mu, nu);
    
    % Average over replicates.
    n_output = squeeze(mean(n_output, 3));
    
    % Process.
    if two_cell_experiment == 0
                       
        n_output = permute(n_output, [1 3 2]);
        
        for i = 1:num_group_types
            
            for y = 1:(y_end + 1)
                
                n_y(y, i, :) = n_output(y, i, :);
                N_y(y, i) = squeeze(sum(n_output(y, i, :), 3));
                x_y(y, i, :) = n_output(y, i, :) / N_y(y, i);
                z_bar_y(y, i, :) = Z' * squeeze(x_y(y, i, :));
                
%                 for k = 1:num_types
%                     
%                     z_social_y(y, i, k, :) = Z(k, :);
%                     z_prime_social_y(y, i, k, :) = z_bar_y(y, i, :); % In continuous models social partner is just the average in the group
%                     
%                 end
                
            end
            
        end
        
    else        
        
        idx_group_types = tril(reshape(1:num_types^2, [num_types, num_types]), 0);
        idx_group_types = idx_group_types(:);
        idx_group_types = idx_group_types(idx_group_types > 0);
        
        [idx_a, idx_b] = ind2sub([num_types, num_types], idx_group_types);
        
        n_output = permute(n_output, [1 3 4 2]);
        
        for i = 1:num_group_types
            
            for y = 1:(y_end + 1)
                
                n_y(y, i, :) = squeeze(n_output(y, idx_a(i), idx_b(i), :));
                N_y(y, i) = squeeze(sum(n_y(y, i, :), 3));
                x_y(y, i, :) = n_y(y, i, :) / N_y(y, i);
                z_bar_y(y, i, :) = Z' * squeeze(x_y(y, i, :));
                
                % Only used for social analysis
                
%                 for k = 1:num_types
%                     
%                     z_social_y(y, i, k, :) = Z(k, :);
%                     z_prime_social_y(y, i, k, :) = z_bar_y(y, i, :); % in continuous models social partner is just the average in the group
%                     
%                 end
                
            end
            
        end
        
    end
    
    % Save the processed data.
    filename_within_output = sprintf('processed_cell_results/two%d_coop%d_rep%d_tend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', two_cell_experiment, cooperation_experiment, replicates, y_end, rho, K, s, phi, mu, nu);
    filename_within_output = strrep(filename_within_output, '.', '-');
    filename_within_output = strcat(filename_within_output, '.mat');
    
    save(filename_within_output, 'n_y', 'N_y', 'x_y', 'z_bar_y', 'fe_count_output', 'fe_b_output', 'fe_d_output', '-v7.3');
    
end


end