
function [ n_group_t_y, z_bar_group_t, n_t, z_bar_t, n_y, z_bar_y ] = between_group_dynamics(two_cell_experiment, lambda, K_g, s_g, c_g, alpha, cooperation_experiment, t_0, t_by, t_end, rho, K, s_c, phi, mu, nu, num_traits, num_types, num_group_traits, num_group_types, Z)

%% Load within group files

listing = dir('processed_cell_results/');

file_list = {listing.name};
file_list = file_list(3:end);

search1 = sprintf('two%d_coop%d', two_cell_experiment, cooperation_experiment);
search1 = strrep(search1, '.', '-');

search2 = sprintf('tend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', t_end, rho, K, s_c, phi, mu, nu);
search2 = strrep(search2, '.', '-');

idx_files1 = ~cellfun(@isempty, strfind(file_list, search1));
idx_files2 = ~cellfun(@isempty, strfind(file_list, search2));

idx_files = logical(idx_files1 .* idx_files2);

file_list = file_list(idx_files);

file = load(['processed_cell_results/', file_list{1}], 'n_y', 'N_y', 'x_y', 'z_bar_y');
n_y = file.n_y;
N_y = file.N_y;
x_y = file.x_y;
z_bar_y = file.z_bar_y;

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

%% B-group calculated parameters

rho = 1 / lambda;

if two_cell_experiment   
    idx_group_types = tril(reshape(1:num_types^2, [num_types, num_types]), 0);
    idx_group_types = idx_group_types(:);
    idx_group_types = idx_group_types(idx_group_types > 0);
end

%% Setting up grid

% Centres

% Age dimensions
y_0 = t_0; % begginging of space
y_end = t_end; % end of space
dy = t_by; % change in space

% Time dimensions
t_0 = 0; % begginging of time
t_end = lambda * 2500; % 2500 generations maximum end.
dt = 0; % change in time (to be determined later by CFL condition)
C = 0.95; % Courant number - proportion of maximum allowable time step

% Number of cells
num_centres = floor(y_end / dy + 1);

% Number of edges
num_edges = num_centres + 1;

% Number of ghost cells
num_ghosts = 2;

% Composition space cell centres
y_centres = (y_0:dy:y_end)'; % composition space cells

% Add 2 ghost cells
tmp = zeros(num_ghosts + num_centres + num_ghosts, 1);
tmp(num_ghosts + 1:num_centres + num_ghosts) = y_centres;
y_centres = tmp;
clear tmp;

% Indices of cell centres
idx_centres = (num_ghosts + 1:num_centres + num_ghosts)';

% Edges

% Composition space cell edges
y_edges = ((y_0 - dy / 2):dy:(y_end + dy / 2))'; % composition space cells

% Add 2 ghost cells
tmp = zeros(num_ghosts + num_edges + num_ghosts, 1);
tmp(num_ghosts + 1:num_edges + num_ghosts) = y_edges;
y_edges = tmp;
clear tmp;

% Indices of cell edges
idx_edges = (num_ghosts + 1:num_edges + num_ghosts)';

% Domains of interest.

% Indices of initial group density. 
idx_initial = num_ghosts + find(y_centres(idx_centres) >= 0 & y_centres(idx_centres) <= lambda);

% Indices of descendants.
idx_descendants = num_ghosts + find(y_centres(idx_centres) >= 0 & y_centres(idx_centres) < 1);
num_descendantCells = length(idx_descendants);

% Indices of reproductives.
idx_reproductives = num_ghosts + find(y_centres(idx_centres) >= alpha * lambda);

%%

% Pre-allocate matrices.

% Group density at cell centres.
n_group_t_y = zeros(num_ghosts + num_centres + num_ghosts, num_group_types);

% Group birth at cell centres.
b_group_t_y = zeros(num_ghosts + num_centres + num_ghosts, num_group_types);

% Group production at cell centres.
p_group_t_y = zeros(num_ghosts + num_centres + num_ghosts, num_group_types);

% Group death at cell centres.
d_group_t_y = zeros(num_ghosts + num_centres + num_ghosts, num_group_types);

% Group inheritance matrix at cell centres.
H_group_t_y = zeros(num_ghosts + num_centres + num_ghosts, num_group_types, num_group_types);

% Group trait values at cell centres.
z_group_t_y = zeros(num_ghosts + num_centres + num_ghosts, num_group_types, num_group_traits);

% Velocity at cell edges
dy_dt_edges = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);

% Inidicator at cell edges
flip_flop_edges = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);
indicator_edges = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);

% Slope function at cell centres
r_edges = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);

% Slope function at cell centres
phi_edges_donor_cell = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);

% Slope function at cell centres
phi_edges_superbee = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);

% Flux along cell edges
flux_edges = zeros(num_ghosts + num_edges + num_ghosts, num_group_types);

%% Get within-group dynamics and set group inheritance and trait values.

        
% Calculate group trait values.
for i = 1:num_group_types
    for y = 1:num_centres            
        z_group_t_y(y + num_ghosts, i, 1) = z_bar_y(y, i, 1);
        z_group_t_y(y + num_ghosts, i, 2) = z_bar_y(y, i, 2);
        z_group_t_y(y + num_ghosts, i, 3) = z_bar_y(y, i, 3);
        z_group_t_y(y + num_ghosts, i, 4) = z_bar_y(y, i, 1) .* z_bar_y(y, i, 2) .* (1 - c_g .* z_bar_y(y, i, 3));            
    end
end

z_group_t_y(idx_reproductives, :, 5) = 1;
    
% Calculate group trait values.
for i = 1:num_group_types
    for y = 1:num_centres            
        if two_cell_experiment == 0            
            H_group_t_y(y + num_ghosts, i, :) = squeeze(x_y(y, i, :))' * H;                           
        else                        
            temp = (squeeze(x_y(y, i, :))' * H)' * (squeeze(x_y(y, i, :))' * H);
            temp = tril(temp + triu(temp, 1)');
            temp = temp(:);
            temp = double(temp(idx_group_types));            
            H_group_t_y(y + num_ghosts, i, :) = temp;
        end
    end
end

%% Initial conditions

% Initial time
t = t_0;

% Initial group density
n_group_t_y(idx_initial, 1) = K_g;
n_group_t_y = n_group_t_y / sum(n_group_t_y, [1, 2]);

% Initial total group density
n_group_t = sum(n_group_t_y(idx_centres, :)) * dy;

% Initial velocity
dy_dt_edges(idx_edges, :) = 1;
dy_dt_edges(3,:) = 0;
% dn1dtEdges(numGhostCells + numEdges) = dn1dtEdges(numGhostCells + 1); % For periodic boundary conditions, the final edge must be the same as the first edge

% Initial flux limiter function at cell edges
flip_flop_edges(dy_dt_edges >= 0, :) = 1;
flip_flop_edges(dy_dt_edges < 0, :) = -1;
indicator_edges(dy_dt_edges >= 0, :) = 1;
indicator_edges(dy_dt_edges < 0, :) = 0;

%% Output

% Save the group output
n_group_output = zeros(t_end, num_ghosts + num_centres + num_ghosts, num_group_types);
z_bar_group_output = zeros(t_end, num_group_traits);

%% Specify initial conditions

t_step = 1;

% Time step - dependent on CFL condition and Courant number
dt = C * dy / max(dy_dt_edges, [], 'all'); % From Dullemond (2013; Chapter 3; pg. 46)

%% Initial output

n_group_output(t_step, :, :) = n_group_t_y;

% z_bar_group_output(t_step, :) = n_group_t;

while t <= t_end
         
    %% Update population statistics.
        
    % Calculate the total size of the group population.
    N_group_t = sum(sum(n_group_t_y(idx_centres, :)) * dy);
    
    % Calculate average group phenotype.
    z_bar_group_t = squeeze(sum(z_group_t_y(idx_centres, :, :) .* n_group_t_y(idx_centres, :) * dy, [1, 2]) / N_group_t);
              
    %% Update group birth and death rates
        
    % Calculate group birth rates.
    invest = 1 - s_g + s_g * (z_group_t_y(idx_centres, :, 5));
    avgInvest = 1 - s_g + s_g * (z_bar_group_t(5));
    b_group_t_y(idx_centres, :) = invest ./ avgInvest;
    
    % Calculate group production rates.    
    for i = 1:num_group_types                
        germ = sum(sum(b_group_t_y(idx_centres, :) .* n_group_t_y(idx_centres, :) .* squeeze(H_group_t_y(num_ghosts + 1, :, i))));
        no_germ = sum(sum(b_group_t_y(idx_centres, :) .* n_group_t_y(idx_centres, :) .* squeeze(H_group_t_y(idx_centres, :, i))));
        p_group_t_y(idx_descendants, i) = numu * germ + (1 - numu) * no_germ;
    end
    p_group_t_y = p_group_t_y / num_descendantCells;
    
    % Calculate the death rates.
    invest = 1 - s_g + s_g * (1 - z_group_t_y(idx_centres, :, 4));
    avgInvest = 1 - s_g + s_g * (1 - z_bar_group_t(4));
    d_group_t_y(idx_centres, :) = N_group_t / K_g * (invest / avgInvest);
    
    %% Update boundary conditions.
    
    % Calculate the new boundary (outflow boundary conditions).
    n_group_t_y(1, :) = n_group_t_y(num_ghosts + 1, :);
    n_group_t_y(2, :) = n_group_t_y(num_ghosts + 1, :);
    n_group_t_y(num_ghosts + num_centres + 1, :) = n_group_t_y(num_ghosts + num_centres, :); % Outflow boundary from LeVeque (2002; pg 131)
    n_group_t_y(num_ghosts + num_centres + num_ghosts, :) = n_group_t_y(num_ghosts + num_centres, :); % Outflow boundary from LeVeque (2002; pg 131)
    
    %% Update the group densities
    
    % Calculate the new group density from the between-group death dynamics.
    n_group_t_y(idx_centres, :) = n_group_t_y(idx_centres, :) + dt * (- rho * n_group_t_y(idx_centres, :) .* d_group_t_y(idx_centres, :));
    
    % Update slope function
    r_edges(idx_edges, :) = indicator_edges(idx_edges, :) .* ((n_group_t_y(idx_edges - 1, :) - n_group_t_y(idx_edges - 2, :)) ./ (n_group_t_y(idx_edges, :) - n_group_t_y(idx_edges - 1, :))) + ~indicator_edges(idx_edges, :) .* ((n_group_t_y(idx_edges + 1, :) - n_group_t_y(idx_edges, :)) ./ (n_group_t_y(idx_edges, :) - n_group_t_y(idx_edges - 1, :)));
    
    phi_edges_donor_cell(idx_edges, :) = 0; % Upwind difference method  
    phi_edges_superbee(idx_edges, :) = max(0, max(min(1, 2 * r_edges(idx_edges, :)), min(2, r_edges(idx_edges, :)))); % Superbee method
    phi_edges_superbee(isnan(phi_edges_superbee)) = 2; % limiting value when r --> inf   
    phi_method = phi_edges_superbee;
    
    % Update flux
    flux_edges(idx_edges, :) = (1 / 2) .* dy_dt_edges(idx_edges, :) .* ((1 + flip_flop_edges(idx_edges, :)) .* n_group_t_y(idx_edges - 1, :) + (1 - flip_flop_edges(idx_edges, :)) .* n_group_t_y(idx_edges, :)) + (1 / 2) .* abs(dy_dt_edges(idx_edges, :)) .* (1 - abs((dy_dt_edges(idx_edges, :) .* dt) / dy)) .* phi_method(idx_edges, :) .* (n_group_t_y(idx_edges, :) - n_group_t_y(idx_edges - 1, :));
    
    % Total flux should be equal to difference between flux at two
    % endpoints - check these two values are equal.
    totalFlux1 = sum(diff(flux_edges(idx_edges, :)));
    totalFlux2 = flux_edges(num_ghosts + 1, :) + flux_edges(num_ghosts + num_centres + 1, :);
    
    % Calculate the new group density from the within-group dynamics.
    n_group_t_y(idx_centres, :) = n_group_t_y(idx_centres, :) - (dt / dy) * diff(flux_edges(idx_edges, :));
    
    % Calculate the new group density from the between-group production dynamics.
    n_group_t_y(idx_centres, :) = n_group_t_y(idx_centres, :) + dt * (rho * p_group_t_y(idx_centres, :));
    
    % Increment time.
    t = t + dt;
    
    %% Update the results.
    
    % Store outputs.
    if t >= t_step
        
        % Next time step.
        t_step = t_step + dy;
        
        % Store group density outputs.
        n_group_output(t_step, :, :) = n_group_t_y;
        
        % Store average trait expression outputs.
        z_bar_group_output(t_step, :) = z_bar_group_t;
        
    end        
    
    % Break if system reaches steady state.   
    if t_step > 5000 && t_step > lambda * 10 && sum(sum(abs(n_group_output(t_step, :, :) - n_group_output(t_step - lambda * 10, :, :)) <= 1e-6) == num_centres + 2 * num_ghosts) == num_group_types
        disp(['Final time: ', num2str(t_step)]);
        break
    end
    
%     % Print outputs for debugging.
%     if mod(t_step, round(lambda_group)) == 0
%         
%         t
%         N_group_t
%         n_group_t = sum(n_group_t_y(idx_centres, :)) * dy
%         z_bar_group_t
% 
%     end
    
end

% Prepare output.
n_group_t_y = n_group_output(sum(n_group_output, [2, 3]) > 0, idx_centres, :);
t_final = size(n_group_t_y, 1);
z_bar_group_t = z_bar_group_output(1:t_final, :);

% Get individual perspective.
n_t = zeros(t_final, num_types);
z_bar_t = zeros(t_final, num_traits);

for k = 1:t_final

    % Calculate individual quantities
    for i = 1:num_group_types % every group type
        for y = 1:num_centres % every age class                 
            for j = 1:num_types % every cell type
                n_t(k, j) = n_t(k, j) + n_group_t_y(k, y, i) .* n_y(y, i, j);
                z_bar_t(k, :) = z_bar_t(k, :) + n_group_t_y(k, y, i) .* n_y(y, i, j) * Z(j, :);
            end
        end
    end
        
    z_bar_t(k, :) = z_bar_t(k, :) / squeeze(sum(n_t(k, :), 2));
            
end


end