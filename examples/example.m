%% Parameters

% within-group

verbose = 1; % describe which file is running

cooperation_experiment = 1; % assume private + private triats (cooperation_experiment = 0) or private + public (cooperation_experiment = 1)
two_cell_experiment = 0; % assume 1 founding cell per group (two_cell_experiment = 0) or 2 (two_cell_experiment = 1)

replicates = 3; % number of replicate stochastic simulations

t_0 = 0; % start time within a group (age 0)
t_by = 1; % resolution for recording within-group state
t_end = 600; % end time within a group (maximum age)

rho = 1; %
K = [100, 200, 500]; % group size
s_c = 0.95; % strength of within-group selection
phi = [0:0.1:1]; % strength of pleiotropy

mu = 0.0001; % loss-of-function mutation rate of traits
nu = 0.01; % relative rate of gain-of-function mutations to loss-of-function mutations

num_traits = 3; % number of binary traits
num_types = 2^num_traits; % number of cell types resulting from binary traits

Z = [
    0, 0, 0;
    0, 0, 1;
    0, 1, 0;
    0, 1, 1;
    1, 0, 0;
    1, 0, 1;
    1, 1, 0;
    1, 1, 1
    ];

% between-group

num_group_traits = 5;
num_group_types = num_types;
if two_cell_experiment == 1
    num_group_types = num_types * (num_types - 1) / 2 + num_types;
end

lambda = [10:5:50]; % parameter lambda in the text, affecting expected lifespan of groups
K_g = 1; % carrying capacity of group popularion (normalised to 1 in this workflow)
s_g = 0.95; % strength of between-group selection
c_g = 0.0; % cost of pleiotropy
gamma = 0; % parameter gamma in the text, affecting germ line vs. somatic cell
alpha = 0.0; % parameter alpha in the text, affecting age of reproductive maturity.

%% 1. Run within-group stochastic dynamics

% for i = 1:length(K) % different group sizes
%     for j = 1:length(phi) % different strengths of pleiotropy
%         
%         % run the stochastic simulations
%         if two_cell_experiment == 0
%             [n_output, fe_count_output, fe_b_output, fe_d_output] = within_group_dynamics_one_cell(verbose, cooperation_experiment, replicates, t_0, t_by, t_end, rho, K(i), s_c, phi(j), mu, nu, num_traits, num_types, Z);
%         else
%             %             [n_output, fe_count_output, fe_b_output, fe_d_output] = within_group_dynamics_two_cell(verbose, cooperation_experiment, replicates, t_0, t_by, t_end, rho, K(i), s_c, phi(j), mu, nu, num_traits, num_types, Z);
%         end
% 
%         % process the stochastic simulation outputs by averaging
%         process_raw_within_group_dynamics(verbose, two_cell_experiment, cooperation_experiment, t_end, rho, K(i), s_c, phi(j), mu, nu, num_traits, num_types, num_group_types, Z);
% 
%     end
% end

%% 2. Run between-group deterministic dynamics

for i = 1:length(K) % different group sizes
    for j = 1:length(phi) % different strengths of pleiotropy
        for k = 1:length(lambda) % different group lifespans

            % save the between-group dynamics results
            filename_between_output = sprintf('group_results/lamg%g_Kg%d_sg%g_alphg%g_gamg%g_costg%g_two%d_coop%d_yend%d_rho%g_K%d_s%g_p%g_mu%g_nu%g', ...
                lambda(k), K_g, s_g, alpha, gamma, c_g, two_cell_experiment, cooperation_experiment, t_end, rho, K(i), s_c, phi(j), mu, nu);
            filename_between_output = strrep(filename_between_output, '.', '-');
            filename_between_output = strcat(filename_between_output, '.mat');

            if ~isfile(filename_between_output)

                % run the pde solver
                [ n_group_t_y, z_bar_group_t, n_t, z_bar_t ] = between_group_dynamics(two_cell_experiment, lambda(k), K_g, s_g, c_g, alpha, cooperation_experiment, t_0, t_by, t_end, rho, K(i), s_c, phi(j), mu, nu, num_traits, num_types, num_group_traits, num_group_types, Z);
                save(filename_between_output, 'n_group_t_y', 'z_bar_group_t', 'n_t', 'z_bar_t', '-v7.3');

            end
        end
    end
end