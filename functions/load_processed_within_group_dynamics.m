function [ n_y, N_y, x_y, z_bar_y, fe_b_output, fe_d_output, fe_count_output ] = load_processed_within_group_dynamics(two_cell_experiment, cooperation_experiment, t_end, rho, K, s_c, phi, mu, nu)

fe_b_output = [];
fe_d_output = [];
fe_count_output = [];

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

% Not every processed file will have mutation data.
try
    file = load(['processed_cell_results/', file_list{1}], 'fe_b_output', 'fe_d_output', 'fe_count_output');
    fe_b_output = file.fe_b_output;
    fe_d_output = file.fe_d_output;
    fe_count_output = file.fe_count_output;
catch
    
end

end