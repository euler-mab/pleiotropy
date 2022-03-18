function [n_output, fe_b_output, fe_count_output, fe_d_output, replicates] = load_raw_within_group_dynamics(verbose, two_cell_experiment, cooperation_experiment, t_end, rho, K, s_c, phi, mu, nu)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n_output = [];
fe_b_output = [];
fe_d_output = [];
fe_count_output = [];
replicates = 0;

desired_replicates = 10000;
if phi == 0
    desired_replicates = 10000;
end

listing = dir('cell_results/');

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

if length(file_list) > 0

    replicates = 0;

    for i = 1:length(file_list)

        % Only need a small number of examples to compute fitness effects.
        if i == 1
            if verbose == 1
                file1 = load(['cell_results/', file_list{i}], 'n_output', 'fe_b_output', 'fe_count_output', 'fe_d_output');
                n_output = file1.n_output;
                fe_b_output = file1.fe_b_output;
                fe_count_output = file1.fe_count_output;
                fe_d_output = file1.fe_d_output;
            else
                file1 = load([ 'cell_results/', file_list{i}], 'n_output');
                n_output = file1.n_output;
            end
        else
            filen = load([ 'cell_results/', file_list{i}], 'n_output');
            n_output = cat(3, n_output, filen.n_output);
        end

        % Store the number of replicates.
        replicates = size(n_output, 3);
        if replicates >= desired_replicates
            n_output = n_output(:, :, 1:desired_replicates, :, :);
            replicates = size(n_output, 3);
            break;
        end

    end

else

    disp(strcat('No files loaded!'))

end

end
