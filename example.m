





verbose=1
cooperation_experiment=1
replicates=1000
t_0=0
t_by=1
t_end=600
rho=1
K=100
s=0.95
p=0
mu=0.0001
nu=0.01

[n_output, fe_count_output, fe_b_output, fe_d_output] = withinGroupDynamicsOneCell(verbose, cooperation_experiment, replicates, t_0, t_by, t_end, rho, K, s, p, mu, nu)