# pleiotropy

This repository holds the code for the paper: Bentley et al (2022) Pleiotropic constraints promote the evolution of cooperation in cellular groups.

# Introduction

The space of possible states in multi-level selection models is often very large (many cells in every group \* many groups = large state space). This renders all but the most heavily simplified methods used to study multi-level selection intractable to analyse.

Bentley et al (2022) outlines a novel approach to study multi-level selection by using stochastic simulations embedded inside a partial differential equation (PDE). This gives rise to an age-structured model of multi-level selection, where within-group selection (among cells) is modelled using Gillespie's direct method, and between-group selection (among groups founded by different types of cells) is modelled as the numeric solution to a PDE.

The within-group dynamics are modelled using stochastic simulations because the invasion of mutant cell lines in small populations is stochastic process that is poorly characterised by deterministic methods. A key feature of stochastic simulations in the present context is that they capture the extended _time_ it takes for invasion to happen in some populations over others. This is a crucial feature to capture in models of the evolution of pleiotropy because it is this delayed _rate_ with which invasion happens in some groups over others which helps explain why natural selection can favour the evolution of pleiotropy, and with it, cooperation.

# Minimal Example

The approach outlined in the paper above simplifies greatly the numerical requirements needed to study complete multi-level selection dynamics, but the process of demonstrating the potential for the evolution of pleiotropy in multi-level seleciton models still requires a large amount of computation. For this reason, the steps outlined below run the user through a simplified, minimal example that can be run on a desktop computer in a few hours. A more complete analysis requires access to a super computer.

**Important Note**: the results from this example are subject to stochastic variation across runs. This is because the scale of the analysis has been reduced in order to allow it to run on a desktop computer (the figures in the MS are generated with a super computer). In practice, a large number of replicates are needed to capture fully the extent to which mutant invasion happens at different rates in different types of groups.

The example focusses on a sweep of three parameters from the paper: the strength of pleiotropy, `phi=[0:0.5:1]`, the group size, `K=[100,200,500]`, and the group lifespan, `lambda = [10:5:50]`. This minimal set of parameters allows us to recreate the three main types of figures in the paper. We will create an equivalent of:

- Fig 3, showing the average stochastic dynamics of different types of cells and different frequencies of traits within a group over time for a given set of parameters;
- Fig 4, showing a heatmap of the final trait values natural selection favours in the global population of cells when both within-group and between-group dynamics are in play for a sweep of parameter settings;
- Fig 5, showing the change in global density of different types of cells and different frequencies of traits over time for a given set of parameters.

# Running the Example

The example needs running in several stages. First, the within-group stochastic simulations need running and processing (averaging). Second, the between-group numerical solutions of the PDE need running. Finally, the results can be plotted.

## Within-group Dynamics

The functions `within_group_dynamics_one_cell` and `within_group_dynamics_two_cell` are used to simulate stochastic within-group dynamics when a group is founded by a single cell or two cells, respectively. There are three binary cell trait types in the model, which gives rise to 8 different cell types. When groups are founded by a single cell, there are also 8 group types. When groups are founded by two cells, there are 8 \* (8 - 1) / 2 + 8 = 36 group types (one for each unique pair of cell types).

The functions above simulate stochastic dynamics for each group type for a specified number of replicate runs. In production of the paper, 10,000 replicates are used. Here, we use only 3 replicates. This means the example can be run in only a few hours, but that its results will be variable between runs, and only a subset of the data needed for the full analysis shown in the paper.

In `example.m`, we assume sweep over a group size of `K=[100,200, 500]`, and a strength of pleiotropy `phi=[0:0.5:1]`. We assume a single cell founder for each type of group. Thus, in our demo, we do a total of 3 (replicates) \* 8 (types) \* 3 (group sizes) \* 11 (strengths of pleiotropy) = 792 simulations. In verbose mode, the console output looks like this:

    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (1) Replicate: 1
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (1) Replicate: 2
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (1) Replicate: 3
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (2) Replicate: 1
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (2) Replicate: 2
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (2) Replicate: 3
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (3) Replicate: 1
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (3) Replicate: 2
    ID: 802cd709-a120-4f50-a6dc-249ef249947e  Type: (3) Replicate: 3

... and so on. The ID is a unique identifier for each run that is included in the output filename along with a summary of the parameters used e.g. `two0_coop1_rep3_tend600_rho1_K100_s0-95_p0_mu0-0001_nu0-01_id_c24535a7-6a73-42b8-bf1c-202da47be97c.mat`. It allows a user to run identical parameter sweeps on different computers, and combine the files later during processing.

The outputs from simulations are stored in matrix files in the `cell_results` directory. They are then averaged together in a processing step, to give a summary of the expected dynamics within each type of group. These results are stored in the `processed_cell_results` directory. The processing step is performed by a script `process_raw_within_group_dynamics`. It uses another script `load_raw_within_group_dynamics` to load outputs with the same parameter settings into a single larger object, and then combines them by averaging over replicates.

## Between-group Dynamics

The function `between_group_dynamics` is used to solve a system of PDEs corresponding to the age-structured model of between-group competition. It loads the data from `processed_cell_results` which it then uses to characterise the typical development of a particular type of group as it ages. It then explores how group-level selection favours different types of groups, given the initial conditions and the expected behaviour within different types of groups as they age.

As a general rule, we find that group-level selection in the PDE favours groups founded by cells with pleiotropy because those groups get invaded by mutant
cheater cell lineages than their competitor groups without pleiotropy. Thus, pleiotropy gives groups a competitive advantage and can be favoured to stabilise cooperation over evolutionary timescales in conditions in which it might otherwise break down.

## Plotting the Results

The function `plot_within_group_dynamics` is used to produce figures like Fig. 3 in the main paper. It produces a plot showing the expected average of the stochastic dynamics of different types of cells and different frequencies of traits within a group over time.

The function `plot_between_group_heatmap` is used to produce figures like Fig. 4 in the main paper. It uses a sweep of a key target variable - the group size, `K`, in the example - to produce heatmaps showing the final trait values selected for when both within-group and between-group dynamics are in play, in different regions of parameter space.

The function `plot_between_group_dynamics` is used to produce figures like Fig. 5 in the main paper. It plots the change in global density of different types of cells and different frequencies of traits over time.
