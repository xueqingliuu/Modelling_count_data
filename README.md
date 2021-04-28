# Modelling_count_data

Simu_count data.R is the main R script for conducting the simulation study.

## Parameters
- simu: number of simulations
- I: number of participants
- T0: number of study days
- sigma.a: parameter for the subject-specific intercept
- alpha: subject-specific intercept that follows N(0, sigma.a^2)
- beta.t: coefficient of the study day
- beta.M: coefficients of the motivational messages
- beta.F: coefficients of the feedback messages
- true_para: parameter vector
- a: extra-dispersion parameter

##Data Simulation
- Reward: daily step count, each column indicates a simulated dataset
- Mot: motivational messages, categorical
- Feed: feedback messages, categorical
- t0: study day
- id: subject
- Data: a list containing the simulated data


##Analysis
- Estresults: a function for calculating estimates, bias, RMSE, SE and relative bias.
- GetAnalyses: a function that helps obtain all performance results
- results: a list of comparison results