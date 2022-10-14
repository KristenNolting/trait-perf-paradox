Scripts for re-running and recreating figures for the manuscript titled: "When traits vary across species but performance doesn't: One solution to the paradox"

All figures saved in \Output\Figures

To recreate figures and results in order they are presented in the manuscript:

1. Data_Wrangling.R - this will import the original data file used in Nolting et al. (2021) and trim and clean it to create a Protea_data object that will be used in subsequent analyses for the current manuscript. 

2. multiple_regressions.R - this script will run through the models to fit the structural trait-physiological performance trait multiple regressions, demonstrating that structural traits explain variation in performance well. *Figure 1* in the manuscript is generated from this script (legend added manually). All model fit files have been saved for reproducibility of results.

3. variance_partitioning.R - this script will run through the models to perform a variance partitioning on all structural traits and on all physiological performance traits. *Figure 2* in the manuscript is generated from this script. All model fit files have been saved for reproducibility of results.

4. two-trait-model_simulation.R - this script runs through the simple two-trait model to simulate trait data under scenario 1 and scenario 2, as presented in the manuscript. *Figure 3* in the manuscript is generated from this script.

5. pc-functional-axis-alignment.R - this script runs through fitting multiple regression models as in 2 above, but without the random effects (for ease of comparison with the generalized model presented in the text). Using this model output, the script then runs through calculating the alignment between each structural trait PC axis and the functional axis corresponding to each physiological performance trait. *Figure 4* in the manuscript is generated from this script.

6. pca-performance_plots.R - this script runs through plotting the Protea_data in structural trait space, and coloring points based on physiological performance. *Figure 5* in the manuscript is generated from this script.



## All session info reported below: