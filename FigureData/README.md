# Numerical values represented in manuscript figures

## Files included

| Figure  | Panel | File name                   | Description                                                                                                             |
|---------|-------|-----------------------------|-------------------------------------------------------------------------------------------------------------------------|
| Fig 1   | A     | `fig1_a.csv`                | AUC values for replicate train/calibrate/test splits                                                                    |
|         | B     | `fig1_b_replicates.csv.zip` | Predictions from 1000 replicate train/calibrate/test splits, used to generate AUC curves                                 |
|         |       | `fig1_b_bagged.csv`         | Predictions from the bagged model, used to generate the summary AUC curve                                                |
|         | C     | `fig1_c.csv`                | Confusion matrix for the best model, derived from S1 Table / bagged predictions                                         |
|         | D     | `fig1_d.csv`                | Viruses encountered and cumulative success when screening viruses in the order recommended by the best model            |
| Fig 2   | A     | `fig2_a.csv`                | SHAP values for each virus-genome feature combination, used to cluster viruses by explanation simililarity              |
|         | B     | `fig2_b.csv`                | Mean absolute SHAP values across all viruses, as a measure of feature importance                                        |
|         | C     | `fig2_c.csv`                | Feature importance and relative ranks for features measuring the same quantity in both referenced and unreferenced form |
|         | D     | `fig2_d_clusters.csv`       | Relative importance of different clusters of correlated features (bar heights and error-bars in Fig 2D)                 |
|         |       | `fig2_d_featuresets.csv`    | Composition of clusters (colours in Fig 2D)                                                                             |
| Fig 3   |       | --                          | See data for S11 Fig; more detailed host/vector information recorded in [InternalData/NovelVirus_Hosts_Curated.csv](../InternalData/NovelVirus_Hosts_Curated.csv) |
| Fig 4   |       | --                          | These plots show partial residuals from a model fit - see the plotting script linked below                              |
| Fig 5   | A     | --                          | Values included in S1 Table                                                                                             |
|         | B     | `fig5_b.csv`                | Predictions for different *Sarbecovirus* genomes                                                                        |
| S1 Fig  |       | `s1_fig.csv`                | List of viruses used for model training, including taxonomy                                                             |
| S2 Fig  | A     | `s2_fig_a.csv`              | Predictions from the taxonomy-based model                                                                               |
|         | B     | `s2_fig_b.csv`              | Predictions from the phylogenetic neighbourhood-based model                                                             |
| S3 Fig  |       | `s3_fig.csv`                | Family-specific AUC values, measuring the probability of accurately ranking known human-infecting viruses above other viruses from the same family |
| S4 Fig  |       | --                          | Values included in S1 Table                                                                                             |
| S5 Fig  |       | --                          | Values included in S1 Table                                                                                             |
| S6 Fig  |       | `s6_fig_taxonomy.csv`       | Cumulative success when screening viruses in the order recommended by the taxonomy-based model                          |
|         |       | `s6_fig_pn.csv`             | Cumulative success when screening viruses in the order recommended by the phylogenetic neighbourhood-based model        |
|         |       | *                           | Note that this figure also displays values from Fig 1D (`fig1_d.csv`) for comparison                                    |
| S7 Fig  |       | --                          | Clustering based on data in S10 Fig, see plotting script linked below                                                                                        |
| S8 Fig  |       | `s8_fig.csv`                | Clustering similarity comparisons using the Fowlkes-Mallows (Bk) index                                                  |
| S9 Fig  |       | `s9_fig.csv`                | Membership of correlated-feature clusters and coordinates in 2 dimensions obtained via multi-dimensional scaling        |
| S10 Fig |       | `s10_fig.csv.zip`           | SHAP values for each virus-genome feature combination, along with the underlying feature values                         |
| S11 Fig |       | `s11_fig.csv`               | Predictions for holdout viruses (also included in S1 Table)                                                             |
| S12 Fig |       | --                          | See data for S10 Fig                                                                                                    |
| S13 Fig |       | `s13_fig.csv`               | AUC values for replicate training rounds with varying numbers of input features                                         |


## Interpretation
For the interpretation of these values, see the manuscript and the underlying code. All panels of a given figure are produced by the same script, as detailed below:

| Figure  | Code                                                                                                        | Command to reproduce figure and values<sup>[1](#footnote1)</sup> |
|---------|-------------------------------------------------------------------------------------------------------------|----------------------------------------|
| Fig 1   | [Scripts/Plotting/MakeFigure1.R](../Scripts/Plotting/MakeFigure1.R)                                         | `make Plots/Figure1.pdf`               |
| Fig 2   | [Scripts/Plotting/MakeFigure2.R](../Scripts/Plotting/MakeFigure2.R)                                         | `make Plots/Figure2.pdf`               |
| Fig 3   | [Scripts/Plotting/MakeFigure3.R](../Scripts/Plotting/MakeFigure3.R)                                         | `make Plots/Figure3.pdf`               |
| Fig 4   | [Scripts/Plotting/MakeFigure4.R](../Scripts/Plotting/MakeFigure4.R)                                         | `make Plots/Figure4.pdf`               |
| Fig 5   | [Scripts/Plotting/MakeFigure5.R](../Scripts/Plotting/MakeFigure5.R)                                         | `make Plots/Figure5.pdf`               |
| S1 Fig  | [Scripts/Plotting/MakeSupplementaryFigure_RawData.R](../Scripts/Plotting/MakeSupplementaryFigure_RawData.R) | `make Plots/Supplement_RawData.pdf`    |
| S2 Fig  | [Scripts/Plotting/MakeSupplement_RelatednessModelRanks.R](../Scripts/Plotting/MakeSupplement_RelatednessModelRanks.R) | `make Plots/Supplement_RelatednessModelRanks.pdf` |
| S3 Fig  | [Scripts/Plotting/MakeSupplement_FamilyAUC.R](../Scripts/Plotting/MakeSupplement_FamilyAUC.R)               | `make Plots/Supplement_family_auc.pdf` |
| S4 Fig  | [Scripts/Plotting/MakeSupplement_TrainingSetRanks.R](../Scripts/Plotting/MakeSupplement_TrainingSetRanks.R) | `make Plots/Supplement_TrainingSetRanks.pdf` |
| S5 Fig  | [Scripts/Plotting/MakeSupplement_HighPriority_MissedZoonoses.R](../Scripts/Plotting/MakeSupplement_HighPriority_MissedZoonoses.R) | `make Plots/Supplement_HighPriority_MissingZoonoses.pdf` |
| S6 Fig  | [Scripts/Plotting/MakeSupplement_ScreeningSuccessRate.R](../Scripts/Plotting/MakeSupplement_ScreeningSuccessRate.R) | `make Plots/Supplement_ScreeningSuccessRate.pdf` |
| S7 Fig  | [Scripts/Plotting/MakeSupplementaryFigure_ClustersVsTaxonomy.R](../Scripts/Plotting/MakeSupplementaryFigure_ClustersVsTaxonomy.R) | `make Plots/Combine_tanglegrams.pdf` |
| S8 Fig  | [Scripts/Plotting/MakeSupplementaryFigure_ClustersVsTaxonomy.R](../Scripts/Plotting/MakeSupplementaryFigure_ClustersVsTaxonomy.R) | `make Plots/Supplement_bk_plots.pdf` |
| S9 Fig  | [Scripts/Plotting/MakeSupplementaryFigure_FeatureClusters.R](../Scripts/Plotting/MakeSupplementaryFigure_FeatureClusters.R) | `make Plots/SupplementaryFigure_FeatureClusters.pdf` |
| S10 Fig | [Scripts/Plotting/MakeSupplementaryFigure_EffectDirection.R](../Scripts/Plotting/MakeSupplementaryFigure_EffectDirection.R) | `make Plots/SupplementaryFigure_EffectDirection.pdf` |
| S11 Fig | [Scripts/Plotting/MakeSupplement_NovelVirusHosts.R](../Scripts/Plotting/MakeSupplement_NovelVirusHosts.R)   | `make Plots/Supplement_NovelVirus_Hosts.pdf` |
| S12 Fig | [Scripts/Plotting/Supplement_IllustrateDerivedGenomeFeatureCalcs.R](../Scripts/Plotting/Supplement_IllustrateDerivedGenomeFeatureCalcs.R) | `make Plots/Supplement_methods_derived_genome_features.pdf` |
| S13 Fig | [Scripts/Plotting/MakeSupplement_FeatureSelection.R](../Scripts/Plotting/MakeSupplement_FeatureSelection.R) | `make Plots/Supplement_FeatureSelection.pdf` |


<sup name="footnote1">1</sup>See the [main README](../README.md#requirements) for setup requirements
