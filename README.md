All code involving the statistical models and tests is in this repository. Each directory is summarized below, and there is an additional readme in each directory explaining each file.

Methodologies:
clustering_kmeans/ - Python code for k-means clustering (all other code is R). Includes the silhouette score in fig_s1. 
exampleAnalysis/ - for those wanting to analyze their own data. 

Main text figures:
fig_3/ - simulated data and analysis for the gates derived from k-means clustering. 
fig_4/ - real data for the two-input gates. Also the data for fig_s5, a re-representation of fig_4.  
fig_5/ - multi-input cases. 

Supporting figures:
fig_s2_s3/ - simulated data and analysis for all eight two-input logics, when the true beta is 0.45, and Monte Carlo power analysis for these conditions. 
fig_s4/ - statistics on additional real two-input gate data. 
fig_s6_s7_s8_s9/ - Monte Carlo power analyses for the 6 k-means example simulated gates. 
fig_s10/ - SI figure showing that rerunning the data in fig3/ doesn't change the stats. 