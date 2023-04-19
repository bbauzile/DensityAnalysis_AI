# DensityAnalysis_AI

These codes are the online material for the article titled:
Impact of palmiped farm density on the resilience of the poultry sector to highly pathogenic avian influenza H5N8 in France

Billy Bauzile; Benoit Durand; Sébastien Lambert; Séverine Rautureau; Lisa Fourtune; Claire Guinat; Alessio Andronico; Simon Cauchemez; Mathilde Paul; Timothee Vergne

Correspondance author : timothee.vergne@envt.fr

## Purpose
We analysed the interplay between palmiped farm density and the vulnerability of the production system to highly pathogenic avian influenza (HPAI) H5N8. To do so, we used a spatially-explicit transmission model, which was calibrated to reproduce the observed spatio-temporal distribution of outbreaks in France during the 2016-2017 epidemic of HPAI.  In silico simulations suggested that reducing palmiped farm density, even slightly, in the densest municipalities was expected to decrease substantially the number of affected farms and therefore provide an intersectorial benefit. However, they also suggest that it would not have been sufficient, even in combination with the intervention measures implemented during the 2016-2017 epidemic, to completely prevent the virus from spreading. Therefore, the effectiveness of alternative structural preventive approaches needs to be assessed, including flock size reduction and targeted vaccination.

## This repository has two main parts:
### Model codes
Initially created for the analysis described in Andronico et al. (2019). In this paper, a mathematical model was developed to analyze the spatiotemporal
evolution of the 2016-2017 HPAI H5N8 epidemic in France and evaluate the impact of control strategies. 
In our analysis, the model was adjusted to investigate the impact of reducing the density of palmiped farms in municipalities with the highest duck farm density and time-varying effective reproduction numbers.
#### Makefile
To optmize the cpp simulation

### Figures
Map density, R0 and Re contain the codes that generate the figures in the article.
