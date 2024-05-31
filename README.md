# Multi-omics framework to reveal the molecular determinants of fermentation performance in wine yeast populations

## About the work
### Background
Connecting the composition and function of industrial microbiomes is a major aspiration in microbial biotechnology. Here, we address this question in wine fermentation, a model system where the diversity and functioning of fermenting yeast species is determinant of the flavor and quality of the resulting wines.
### Results
First, we surveyed yeast communities associated with grape musts collected across wine appellations, revealing the importance of environmental (i.e., biogeography) and anthropic factors (i.e., farming system) in shaping community composition and structure. Then, we assayed the fermenting yeast communities in synthetic grape must under common winemaking conditions. The dominating yeast species defines the fermentation performance and metabolite profile of the resulting wines, and it is determined by the initial fungal community composition rather than the imposed fermentation conditions. Yeast dominance also had a more pronounced impact on wine meta-transcriptome than fermentation conditions. We unveiled yeast-specific transcriptomic profiles, leveraging different molecular functioning strategies in wine fermentation environments. We further studied the orthologs responsible for metabolite production, revealing modules associated with the dominance of specific yeast species. This emphasizes the unique contributions of yeast species to wine flavor, here summarized in an array of orthologs that defines the individual contribution of yeast species to wine ecosystem functioning. 
### Conclusions
Our study bridges the gap between yeast community composition and wine metabolite production, providing insights to harness diverse yeast functionalities with the final aim to producing tailored high-quality wines.

Currently available as a preprint at [bioRxiv](https://doi.org/10.1101/2023.12.02.569693)

## About the sampling effort and data analysis

### Sampling design

<img src="/Figures/Map.png" width="350" align="left"> </img>  We surveyed five distinct Spanish wine appellations, sampling grapes from vineyards under conventional and organic management. We also sampled at lower geographic scales to finely disentangle whether grape must and fungal variability was distance dependent. The file `0.Sampling_design.R` describe how to draw the sampling map presented. We collected grapes from five different grapevine plants, making a composite sample from each sampling point which were pressed under sterile conditions and the resulting grape must was dispensed into sterile bottles. 

<br clear="left"/>

### Observational study

<p align="center">
<img src="/Figures/GM.png" height=300 align="center">
</p>
<br clear="left"/>

These fermentations were sampled immediately after being dispensed, at the tumultuous stage of fermentation (between 5-50% sugars consumed) and at the end of fermentations (weight loss remained consinstently below 0.01g per day). We used each sampled fermentation point for different purposes:
  - The initial stage served as a comparison of yeast communities (**Supplementary Table S3**) across different biogeographical origins, in vineyards subject to different farming systems as represented in **Supplementary Figure S2**. Besides, we also studied grape must composition, as it widely variable across different vineyards (**Supplementary Table S1**, **Supplementary Table S2**). A through representation of these results can be found in **Figure 1**.
  - The tumultuous stage of fermentation served as a selection of fermenting yeast communities to later inoculate the laboratory fermentations.
  - With the final stage we aimed to assess which factor was more important for yeast community and wine metabolite composition at the end of fermentations. Thus, our sampling effort was focused on ITS amplicon sequencing and metabolite composition assessment (**Supplementary Figure S4**, **Supplementary Table S4**).

The files `1.GrapeMust_diversity.R` and `2.FermentedMust_metabolite.R` contains the necessary code to conduct these analysis. In addition, in the `Raw_sequences_analysis` the RScripts needed to analyse ITS raw sequences are can be found. 


### Laboratory fermentations

<p align="center">
<img src="/Figures/SGM.png" height=300 align="center">
</p>
<br clear="left"/>

These fermentations were subjected to the same fermentative conditions than the observational study. Besides, the yeast community of each control fermentation was used to inoculate this experiment. Here, we focused on the relationship between yeast community and metabolite production under homogenized media and environmental conditions. We conducted and RNAseq experiment to determine the (meta-)transcriptomic expression of each fermentation, yielding different metabolite profiles at the end of fermentations. It is important to note that we consider a fermentation ended when the community no longer consumes sugar. Even though from an enological point of view this comparison cannot be made, we are interested on the ecological significance of why some communities are able to consume all the sugars present and the metabolite profile associated. Thus, we sampled:
 
  - At the tumultuous stage to evaluate the yeast community composition (ITS sequencing) and (meta-)transcriptomic profiles (**Supplementary Figure S3**, **Supplementary Figure S7**, **Supplementary Figure S9**). Sugar content was also measured to assure that communities were at this fermentation stage (**Supplementary Table S1**).
  - At the end of fermentations to evaluate the metabolite profiles (**Supplementary Figure S5**, **Supplementary Figure S6**, **Supplementary Figure S8**).

**Figure 2** represents the yeast community compositon (assessed via RNA sequencing, represented in **Figure 3** and **Figure 4**) at the tumultuous stage already dominated by a handfull of species, and differences in the metabolite composition at the end of fermentation explained by the dominant fermentative species. **Figure 5** represents the relationship between transcriptomic activity of fermentative communities and the final metabolite composition.

The files `3.SyntheticMust_metabolite.R`, `4.Transcriptional_profiles.R`, `5.Transcriptional_profiles-Saccharomyces.R`, `6.Metabolites_transcriptomesR.R` contain all the code used for these analysis, and in the folder `Raw_sequences_analysis` you can find additional code used to analyse raw ITS and RNA sequences.

## Contact

Microbial Interactions and Ecology Lab - [MINElab](http://minelab.bioucm.es/)

[@MigueldCe](https://twitter.com/MigueldCe) - migueldc@ucm.es

[https://github.com/Migueldc1/Wineteractions](https://github.com/Migueldc1/Wineteractions/)
