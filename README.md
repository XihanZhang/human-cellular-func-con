# human-cellular-func-con

This repo provide analysis-ready cell type abundance maps and analysis code for [The cell-type underpinnings of human functional cortical connectome (Nature Neuroscience 2024)](https://doi.org/10.1038/s41593-024-01812-2) 

Zhang, X. H.(张喜寒), Anderson, K. M., Dong, H. M.（董昊铭）, Chopra, S., Dhamala, E., Emani, P. S., Grestein, M. B., Margulies, D. S., & Holmes, A. J. (2024). The cell-type underpinnings of human functional cortical connectome. (Nature Neuroscience 2024)

Keywords: Functional gradients; Functional networks; Cortical organization; fMRI; Genetics; Cell-type imputation; Inhibitory interneurons; Excitatory neurons; Gene transcription

## Analysis-ready cell maps
Please cite the above paper for the use of cell maps and code:

Below is the structure of the cell maps:
```
--cell maps
      |--vertex_level
              |--'Jorstad_cell_type_fractions_vertex_level.csv': cell fractions per AHBA sample were imputed from Jorstad dataset. 'well_id' is the AHBA sample ID.
              |--'sample_info_vertex_reannot_mapped_0.3.csv': AHBA sample information are stored in this table. Use the column 'well_id' to match the other information like MNI coordinates, or vertex id.
      |--Schaefer_atalas
                |--Jorstad
                      |--donor_level_cell_maps: cell fractions imputed from Jorstad dataset and each individual AHBA donor are aggerated in Schaefer atalses.
                      |--cell fractions imputed from Jorstad dataset and combining all AHBA donors are aggerated in Schaefer atalses (100-1000 parces).
                |--Lake
                      |--donor_level_cell_maps: cell fractions imputed from Lake dataset and each individual AHBA donor are aggerated in Schaefer atalses.
                      |--cell fractions imputed from Lake dataset and combining all AHBA donors are aggerated in Schaefer atalses (100-1000 parces).
```
### Cell-type abundaces in Schaefer 100/200/300/400/500/600/700/800/900/1000 parcellation are provided in `/cell_maps/Schaefer_atalas/Jorstad` and `/cell_maps/Schaefer_atalas/Lake`.
e.g. for group-aggregated cell type abundances in schaefer 400 and 7 network parcellation, files are:
- `schaeffer_Jorstad_400_7Net_expr_mat_new_NormZscore0.3.csv`, imputed from [Jorstad 2023](https://www.science.org/doi/10.1126/science.adf6812) and AHBA processed with `ibf=0.3`, `normalization=zscore`.
- `schaeffer_LAKE_DFC_400_7Net_expr_mat_new_NormZscore0.3.csv`, imputed from Lake_DFC and AHBA processed with `ibf=0.3`, `normalization=zscore`.
- `schaeffer_LAKE_VIS_400_7Net_expr_mat_new_NormZscore0.3.csv`, imputed from Lake_VIS and AHBA processed with `ibf=0.3`, `normalization=zscore`.
- `schaefer_400_7Net_labels.csv`, for the corresponding labels of the networks.
- `schaeffer400_7_numDonorsInEachParcel_abagen_NormZscore0.3.csv`, the number of donors contributing to each parcel.
- `schaeffer400_7_maxDonorPresenceInParcel_abagen_NormZscore0.3.csv`, the dominant donor within each parcel.
and the donor level cell type abundances in schaefer 400 and 7 network parcellations are within `/cell_maps/Schaefer_atalas/Jorstad/donor_level_cell_maps` and `/cell_maps/Schaefer_atalas/Lake/donor_level_cell_maps`
Toy code to visualize the group-aggregated cell type abundances in schaefer 400 and 7 network parcellation, and spin-test are provided for quick understanding.
- `1_CellAbundance_Visulization.ipynb`: for visualization.
- `2_spintesting-cleaned.ipynb`: for spin-testing.

## Analysis code for the paper
Please cite the above paper for the use of these code:
### Part A: processing and aligning AHBA bulk tissue samples
#### `00_aggregate_ahba.R`
- Extract and aggregate necessary information from tables downloaded from [AHBA](http://human.brain-map.org/)
- Partial data of the AHBA downloaded by abagen are provided in `/microarray/normalized_microarray_donor9861`, as an example.
- Input: downloaded AHBA raw data (specifically, the `top_level_categories.csv`, and `SampleAnnot.csv`, `Ontology.csv` for each donor)
- Output: `ahba_sampleInfo.csv`

#### `01_tutorial.ipynb` Part 1.
- Using abagen to correct the old MNI coordinates
- Input: `ahba_sampleInfo.csv`
- Output: `ahba_sampleInfo_reannot.csv`
  
#### `01_1_project_freesurfer.bash` (Also in `01_tutorial.ipynb` Part 2 Step 1)
- Project fsLR parcellations to individual AHBA freesurfer space
- Input and output files are in the `mri`, `surf`, `label` folders under each donor dirs that sit wherever the abagen download dir sit, within `/abated-data/freesurfer`.

(e.g. abagen by default download the files to `/Users/zhangxihan/abated-data/freesurfer/donor9861/mri`, but I moved the entire freesurfer folder to `/gpfs/milgram/project/holmes/xz555/gradient_shift/data/ahba_fs`)

#### `01_2_map_ahba_ctx_to_surface.py` (Also in `01_tutorial.ipynb` Part 2 Step 2)
- Map native space AHBA sample coordinates (x,y,z) onto a 32k midthickness file that is spatially aligned with native cortical geometry (vertices are aligned across individuals)
- Input: `ahba_sampleInfo_reannot.csv` + files generated from `01_1_project_freesurfer.bash`.
- Output: `sample_info_vertex_reannot_mapped.csv`

#### `01_3_preprocess_ahba_abagen.py` (Also in `01_tutorial.ipynb` Part 3 Step 1)
- Pre-process the gene under different combinations of parameters of abagen toolbox.

#### `01_4_GeneSymbolToEntrezID.R` (Also in `01_tutorial.ipynb` Part 3 Step 2)
- Some Entrez ID are missing in the AHBA tables, so in the earlier works, the probes survived from the QC but missing Entrez ID were removed. Those probes have a gene symbol with them, so in this current work, we used these gene symbol as index and pull out their corresponding Entre ID from the human genome wide annotation `org.Hs.eg.db`.

#### `01_5_make_probe_table.py` (Also in `01_tutorial.ipynb` Part 3 Step 3)
- Assemble the probe information into a table for later analysis.

#### `01_tutorial.ipynb` Part 4.
Compares the processed gene lists from different parameter combo of abagen with the gene list generated from previous published works by Anderson et al. Which validated the processed AHBA gene expressions and the imputed cell types that based upon it.
[Anderson, KM, Collins, MA, Chin, R, Ge T, Rosenberg MD, Holmes AJ. (2020). Transcriptional and imaging-genetic association of cortical interneurons, brain function, and schizophrenia risk. Nature Communications, 11, 2889](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7280213/pdf/41467_2020_Article_16710.pdf)

#### `01_6_make_abagen_object.R`
- Convert the table into R object for later analyses that run in R.

#### `02_vertex_to_schaeffer_parcel_gene.R`
- This script is not directly contributing to the current work, as the cell type deconvolution should be done at the sample level and within-donor to control for batch effect. We provide this script as a BONUS for people who want to parcel the gene expression into Schaefer atlas.
- Makes parcel-wise gene expression averages for each ROI in the 100/200/300/400/500/600/700/800/900/1000 7/17 Network parcellations of Scheaffer et al (2018).
- You can easily change the ROI and network resolution by changing the number in the script. Just remember to download the parcellation from [Git repo of Thomas Yeo group](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti), and cite [Schaeffer et al (Cerebral Cortex 2018)](https://doi.org/10.1093/cercor/bhx179). (e.g. `Schaefer2018_300Parcels_7Networks_order.dscalar.nii` and `Schaefer2018_300Parcels_7Networks_order_info.txt` for 300 ROIs of 7 networks).
- This script allows gene expression aggregation across all donors and within each donor.
- The number of donor contributing to and the dominant donor within each parcel is reported in csv files, which help you to check the donor dominance issue. (e.g. `schaeffer300_7_numDonorsInEachParcel_abagen_NormZscore0.3.csv` and `schaeffer300_7_maxDonorPresenceInParcel_abagen_NormZscore0.3.csv` for 300 ROIs of 7 networks)

#### `03_preprocess_lake_data.R`
Preprocesses raw single-cell UMI visual and frontal cortex data from [Lake et al. 2018](https://www.nature.com/articles/nbt.4038). Freely available for download [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930).
- Applies basic preprocessing steps using the Seurat package.
(1) identification of cell and gene expression outliers.
(2) log-normalization.
(3) scaling and regression of number of detected genes per cell.
- Create gene name dictionary to match AHBA genes to sn-DropSeq genes.
- Write *seurat_processed.Rdata data objects for later reading.

#### `04_cibersortx_prep.R`
For visual and frontal cortex sn-DropSeq data, subset to those genes that are present in both.
- Transform data out of log-space (format expected by CIBERSORTx).
- To reduce collinearity among gene signature matrix, collapse cell subtypes into overarching categorie (e.g. In6a and In6b cell become just PVALB).
- Write single-nuclei expression matrix, and AHBA mixture files per donor, as the cell type deconvolution should be done at the sample level and within-donor to control for batch effect. That's why we ask abagen to output the result at the sample level, instead of the parcel level. Abagen toolbox did the sample and gene normalization within each donor, so it suits our purpose well in this case.

#### Cell type deconvolution is performed at [CIBERSORTx](https://cibersortx.stanford.edu/)
- Files to feed CIBERSORTx:
- 1) donor-level `Mixture` matrix (e.g. for donor 9861): `FrontalCortex_ahba_9861_mixture_new_NormZscore0.3.txt`, `VisualCortex_ahba_9861_mixture_new_NormZscore0.3.txt`
- 2) `Single Cell Reference` matrix: `Lake_FrontalCortex_ahba_matched_sc_sample_new_NormZscore0.3.txt`, `Lake_VisualCortex_ahba_matched_sc_sample_new_NormZscore0.3.txt`
- Procedures in CIBERSORTx:
- 1) Upload the above files: `Menu` -> `Upload Files`
- 2) Calculated expression signature for Lake_DFC and Lake_VIS (taking Lake_VIS as example):
     `Menu` -> `Run CIBERSORTx` -> `Create Signature Matrix` -> `Custom` -> `scRNA-Seq` -> Single cell reference matrix file: `Lake_VisualCortex_ahba_matched_sc_sample_new_NormZscore0.3.txt` -> check `Disable quantile normalization` -> Unfold `Single Cell Input Options` -> `Min. Expression`=0.5 (or 0) -> `Run`
- 3) Impute cell type fractions (taking donor 9861 with Lake_VIS as example):
     `Menu` -> `Impute Cell Fractions` -> `Custom` -> Signature matrix file: the file derived from 2) above -> Mixture file: `VisualCortex_ahba_9861_mixture_new_NormZscore0.3.txt` -> check `Enable batch correction` -> `S-mode` -> Single cell reference matrix file: `Lake_VisualCortex_ahba_matched_sc_sample_new_NormZscore0.3.txt` -> check `Disable quantile normalization` -> check `Run in absolute mode` -> Permutations for significance analysis: None
- Download the reults and figures.

#### `05a_vertex_to_schaeffer_parcel_cell.R`
- Makes parcel-wise gene expression averages for each ROI in the 200/300/400/1000 7/17 Network parcellations of Scheaffer et al (2018).
- In our work, we used 400 parcels and 7 networks resolution. You can easily change the ROI and network resolution by changing the number in the script. Just remember to download the parcellation from [Git repo of Thomas Yeo group](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti), and cite [Schaeffer et al (Cerebral Cortex 2018)](https://doi.org/10.1093/cercor/bhx179). (e.g. `Schaefer2018_300Parcels_7Networks_order.dscalar.nii` and `Schaefer2018_300Parcels_7Networks_order_info.txt` for 300 ROIs of 7 networks).
- Also plot the gene expression signature matrices, gene expression signature correlation matrices, and the cortical cell type fraction correlation matrices of Lake_DFC and Lake_VIS. 

#### `05b_cibersortx_compare_plots.R`
- This script is the plotting part of the `05a_vertex_to_schaeffer_parcel_cell.R`. It mainly serves the purpose of making plots for the cell-type fractions imputed from AHBA processed under different combinations of parameters. You can skip this script if you just want to use the final combo for your project.

### Part B: Parceling the gradients and comparing with cell type abundance distribution
#### `06_hcp_gradient_giitonii.sh`
- Convert the gradients stored in gifti to nifti.

#### `07_SpintestHCPGradient.m`
- Spins the gradient 1 and 2 values on fsaverage5.
- Parcel the gradient values in Schaefer 400.
- Calculate the parcel-level correlation between gradient values and cell type abundances, and calculate the p-value based on the null gradients generated from spin-test.
- Requires matlab package download from [Alexander et al.](https://github.com/spin-test/spin-test).

#### `07a_VisualizationFig1_GradientNetwork.ipynb`
- Make plots for Figure 1.

#### `07b_VisualizationFig2Supp1-2_CellGradientCorrelation.ipynb`
- Make plots for Figure 2 and supplement figures 1-2.

#### `08a_PermCCA_GradientVarbyAllCell.m`
- Permulational CCA between gradients and all cell types common between Lake_DFC and Lake_VIS.
- Require Matlab software download from [Winkler et al.](https://github.com/andersonwinkler/PermCCA).

#### `08b_PermCCA_GradualRemoval.m`
- Permulational CCA between gradients and cell types with different removal combinations of gradient-correlated cell types.

#### `08c_VisualizationFig3Supp3-4_PermCCA.ipynb`
- Make plots for Figure 3 and supplement figures 3-4.

### Part C: Network-level cell type enrichment
#### `09a_SpintestCellTypeNull_DFC.py` and `09a_SpintestCellTypeNull_VIS.py`
- Prepare the cell type abundances map and use spin-test to generate null abundance distribution of each cell type.
- `09a_SpintestCellTypeNull_DFC.sh` and `09a_SpintestCellTypeNull_VIS.sh` are the scripts to run these analysis in Milgram cluster.

#### `09b_VisualizationFig4_Network_Cell_Enrichment.ipynb`
- Calculate p-values and make plots for Figure 4.

### Part D: Predicting functional network from cell type abundance
#### `10a_Spin_Network_Labels.ipynb`
- Using spin-test to generate 1000 null network labels fro each parcel.
- Also provide code to spin the cell type abundance at the parcel level, which is already included in `09a_SpintestCellTypeNull_DFC.py`. Here we include it just in case someone want to do both in one script.

#### `10b_Classifying pipeline.ipynb`
- Train the classifiers from the empirical data and the null data generated in the previous step.

#### `10c_Classifying pipeline.ipynb`
- Make plots for figure 5 and supplement figures 5-10.


