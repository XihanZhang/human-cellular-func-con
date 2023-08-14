# human-cellular-func-con

Analysis code for [The cellular underpinning of human cortical functional connectome (under review)](https://www.biorxiv.org/content/10.1101/2023.07.05.547828v1.abstract)

## Part A: processing and aligning AHBA bulk tissue samples
### `00_aggregate_ahba.R`
- Extract and aggregate necessary information from tables downloaded from [AHBA](http://human.brain-map.org/)
- Partial data of the AHBA downloaded by abagen are provided in `/microarray/normalized_microarray_donor9861`, as an example.
- Input: downloaded AHBA raw data (specifically, the `top_level_categories.csv`, and `SampleAnnot.csv`, `Ontology.csv` for each donor)
- Output: `ahba_sampleInfo.csv`

### `01_tutorial.ipynb` Part 1.
- Using abagen to correct the old MNI coordinates
- Input: `ahba_sampleInfo.csv`
- Output: `ahba_sampleInfo_reannot.csv`
  
### `01_1_project_freesurfer.bash` (Also in `01_tutorial.ipynb` Part 2 Step 1)
- Project fsLR parcellations to individual AHBA freesurfer space
- Input and output files are in the `mri`, `surf`, `label` folders under each donor dirs that sit wherever the abagen download dir sit, within `/abated-data/freesurfer`.

(e.g. abagen by default download the files to `/Users/zhangxihan/abated-data/freesurfer/donor9861/mri`, but I moved the entire freesurfer folder to `/gpfs/milgram/project/holmes/xz555/gradient_shift/data/ahba_fs`)

### `01_2_map_ahba_ctx_to_surface.py` (Also in `01_tutorial.ipynb` Part 2 Step 2)
- Map native space AHBA sample coordinates (x,y,z) onto a 32k midthickness file that is spatially aligned with native cortical geometry (vertices are aligned across individuals)
- Input: `ahba_sampleInfo_reannot.csv` + files generated from `01_1_project_freesurfer.bash`.
- Output: `sample_info_vertex_reannot_mapped.csv`

### `01_3_preprocess_ahba_abagen.py` (Also in `01_tutorial.ipynb` Part 3 Step 1)
- Pre-process the gene under different combinations of parameters of abagen toolbox.

### `01_4_GeneSymbolToEntrezID.R` (Also in `01_tutorial.ipynb` Part 3 Step 2)
- Some Entrez ID are missing in the AHBA tables, so in the earlier works, the probes survived from the QC but missing Entrez ID are removed. Those probes have a gene symbol with them, so we used these gene symbol as index and pull out their corresponding Entre ID from the human genome wide annotation `org.Hs.eg.db`.

### `01_5_make_probe_table.py` (Also in `01_tutorial.ipynb` Part 3 Step 3)
- Assemble the probe information into a table for later analysis.

### `01_tutorial.ipynb` Part 4.
Compares the processed gene lists from different parameter combo of abagen with the gene list generated from previous published works by Anderson et al. Which validated the processed AHBA gene expressions and the imputed cell types that based upon it.
[Anderson, KM, Collins, MA, Chin, R, Ge T, Rosenberg MD, Holmes AJ. (2020). Transcriptional and imaging-genetic association of cortical interneurons, brain function, and schizophrenia risk. Nature Communications, 11, 2889](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7280213/pdf/41467_2020_Article_16710.pdf)

### `01_6_make_abagen_object.R`
- Convert the table into R object for later analyses that run in R.
