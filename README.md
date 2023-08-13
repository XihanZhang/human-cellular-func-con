# human-cellular-func-con

Analysis code for [The cellular underpinning of human cortical functional connectome (under review)](https://www.biorxiv.org/content/10.1101/2023.07.05.547828v1.abstract)

## Part A: processing and aligning AHBA bulk tissue samples
### `00_aggregate_ahba.R`
- Extract and aggregate necessary information from tables downloaded from [AHBA](http://human.brain-map.org/)
- Input: downloaded AHBA raw data (specifically, the `top_level_categories.csv`, and `SampleAnnot.csv`, `Ontology.csv` for each donor)
- Output: `ahba_sampleInfo.csv`
### `01_tutorial.ipynb` Part 1.
- Using abagen to correct the old MNI coordinates
- Input: `ahba_sampleInfo.csv`
- Output: `ahba_sampleInfo_reannot.csv`
### `01_1_project_freesurfer.bash`
- Project fsLR parcellations to individual AHBA freesurfer space
- Input and output files are in the `mri`, `surf`, `label` folders under each donor dirs that sit wherever the abagen download dir sit, within `/abated-data/freesurfer`.

(e.g. abagen by default download the files to `/Users/zhangxihan/abated-data/freesurfer/donor9861/mri`, but I moved the entire freesurfer folder to `/gpfs/milgram/project/holmes/xz555/gradient_shift/data/ahba_fs`)
