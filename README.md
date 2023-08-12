# human-cellular-func-con

Analysis code for [The cellular underpinning of human cortical functional connectome (under review)](https://www.biorxiv.org/content/10.1101/2023.07.05.547828v1.abstract)

## Part A: processing and aligning AHBA bulk tissue samples
### `00_aggregate_ahba.R`
- Extract and aggregate necessary information from tables downloaded from [AHBA](http://human.brain-map.org/)
- Input: downloaded AHBA raw data (specifically, the `top_level_categories.csv`, and `SampleAnnot.csv`, `Ontology.csv` for each donor)
- Output: `ahba_sampleInfo.csv`
