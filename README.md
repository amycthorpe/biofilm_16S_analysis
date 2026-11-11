# River biofilm bacterial community analysis

R scripts to process and analyse bacterial 16S rRNA gene sequencing data from river biofilms.  

### Initial processing
Raw reads have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB90117
- Primers are removed from raw sequencing reads with [cutadapt](https://cutadapt.readthedocs.io/en/stable/) using DADA2/cutadapt_script.R
- Reads are processed according to the [DADA2 workflow for big data](https://benjjneb.github.io/dada2/bigdata.html) using DADA2/dada2_script.R
- Taxonomy is assigned according to the [SILVA v138.2](https://www.arb-silva.de/) reference database using DADA2/merging_and_taxonomy_script.R
- Rarefaction and filtering is performed using DADA2/rarefaction_filtering_script.R

### Analysis and visualisation
Outputs (available in the [Zenodo archive](https://zenodo.org/records/18668846)) can then be analysed with the provided scripts to:
   * Assess the diversity, community structure, and composition of the core vs non-core communities using analysis_scripts/community_structure_script.R
   * Assign predicted functions using [FAPROTAX](https://pages.uoregon.edu/slouca/LoucaLab/archive/FAPROTAX/lib/php/index.php) and analysis_scripts/functions_script.R
   * Run ordination and correlation analyses to identify environmental drivers using analysis_scripts/ordination_script.R and analysis_scripts/correlations_script.R
   * Construct co-occurrence networks using analysis_scripts/module_network_script.R
   * Identify environmental thresholds and indicator taxa using analysis_scripts/titan_run.R, visualised using analysis_scripts/titan_plots.R

### Citation
Thorpe et al. (2026) River Microbiomes as Sentinels of National-Scale Freshwater Ecosystems. Global Change Biology 32, no. 3: e70809. https://doi.org/10.1111/gcb.70809.
