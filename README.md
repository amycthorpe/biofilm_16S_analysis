# River biofilm bacterial community analysis

R scripts to process and analyse bacterial 16S rRNA gene sequencing data from river biofilms.  

### Initial processing
Raw sequencing reads are first processed using DADA2 according to XXX_script.R
- Raw reads have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB90117

### Analysis and visualisation
Outputs (available in the Zenodo archive) can then be analysed with the provided R scripts to:
   * Assess the diversity, community structure, and composition of the core community using community_structure_script.R
   * Run ordination and correlation analyses to identify environmental drivers using ordination_script.R and correlation_script.R
   * Construct co-occurrence networks using module_network_script.R
   * Identify environmental thresholds and indicator taxa using titan_run.R, visualised using titan_plots.R

### Citation
https://www.biorxiv.org/content/XXX
