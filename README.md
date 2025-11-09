
# JUMPt Version 1.1.1

> **✅ Latest:** This repository is for **JUMPt Version 1.1.1** (deterministic & Python included).  
> Looking for the previous release? **JUMPt Version 1.1.0** → **[Go to the older repository](https://github.com/abhijitju06/JUMPt-Version-1.1.0)**.

# JUMPt Version 1.1.1
## Contents <br>
<div align="justify"> 
• Introduction <br>
• Release notes (Version 1.1.1) <br>
• Software and Hardware Requirements <br>
• Installation <br>
• Input File Preparation <br>
• Update the parameter file <br>
• Run the JUMPt program (Demo data sets are provided) <br> 
• Output file <br> 
• Authors <br>
• Maintainers <br>
• Acknowledgments <br>
• References <br>
</div>

## Introduction <br>
<div align="justify"> 
JUMPt (Jumbo Mass Spectrometry-based Proteomics-turnover) software determines protein turnover rates in pulse SILAC-labeled animals using mass spectrometry (MS) data. JUMPt uses a recycling-aware differential equation model that jointly fits the dynamics of unlabeled free Lys and protein-bound Lys to estimate accurate, per-protein half-lives (“Corrected half-life”). The program can also compute apparent half-lives via exponential fits for quick screening. <br>
Version 1.1.1 focuses on <b>deterministic, reproducible results</b> and introduces a <b>robust, fast Python implementation</b> alongside MATLAB. However, the current version supports only *in vivo* data where pulse time is in days. It does not work for *in vitro* data where pulse time is in hours.
</div>

## Release notes (Version 1.1.1) <br>
<div align="justify"> 
In this version: <br>
1. <b>Deterministic initialization</b> — fixed random seeds and fixed start-point sets remove run-to-run variability; identical inputs now yield identical results. <br>
2. <b>Robust optimization path</b> — bounded global fit in γ-space followed by refine in log-γ space; confidence intervals are computed in log-γ and mapped to half-life. <br>
3. <b>Physics & weighting controls</b> — parameters expose the “physics filter” (clip/drop/disable) and <b>free-Lys weighting</b> modes (e.g., MATLAB-like <code>wF = M/5</code>) to replicate legacy behavior or enforce stricter physical constraints. <br>
4. <b>Apparent T50 fix</b> — positivity-constrained exponential fits eliminate negative/invalid apparent half-lives. <br>
5. <b>I/O hardening</b> — more robust Excel parsing (auto-detects time row and free-Lys row), version-stamped outputs, and a full parameter provenance sheet saved with results. <br>
6. <b>Python implementation</b> — a fast, deterministic CLI mirroring MATLAB Setting-2 behavior using exact LTI propagation (matrix exponential) and the same optimization workflow. <br>
7. <b>Compatibility</b> — supports 3, 4, or 5 time points (including 0 days) and only Setting-2.
</div>

## Software and Hardware Requirements <br>
<div align="justify"> 
The MATLAB implementation runs on Linux, macOS, or Windows with MATLAB R2014 or newer (tested on recent versions). The program benefits from multi-core CPUs and ≥16 GB RAM for large datasets. <br>
<b>MATLAB toolbox needed:</b> <br>
- Global Optimization Toolbox (plus standard base toolboxes). <br><br>
<b>Python (optional alternative):</b> Python ≥ 3.9 with the following packages: <code>numpy</code>, <code>pandas</code>, <code>scipy</code>, <code>openpyxl</code>.
</div>

## Installation <br>
<div align="justify"> 
<b>MATLAB:</b> No installation required. Download/clone this repository to any working directory (e.g., <code>/home/usr/JUMPt</code>). <br><br>
<b>Python (optional):</b> <br>
• Option A (run in place): install dependencies and run the CLI script with a params file. <br>
• Option B (package mode): if a Python package is provided, install with <code>pip install .</code> and use the <code>jumpt</code> CLI. 
</div>

<!-- (Optional figure block — add/update image links if desired)
![Figure1](<ADD_LINK_IF_NEEDED>)
<p align="center">Figure 1</p>
-->

## Input File Preparation <br>
<div align="justify"> 
Testing datasets are provided for different time points. Prepare input data similarly, including: <br>
• pSILAC proteins (mandatory) <br>
• pSILAC free (unbound) Lys (required for Setting-2) <br>
By default, <b>Setting-2</b> is used (free-Lys labeling timecourse + protein labeling timecourses).
</div>

## Update the parameter file <br>
<div align="justify"> 
The program reads a parameter file <code>JUMPt.params</code> (MATLAB & Python share the same keys). Specify: <br>
• Input file name (full path) <br>
• Bin size (<code>bin_size</code>) <br>
• Number of time points <br>
• Purity of SILAC food <br>
• <code>alphaUB</code> (upper bound for recycling coupling α) <br>
• Half-life bounds: <code>hl_tmin_days</code> / <code>hl_tmax_days</code> (proteins), <code>hl_A_min_days</code> / <code>hl_A_max_days</code> (free Lys) <br>
• Physics controls: <code>disable_physics_filter</code> (0/1), <code>physics_filter_mode</code> (<code>clip</code>/<code>drop</code>) <br>
• <code>min_observations_per_protein</code> <br>
• Weighting & IC: <code>free_lys_weight_mode</code> (<code>matlab</code> for wF=M/5; <code>auto</code> for wF=M; or numeric), <code>initial_condition_mode</code> (<code>ones</code>/<code>data</code>) <br>
• Ordering: <code>sort_by_apparent_t50</code> (0/1) <br>
• Reproducibility: <code>random_seed</code>, <code>n_starts</code>, <code>max_nfev_global</code>, <code>max_nfev_refine</code> <br>
</div>

## Run the JUMPt program (Demo data set) <br>
<div align="justify"> 
<b>MATLAB:</b> Launch MATLAB and open the main program file (e.g., <code>Run_Main_File.m</code>). Press “Run” to start. The console displays bin progress and completion. <br><br>
<b>Python (CLI):</b> <br>
• Prepare <code>JUMPt_python.params</code> (same keys as above). <br>
• Run: <code>python jumpt_python.py --params JUMPt_python.params</code> <br>
The Python CLI produces the same outputs (including a <i>parameter_file</i> sheet for provenance) and mirrors MATLAB Setting-2 behavior with deterministic results. <br><br>
<b>Note:</b> Nonlinear ODE fitting is computationally intensive for large proteomes; JUMPt bins proteins (e.g., 10–100 per batch) to manage runtime and memory.
</div>

<!-- (Optional figure blocks — add/update image links if desired)
![Figure2](<ADD_LINK_IF_NEEDED>)
<p align="center">Figure 2</p>

![Figure3](<ADD_LINK_IF_NEEDED>)
<p align="center">Figure 3</p>
-->

## Output file <br>
<div align="justify"> 
Two output Excel files are generated with the prefix <code>results_Corrected_T50</code> and <code>results_Apparent_T50</code> appended to the input file name. They are written in the same folder as the input file. <br>
<b>Corrected T50:</b> per-protein corrected half-lives (days), confidence intervals, residuals, extremes flags, and a <i>parameter_file</i> sheet with all parameters used and the program version (1.1.1). <br>
<b>Apparent T50:</b> apparent half-lives from positivity-constrained exponential fits.
</div>

## Authors <br>
<div align="justify">
<b>Dr. Abhijit Dasgupta</b> and <b>Abhisek Bakshi</b>.
</div>

## Maintainers <br>
<div align="justify">
For bug reports and feature suggestions, please open a GitHub issue on this repository.  
You may also contact <b>Dr. Abhijit Dasgupta</b> at <code>abhijitju06@gmail.com</code>.
</div>

## Acknowledgment <br>
<div align="justify"> 
We acknowledge St. Jude Children's Research Hospital, ALSAC (American Lebanese Syrian Associated Charities), and the National Institute of Health for supporting the development of the JUMP Software Suite.
</div>

## References <br>
<div align="justify"> 
1. Chepyala et al., JUMPt: Comprehensive protein turnover modeling of in vivo pulse SILAC data by ordinary differential equations. Analytical Chemistry, 2021. 93(40): 13495–13504. <br>
2. Wang, X., et al., JUMP: a tag-based database search tool for peptide identification with high sensitivity and accuracy. Molecular & Cellular Proteomics, 2014. 13(12): 3663–3673. <br>
3. Wang, X., et al., JUMPm: A Tool for Large-Scale Identification of Metabolites in Untargeted Metabolomics. Metabolites, 2020. 10(5): 190. <br>
4. Li, Y., et al., JUMPg: an integrative proteogenomics pipeline identifying unannotated proteins in human brain and cancer cells. Journal of Proteome Research, 2016. 15(7): 2309–2320. <br>
5. Tan, H., et al., Integrative proteomics and phosphoproteomics profiling reveals dynamic signaling networks and bioenergetics pathways underlying T cell activation. Immunity, 2017. 46(3): 488–503. <br>
6. Peng, J., et al., Evaluation of multidimensional chromatography coupled with tandem mass spectrometry (LC/LC-MS/MS) for large-scale protein analysis: the yeast proteome. Journal of Proteome Research, 2003. 2(1): 43–50. <br>
7. Niu, M., et al., Extensive peptide fractionation and y1 ion-based interference detection method for enabling accurate quantification by isobaric labeling and mass spectrometry. Analytical Chemistry, 2017. 89(5): 2956–2963.  <br>
8. Li, Wenxue, et al. "Turnover atlas of proteome and phosphoproteome across mouse tissues and brain regions." Cell 188.8 (2025): 2267-2287.  <br>
9. Yarbro, Jay M., et al. "Human and mouse proteomics reveals the shared pathways in Alzheimer’s disease and delayed protein turnover in the amyloidome." Nature communications 16.1 (2025): 1533.
</div>
