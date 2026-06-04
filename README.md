# CyanoPeptideMatching (CPM) Python Version

This project provides a comprehensive analytical workflow for analyzing LC-MS/MS mass spectrometry fies (mzML format) for the rapid identification and classification of cyanobacterial secondary metabolites. CPM utilizes diagnostic-ion-guided MS/MS querying for specific cyanopeptide classes, precursor-level MS1 consolidation, adduct-feature relationship mapping, and database-assisted matching to CyanoMetDB. The workflow employs the MassQL query language to perform flexible, reproducible searches across MS2 spectra based on class-specific diagnostic ions, mass tolerances, and retention time windows. CPM supports relative trend analysis by extracting MS1 chromatographic signal directly from the user-provided LC–MS/MS datasets within the  workflow. Extracted ion chromatograms (EICs) are generated for matched precursor features, and area under the curve (AUC) values are computed to provide relative abundance estimates. A user-defined reference compound enables normalization of MS1 AUC values across samples, facilitating relative comparison of feature abundances while maintaining a streamlined analysis pipeline.

*   Technologies used 
        *   Proteowizard (converting raw mass spectrometry files to .mzML)
        *   Python (Python 3.13)

#-----------------------------------------------------------------------------------------------
# Getting Started 
## Linux/macOS
Open terminal
# Install micromamba once
cd ~ \
curl -Ls https://micro.mamba.pm/install.sh | bash \
source ~/.bashrc

# Create environment
micromamba create -n cpm python=3.13 pip -c conda-forge -y

# Activate it
micromamba activate cpm

# Check Python version (should be 3.13.x)
python --version

# Upgrade install tools
pip install --upgrade pip setuptools wheel

# Install CPM package
pip install --no-cache-dir --force-reinstall \
git+https://github.com/sheffera01/CyanopeptideMatchingPythonVersion.git@package

# For Future Logins
source ~/.bashrc \
micromamba activate cpm

# For newest version pulling from GitHub:
pip install --no-cache-dir --force-reinstall \
git+https://github.com/sheffera01/CyanopeptideMatchingPythonVersion.git@package

# To run: 
cpm \
  --class-tag MC \ \
  --files /path_to_file/filename.mzML \ \
          /path_to_file/filename2.mzML \ \
  --metadata /path_to__metadata/metadata.csv \ \
  --output-root /path_to_result_output/results \ \
  --blank-filter \ \
  --batch-correct \
  
Use cpm --help for usage informnation

# Options:
--class tag: MC, AP, AB, AR, MG, or ALL \
--blank-filter: optional. Remove if not needed \
--batch-correct: optional. Remove if not needed. 

#--------------------------------------------------------------------------------------------------
## PC
Open PowerShell
# Install micromamba once
1. cd $HOME
2. Invoke-WebRequest -Uri https://micro.mamba.pm/api/micromamba/win-64/latest -OutFile micromamba.tar.bz2
3. tar xf micromamba.tar.bz2
4. .\Library\bin\micromamba.exe shell init -s powershell -r "$HOME\micromamba" 

Close PowerShell and reopen it. 
# Create environment
micromamba create -n cpm python=3.13 pip -c conda-forge -y

# Activate it
micromamba activate cpm

# Check Python version (should be 3.13.x)
python --version

# Upgrade install tools
pip install --upgrade pip setuptools wheel

# Install CPM package
pip install --no-cache-dir --force-reinstall "git+https://github.com/sheffera01/CyanopeptideMatchingPythonVersion.git@package"

#-----------------------------------------------------------------------------------------------
# For Future Logins
Open PowerShell and run:
micromamba activate cpm

# For newest version pulling from GitHub:
pip install --no-cache-dir --force-reinstall "git+https://github.com/sheffera01/CyanopeptideMatchingPythonVersion.git@package"

#-----------------------------------------------------------------------------------------------
# To run: 
cpm ` --class-tag MC ` --files ` "C:\path_to_file\filename.mzML" ` "C:\path_to_file\filename2.mzML" ` --metadata "C:\path_to_metadata\metadata.csv" ` --output-root "C:\path_to_output\results" ` --blank-filter ` --batch-correct \

Use cpm --help for usage informnation

# Options:
--class tag: MC, AP, AB, AR, MG, or ALL \
--blank-filter: optional. Remove if not needed \
--batch-correct: optional. Remove if not needed. \
#-----------------------------------------------------------------------------------------------
### Prerequisites

#%pip install pandas (2.3.3) matplotlib (3.10.8) numpy (2.2.6) seaborn (0.13.2) networkx (3.4.2) massql (2025.12.10) openpyxl (3.1.5) pyteomics (4.7.5) lxml (6.0.2) \
need to be running also in massql_env (Python 3.13.X)

#-----------------------------------------------------------------------------------------------
## Areas for users to edit specific parameters
Adding to ion search lists (cpm.py --> edit Ion Lists AND Ion label dictionary) \
Blank Ratio Threshold (cpm.py --> blank_ratio_threshold) \
Tolerance matching to CyanoMetDB (cyanometdb_match.py --> edit tol_da: float = __; cpm.py --> edit tol_da: float = 0.05,) \
Number of Diagnostic Product Ions Needed for matches (cyanometdb_match.py --> "n_diagnostic"] >= 2) \
Adduct additions (adduct_finder.py --> add into default_deltas) \
Adjusting merged summary (summary_builder.py --> merge_tol_mz and merge_tol_rt edit) \
MS1 AUC (summary_builder.py --> tol_mz and rt_pad edit) \
MS/MS diagnostic ion tolerance fragment (massql_utils.py --> tol_mz) \
#-----------------------------------------------------------------------------------------------
## Acknowledgements 

*    We are grateful for support from the National Institute of Environmental Health Sciences (NIEHS) of the NIH under award numbers 5P01ES028939-02 and R21ES033758 (M.J.B.) and the National Science Foundation (NSF) under award number OCE-1840715, T32 GM140223 Pharmacological Sciences Training Program (S. L. H.),  the National Institute of Health (NIH) F31 1F31ES036421-01 (L.N.H.), the National Institute of Health (NIH) F31 1F31AI186432-01 (K.L.L.) We thank the United States Geological Survey (USGS) and NOAA/GLERL for providing access to environmental metabolomics datasets used in this study. 
