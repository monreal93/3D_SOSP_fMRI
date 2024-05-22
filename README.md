# 3D_SOSP_fMRI
This repository houses the sequence and analysis scripts for the paper **Combining the benefits of 3D acquisitions and spiral readouts for VASO fMRI at UHF**

# Sequence (pulseq):
## Intial configuration steps:
1) Clone this repository
   ```shell
   git clone https://github.com/monreal93/3D_SOSP_fMRI.git
   ```
2) clone pulseq toolbox into ./tools/pulseq/
   ```shell
   git clone https://github.com/pulseq/pulseq.git
   ```
3) clone PNS prediction into ./tools/pns_prediction/
   ```shell
   git clone https://github.com/filip-szczepankiewicz/safe_pns_prediction.git
   ```
4) The main file to generate pulseq sequences is **./pulseq/create_pulseq.m**
5) Check system limits, if yours is not included, add an extra option for the params.gen.field_strength field and include the system limits in function ./pulseq/functions/prepare_system_limits.m current systems limits:
   - SIEMENS 7T Magnetom Plus with **SC72** gradients
   - SIEMENS 9.4T with **AC84-II** gradients
   - SIEMENS NextGen 7T with **impulse** gradient

5) If you want to check for PNS, inlcude your gradient .asc files in the folder ***./tools/pns_prediction/gradient_files*** and use the flag ***params.gen.pns_check=1***

6) To generate a sequence run the file ***./pulseq/create_pulseq.m*** , different parameters are to be set in the section "Define parameters":
   - folder_name = "folder_name"  -> Folder where the pulseq and trajectory files will be saved
   - seq_name = "seq_name" -> Sequence name (any name you want to use)
   - params.gen.seq = "sequence_type" At the moment the well tested sequences are 1-VASO and 4-BOLD

7) After the script runs a folder **./data/"folder_name"** will be created, with the following subfolders:
   - acq ->
   - analysis ->
   - ismrmd ->
   - raw ->
   - recon ->
   - sim ->
   - tmp ->
9) Several files will be written into the folder **./data/"folder_name"/acq** , this files will be used to run the sequence at the scanner and to reconstruct the data:
   - "seq_name".seq -> pulseq sequence file
   - "seq_name"_params.mat -> directory with sequence parameters
   - "seq_name"_ks_traj.mat -> nominal trajectory

9) Once you acquire the data, save it in **./data/"folder_name"/raw/twix/**

# RECONSTRUCTION (MRIReco.jl):
## Pre-requisits:
1) [ISMRMD](https://github.com/ismrmrd/ismrmrd)
2) Vendor raw data format to ISMRMRD converter
3) VSCode with remote explorer and julia extensions
   
## Initial configuration steps:
1) Get Julia docker image docker pull julia:1.8.5 and create a container
   ```shell
   docker run -t -d -P -v "path_to_folder":/usr/share/sosp_vaso/ --name julia_mri_recon julia:1.8.5
   docker start julia_mri_recon
   ```
3) I advise using VS code, follow the instructions to attach to a running docker container from [VSCode](https://code.visualstudio.com/docs/devcontainers/attach-container)
4) Install [julia](https://code.visualstudio.com/docs/languages/julia) extension in container 
5) Once attached to contaier, add folder ***/usr/share/sosp_vaso*** to workspace and open file ***./recon/reconstructions_julia_mrd.jl***
6) Open the julia REPL (CTRL+SHIFT+P) `Start julia REPL`
7) Activate julia enviroment in recon folder:
   ```julia
   using Pkg
   Pkg.activate("./recon/")
   ```
9) Instantiate (download and install packages):
   ```julia
   Pkg.instantiate()
   ```

## Reconstruct the sample data
1) Download sample raw data from LINK and save it into **./data/**
2) Convert vendor spiral and GRE-ME (for sensitivity and off-resonance maps) raw data to ISMRMRD format and save it into **./data/sample/ismrmd/3d/**
   ```console
   siemens_to_ismrmrd -f ./data/sample/raw/twix/sample.dat -o ./data/sample/ismrmd/3d/sample.h5
   siemens_to_ismrmrd -f ./data/sample/raw/twix/b0_sample_fieldmap.dat -o ./data/sample/ismrmd/3d/b0_sample_fieldmap.h5
   ```
3) ``` 
  
