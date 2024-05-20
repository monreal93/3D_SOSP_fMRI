# 3D_SOSP_fMRI
This repository houses the sequence and anaylis scripts for the paper "Combining the benefits of 3D acquisitions and spiral readouts for VASO fMRI at UHF"

PULSEQ SEQUENCE GENERATION:
Intial configuration steps:
1) Clone this repository
   git clone ...
2) clone pulseq toolbox into ./tools/pulseq/
   git clone https://github.com/pulseq/pulseq.git
3) clone PNS prediction into ./tools/pns_prediction/
   git clone https://github.com/filip-szczepankiewicz/safe_pns_prediction.git

- Check system limits, if yours is not included, add an extra option for the params.gen.field_strength field and include the system limits in function ./pulseq/functions/prepare_system_limits.m current systems limits:
- 7T Magnetom Plus with SC72 gradients
- 9.4T with AC84-II gradients
- NextGen 7T with impulse gradient

- If you want to check for PNS, inlcude your gradient .asc files in the folder ./tools/pns_prediction/gradient_files and use the flag params.gen.pns_check

To generate a sequence run the file ./pulseq/create_pulseq.m , different parameters are to be set in the section "Define parameters":
- folder_name = "folder_name"  -> Folder where the pulseq and trajectory files will be saved
- seq_name = "seq_name" -> Sequence name (any name you want to use)
- params.gen.seq = "sequence_type" At the moment the well tested sequences are 1-VASO and 4-BOLD

After the script runs several files will be written into the folder ./acq

DATA RECONSTRUCTION:
Initial configuration steps:
1) If Julia is not installed in your system, follow the instructions to install it: https://julialang.org/downloads/
2) 
  
