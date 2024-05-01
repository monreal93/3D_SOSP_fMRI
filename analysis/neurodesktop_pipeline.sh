### ToDO:

# Neurodesktop tools to use:
ml afni
ml fsl
ml laynii

cd /neurodesktop-storage/5T4/Alejandro/sosp_vaso/data

folder="05192023_sv_paper"
scan="cv_01"
r_a_tr=7     # rest/activity TRs paper=sv(8), cv(7)
tr=1.81      # volume TR paper=sv(1.66), cv(1.81)

# Spiral reconstruction options
traj="_nom" && cs="_cs" && b0="_b0" && k0="_k0" && rDORK="_rDORK"

cd ${folder}
mkdir ./analysis/${scan}
cd ./analysis/${scan}

if [ "${scan:0:2}" = "sv" ]; then
    echo "Spiral VASO .."
    # Set path for reconstruction
    v_file=../../recon/${scan}_v${traj}${cs}${b0}${k0}${rDORK}.nii
    b_file=../../recon/${scan}_b${traj}${cs}${b0}${k0}${rDORK}.nii
    gre1=../../tmp/${scan}_1ech.nii
elif [ "${scan:0:2}" = "cv" ]; then
    echo "Cartesian VASO .."
    file=../../recon/${scan}_bv_epi.nii
    v_file=../../recon/${scan}_v_epi.nii
    b_file=../../recon/${scan}_b_epi.nii
    # Split VASO and BOLD
    3dTcat -prefix ${b_file} ${file}'[0..$(2)]' -overwrite
    3dTcat -prefix ${v_file} ${file}'[1..$(2)]' -overwrite
fi

# ToDo: Get vol from recon nifti and write tr in nifti.. (when merging after recon)
# vol=160  # sv_01/02 =140, sv_03=220
# tr=1.66  # sv_01/02 = 1.87 , sv_03=0.93
# r_a_tr=8  # sv_01/02 = 7 , sv_03=11

vol=$(3dinfo -nv ${v_file})
# tr=$(3dinfo -tr ${v_file})
# tr=$(echo $tr*1000 | bc -l)

blocks=$(echo $vol/$r_a_tr/2 | bc -l)
blocks=$(echo ${blocks%.*})
block_dur=$(echo $tr*$blocks*2 | bc -l)
block_dur=$(echo ${block_dur%.*})
block_upsample=$(echo $blocks*2 | bc -l)
block_trs=$(echo $r_a_tr*2 | bc -l)
block_trs=$(echo ${block_trs%.*})

# ToDo: Find an easy way to realign Anatomy as fMRI...
# Center and deoblique datasets
3dinfo -obliquity ${b_file} >> obliquity.txt 
3drefit -xorigin cen -yorigin cen -zorigin cen ${v_file}
3drefit -xorigin cen -yorigin cen -zorigin cen ${b_file}
3drefit -xorigin cen -yorigin cen -zorigin cen ${gre1}
3dWarp -deoblique ${v_file}
3dWarp -deoblique ${b_file}
3dWarp -deoblique ${gre1}

# 1) Creating Mask
#--- with FSL
# bet $gre1 ./tmp -Z -f 0.4 -g 0 -n -m
# mv tmp_mask.nii.gz mask.nii.gz
# gzip -d mask.nii.gz
#--- with AFNI
if [ "${scan:0:2}" = "sv" ]; then
    3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${gre1}
else
    3dAutomask -prefix mask.nii -peels 3 -dilate 2 ${v_file}
fi
3drefit -xdel $(3dinfo -adi ${b_file}) mask.nii
3drefit -ydel $(3dinfo -adj ${b_file}) mask.nii
3drefit -zdel $(3dinfo -adk ${b_file}) mask.nii

# Replacing non-steady state volumes..
# ToDo: How many volumes should I replace?
3dTcat -prefix ${b_file} ${b_file}'[4..7]' ${b_file}'[4..$]' -overwrite
3dTcat -prefix ${v_file} ${v_file}'[4..7]' ${v_file}'[4..$]' -overwrite

# ) Motion correction
# ToDo: Check if mask can be used...
# ToDo: Use same function as Renzo's new scripts... 3dAllineate
echo "Motion correction... BOLD..."
# --- 3dvolreg
3dvolreg -base 4 -heptic -zpad 1 -overwrite -prefix ./${scan}_b_mc.nii -1Dfile motion_b.1D ${b_file}
# --- 3dAllineate
# 3dTstat -mean -overwrite -prefix ./n_ref.nii ${b_file}'[0..3]'  # create reference
# 3dAllineate -1Dmatrix_save  matrix.aff12.1D -1Dparam_save   param.aff12 -cost lpa \
#     -prefix ./${scan}_b_mc.nii -base ${output_dir} ./n_ref.nii -source ${b_file} \
#     -weight ./mask.nii -warp shift_rotate -final wsinc5

echo "Motion correction... VASO..."
# --- 3dvolreg
3dvolreg -base 4 -heptic -zpad 1 -overwrite -prefix ./${scan}_v_mc.nii -1Dfile motion_v.1D ${v_file}
# --- 3dAllineate
# 3dTstat -mean -overwrite -prefix ./n_ref.nii ${v_file}'[0..3]'  # create reference
# 3dAllineate -1Dmatrix_save  matrix.aff12.1D -1Dparam_save   param.aff12 -cost lpa \
#     -prefix ./${scan}_v_mc.nii -base ${output_dir} ./n_ref.nii -source ${v_file} \
#     -weight ./mask.nii -warp shift_rotate -final wsinc5

#### ) Temporal upsampling
3dUpsample -overwrite  -datum short -prefix ./${scan}_b_mc_ups.nii -n 2 -input ./${scan}_b_mc.nii
3dUpsample -overwrite  -datum short -prefix ./${scan}_v_mc_ups.nii -n 2 -input ./${scan}_v_mc.nii
# Updating TR times
3drefit -TR $tr ./${scan}_b_mc_ups.nii
3drefit -TR $tr ./${scan}_v_mc_ups.nii

#### ) Temporal filtering
echo "Temporal High-pass filter..."
bptf=$(echo $block_dur/$tr/2 | bc -l)
bptf=$(echo ${bptf%.*})
fslmaths ./${scan}_b_mc_ups.nii -Tmean tempmean
fslmaths ./${scan}_b_mc_ups.nii -bptf $bptf -1 -add tempmean ./${scan}_b_mc_ups_hpf.nii
fslchfiletype NIFTI ${scan}_b_mc_ups_hpf.nii
fslmaths ./${scan}_v_mc_ups.nii -Tmean tempmean
fslmaths ./${scan}_v_mc_ups.nii -bptf $bptf -1 -add tempmean ./${scan}_v_mc_ups_hpf.nii
fslchfiletype NIFTI ${scan}_v_mc_ups_hpf.nii

#### ) BOLD correction
LN_BOCO -Nulled ./${scan}_v_mc_ups_hpf.nii -BOLD ./${scan}_b_mc_ups_hpf.nii -trialBOCO $block_trs

# #### ) Calculating T1
# echo "calculating T1 ..."
# # 3dcalc -a sample.nii -expr '0' -prefix combined.nii
# # 3dcalc -prefix combined.nii -a combined.nii'[0..$(2)]' -b ${b_file} -expr 'a+b' -overwrite

# NumVol=`3dinfo -nv ./${scan}_b_mc_ups_hpf.nii`
# # 3dcalc -a ./${scan}_b_mc_ups_hpf.nii'[3..'`expr $NumVol - 2`']' -b  ./${scan}_v_mc_ups_hpf.nii'[3..'`expr $NumVol - 2`']' -expr 'a+b' -prefix combined.nii -overwrite
# # ToDo: Optimize combination of vaso and bold
# 3dTcat -overwrite ./${scan}_b_mc_ups_hpf.nii'[0]' -prefix combined.nii
# 3dTcat -overwrite ./combined.nii ./${scan}_v_mc_ups_hpf.nii'[0]' -prefix combined.nii
# for (( i = 1; i <= ${NumVol}; i++ ))
# do
#     3dTcat -overwrite ./combined.nii ./${scan}_b_mc_ups_hpf.nii"[$i]" -prefix combined.nii
#     3dTcat -overwrite ./combined.nii ./${scan}_v_mc_ups_hpf.nii"[$i]" -prefix combined.nii
# done
# 3dTstat -cvarinv -prefix T1_weighted.nii -overwrite combined.nii 
# rm combined.nii

#### ) Quality metrics
echo "calculating Mean and tSNR maps ..."
3dTstat -mean -prefix mean_b.nii ./${scan}_b_mc_ups_hpf.nii'[1..$]' -overwrite
3dTstat -mean -prefix mean_v.nii ./VASO_LN.nii'[1..$]' -overwrite
3dTstat  -overwrite -cvarinv  -prefix tSNR_b.nii ./${scan}_b_mc_ups_hpf.nii'[1..$]'
3dTstat  -overwrite -cvarinv  -prefix tSNR_v.nii VASO_LN.nii'[1..$]'

#### ) Activation maps
block_dur=$(echo $tr*$block_trs | bc -l)
tmp=0
start=1
end=$(echo $blocks-1 | bc -l)
stim_times='1D: 0 '

for (( i=$start; i<=$end; i++))
do
	tmp=$(echo $tmp+$block_dur*2 | bc -l)  # If temporal upsampling
    # tmp=$(echo $tmp+$block_dur | bc -l)
	stim_times=$(echo "$stim_times $tmp ")
done
tmp=$(echo $block_dur | bc -l)
tmp=$(echo ${tmp%%.*})
ublock=$(echo "UBLOCK($tmp,1)")

# Finding rest and activity volumes
tmp=$(echo $block_trs/4 | bc -l)
tmp=$(echo ${tmp%.*})
r1=$(echo $tmp-1 | bc -l)
r1=$(echo ${r1%.*})
r2=$(echo $r1+$tmp | bc -l)
r2=$(echo ${r2%.*})
a1=$(echo $r2+$tmp | bc -l)
a1=$(echo ${a1%.*})
a2=$(echo $a1+$tmp | bc -l)
a2=$(echo ${a2%.*})
tr_av_r=$(echo "[$r1-$r2]")
tr_av_a=$(echo "[$a1-$a2]")

#### VASO based on difference
echo "VASO based on difference..."
3dTstat -mean -prefix VASO_r.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix VASO_a.nii -overwrite  VASO_trialAV_LN.nii"$tr_av_a"
3dcalc -a  VASO_r.nii -b VASO_a.nii -overwrite -expr '(a-b)/a' -prefix delta_VASO.nii

#### VASO GLM
echo "VASO based on GLM..." 
3dDeconvolve -overwrite -jobs 16 -polort a -input VASO_LN.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_VASO.nii \
             -bucket STATS_VASO.nii
             
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a'    -prefix 1_HRF_VASO.nii   -overwrite 
3dcalc -a HRF_VASO.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_VASO.nii   -overwrite 

3dcalc -a STATS_VASO.nii'[0]'  -expr 'a'    -prefix 0_STATS_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_VASO.nii -overwrite
3dcalc -a STATS_VASO.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_VASO.nii -overwrite 
3dcalc -a STATS_VASO.nii'[2]'  -expr 'a'    -prefix 2_STATS_VASO.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii VASO_LN.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_VASO.nii      -expr 'b/a*100' -prefix 1_HRF_percent_VASO.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_VASO.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_VASO.nii

#### BOLD based on difference
echo "BOLD based on difference... " 
3dTstat -mean -prefix BOLD_r.nii -overwrite  BOLD_trialAV_LN.nii"$tr_av_r"
3dTstat -mean -prefix BOLD_a.nii -overwrite  BOLD_trialAV_LN.nii"$tr_av_a"
3dcalc -a  BOLD_r.nii -b BOLD_a.nii -overwrite -expr '(b-a)/a' -prefix delta_BOLD.nii

#### BOLD
echo "BOLD based on GLM..."
3dDeconvolve -overwrite -jobs 16 -polort a -input ./${scan}_b_mc_ups_hpf.nii\
             -num_stimts 1 \
             -TR_times $tr \
             -stim_times 1 "$stim_times" "$ublock" -stim_label 1 Task \
             -tout \
             -x1D MODEL_wm \
             -iresp 1 HRF_BOLD.nii \
             -bucket STATS_BOLD.nii

3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a'    -prefix 1_HRF_BOLD.nii   -overwrite 
3dcalc -a HRF_BOLD.nii'[1]'    -expr 'a*-1'    -prefix 1_HRF_NEG_BOLD.nii   -overwrite 

3dcalc -a STATS_BOLD.nii'[0]'  -expr 'a'    -prefix 0_STATS_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[1]'  -expr '-1*a' -prefix 1_STATS_NEG_BOLD.nii -overwrite
3dcalc -a STATS_BOLD.nii'[2]'  -expr '-1*a' -prefix 2_STATS_NEG_BOLD.nii -overwrite 
3dcalc -a STATS_BOLD.nii'[2]'  -expr 'a'    -prefix 2_STATS_BOLD.nii -overwrite 

3dTstat -mean -overwrite -prefix mean.nii ./${scan}_b_mc_ups_hpf.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_BOLD.nii      -expr 'b/a*100' -prefix 1_HRF_percent_BOLD.nii
3dcalc -a mean.nii -overwrite -b 1_HRF_NEG_BOLD.nii -expr 'b/a*100' -prefix 1_HRF_NEG_percent_BOLD.nii

####### 10) Masking relevant volumes
3dcalc -a T1_weighted.nii -b mask.nii -expr 'a*b' -prefix T1_msk.nii
3dcalc -a mean_v.nii -b mask.nii -expr 'a*b' -prefix mean_v_msk.nii
3dcalc -a mean_b.nii -b mask.nii -expr 'a*b' -prefix mean_b_msk.nii
3dcalc -a 2_STATS_VASO.nii -b mask.nii -expr 'a*b' -prefix VASO_msk.nii
3dcalc -a 2_STATS_NEG_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_msk.nii

3dcalc -a delta_BOLD.nii -b mask.nii -expr 'a*b' -prefix BOLD_delta_msk.nii
3dcalc -a delta_VASO.nii -b mask.nii -expr 'a*b' -prefix VASO_delta_msk.nii

##### ) Cluster activations
# Here, I need to play with the values after -1clip:
# -1clip threshold.. (~1.8), (1.5,1.2,270)
# rmm = cluster connection radius, larger value->remove small clusters
# vmul minimum cluster volume, smaller value->removes small clusters
3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 2 1.4 120 VASO_msk.nii
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 3 1.4 120 BOLD_msk.nii

3dclust -1noneg -overwrite -prefix clustered_VASO.nii -1clip 0.02 1 120 delta_VASO.nii
3dclust -1noneg -overwrite -prefix clustered_BOLD.nii -1clip 0.02 1 120 delta_BOLD.nii

####### ) Mean tSNR and effective tSNR
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./tSNR_b.nii) 
echo -e "BOLD brain mean tSNR: \n $mean_tSNR_b" >> results.txt
mean_tSNR_v=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./tSNR_v.nii) 
echo -e "VASO brain mean tSNR: \n $mean_tSNR_v" >> results.txt

3dcalc -a ./tSNR_b.nii -expr "a/sqrt($tr)" -prefix ./eff_tSNR_b.nii
mean_tSNR_b=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./eff_tSNR_b.nii) 
echo -e "BOLD brain mean effective tSNR: \n $mean_tSNR_b" >> results.txt
3dcalc -a ./tSNR_v.nii -expr "a/sqrt($tr)" -prefix ./eff_tSNR_v.nii
mean_tSNR_v=$(3dBrickStat -mean -non-negative -nonan -mask ./mask.nii ./eff_tSNR_v.nii) 
echo -e "VASO brain mean effective tSNR: \n $mean_tSNR_v" >> results.txt

#######  ) Quantitive fMRI analysis
## 1) Manually align ananotmy to fMRI in ITK-SNAP, and create mask for registration in ANTs
## 2) Run neurodesktop_ants.sh to perform registration of T1 to fMRI space
## 3) Run neurodesktop_freesurfer.sh to perform segmentations...

# # Let's resample the gm and wm mask from freesurfer to fMRI space
# 3dresample -master ./T1_weighted_masked.nii -prefix gm_mask.nii -input gm_mask.nii -overwrite
# 3dresample -master ./T1_weighted_masked.nii -prefix wm_mask.nii -input wm_mask.nii -overwrite

# Let's make them binary mask
3dcalc -a ./"$scan"_gm_msk.nii -expr 'ispositive(a-20)' -prefix "$scan"_gm_msk1.nii -overwrite
3dcalc -a ./"$scan"_wm_msk.nii -expr 'ispositive(a-0.8)' -prefix "$scan"_wm_msk1.nii -overwrite

# Let's mask to a ROI in visual cortex
3dcalc -a "$scan"_gm_msk1.nii -b "$scan"_roi_msk.nii -expr 'a*b' -prefix "$scan"_gm_msk1.nii -overwrite
3dcalc -a "$scan"_wm_msk1.nii -b "$scan"_roi_msk.nii -expr 'a*b' -prefix "$scan"_wm_msk1.nii -overwrite

# Let's split the activations from GM and WM
# ToDo: Not sure if I want to use clustered_BOLD or BOLD_mask
3dcalc -a ./"$scan"_gm_msk1.nii -b ./BOLD_msk.nii -expr 'a*b' -prefix BOLD_gm.nii -overwrite
3dcalc -a ./"$scan"_wm_msk1.nii -b ./BOLD_msk.nii -expr 'a*b' -prefix BOLD_wm.nii -overwrite
3dcalc -a ./"$scan"_gm_msk1.nii -b ./VASO_msk.nii -expr 'a*b' -prefix VASO_gm.nii -overwrite
3dcalc -a ./"$scan"_wm_msk1.nii -b ./VASO_msk.nii -expr 'a*b' -prefix VASO_wm.nii -overwrite

# Let's ommit the first two slices (fold over artifacts)
slices=$(3dinfo -nk "$scan"_gm_msk.nii)
slices=$(($slices-1))
3dZcutup -keep 2 "${slices}" -prefix BOLD_gm.nii BOLD_gm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix BOLD_wm.nii BOLD_wm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix VASO_gm.nii VASO_gm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix VASO_wm.nii VASO_wm.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix "$scan"_gm_msk1.nii "$scan"_gm_msk1.nii -overwrite
3dZcutup -keep 2 "${slices}" -prefix "$scan"_wm_msk1.nii "$scan"_wm_msk1.nii -overwrite

rm BOLD_roc*
rm VASO_roc*
step=0.5
##### WM VASO ROC
echo -e "Threshold \t Mean \t Voxels" >> VASO_roc_wm.txt
3dcalc -overwrite -a "$scan"_wm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_wm_msk1.nii)
echo -e "\n 99   $VASO_roc" >> VASO_roc_wm.txt

thr=0
for (( i=1; i<=10; i++))
do
    3dcalc -overwrite -a VASO_wm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet VASO_wm.nii)
    echo -e "\n $thr   $VASO_roc" >> VASO_roc_wm.txt
    thr=$(echo $thr+0.5 | bc -l)
done

##### GM VASO ROC
echo -e "Threshold \t Mean \t Voxels" >> VASO_roc_gm.txt
3dcalc -overwrite -a "$scan"_gm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_gm_msk1.nii)
echo -e "\n 99   $VASO_roc" >> VASO_roc_gm.txt

thr=0
for (( i=1; i<=10; i++))
do
    3dcalc -overwrite -a VASO_gm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    VASO_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet VASO_gm.nii)
    echo -e "\n $thr   $VASO_roc" >> VASO_roc_gm.txt
    thr=$(echo $thr+0.5 | bc -l)
done

##### WM BOLD ROC
echo -e "Threshold \t Mean \t Voxels" >> BOLD_roc_wm.txt
3dcalc -overwrite -a "$scan"_wm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_wm_msk1.nii)
echo -e "\n 99   $BOLD_roc" >> BOLD_roc_wm.txt

thr=0
for (( i=1; i<=10; i++))
do
    3dcalc -overwrite -a BOLD_wm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet BOLD_wm.nii)
    echo -e "\n $thr   $BOLD_roc" >> BOLD_roc_wm.txt
    thr=$(echo $thr+0.5 | bc -l)
done

##### GM BOLD ROC
echo -e "Threshold \t Mean \t Voxels" >> BOLD_roc_gm.txt
3dcalc -overwrite -a "$scan"_gm_msk1.nii -expr 'step(a-0)' -prefix binarised_output.nii
BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet "$scan"_gm_msk1.nii)
echo -e "\n 99   $BOLD_roc" >> BOLD_roc_gm.txt

thr=0
for (( i=1; i<=10; i++))
do
    3dcalc -overwrite -a BOLD_gm.nii -expr 'step(a-'$thr')' -prefix binarised_output.nii
    BOLD_roc=$(3dROIstats -mask binarised_output.nii -nzvoxels -quiet BOLD_gm.nii)
    echo -e "\n $thr   $BOLD_roc" >> BOLD_roc_gm.txt
    thr=$(echo $thr+0.5 | bc -l)
done