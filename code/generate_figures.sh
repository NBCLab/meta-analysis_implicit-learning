project_directory='/Users/miriedel/Desktop/GitHub/meta-analysis_implicit-learning/'

#make figures for grammatical, ungrammatical, conjunction, contrast analyses
#grammatical
python3 /Users/miriedel/Desktop/GitHub/surflay/make_figures.py \
  -f $project_directory/derivatives/grammatical/grammatical_z_corr-FWE-cluster_thresh-05.nii.gz \
  -o $project_directory/derivatives/figures/ \
  -z 38 30 14 -2 \
  --colormap Reds \
  --opacity 0.75

#ungrammatical
python3 /Users/miriedel/Desktop/GitHub/surflay/make_figures.py \
  -f $project_directory/derivatives/ungrammatical/ungrammatical_z_corr-FWE-cluster_thresh-05.nii.gz \
  -o $project_directory/derivatives/figures/ \
  -z 42 32 18 -4 \
  --colormap Blues \
  --opacity 0.75

#conjunction
python3 /Users/miriedel/Desktop/GitHub/surflay/make_figures.py \
  -f $project_directory/derivatives/conjunction/conjunction.nii.gz \
  -o $project_directory/derivatives/figures/ \
  -z 26 16 0 -6 \
  --colormap Oranges \
  --opacity 0.75

#pooled
python3 /Users/miriedel/Desktop/GitHub/surflay/make_figures.py \
  -f $project_directory/derivatives/pooled/pooled_z_corr-FWE-cluster_thresh-05.nii.gz \
  -o $project_directory/derivatives/figures/ \
  -z 48 38 20 -2 \
  --colormap Greens \
  --opacity 0.75

#make figures for rsfc and macm maps
rois=$(ls $project_directory/derivatives/pooled_analysis/rois/)

for roi in $rois; do
  roi_prefix=${roi%.nii.gz*}
  echo $roi_prefix

  x=${roi_prefix%%_*}
  y=${roi_prefix%_*}
  y=${y#*_}
  z=${roi_prefix##*_}

  #rsfc
  if [ ! -d $project_directory/derivatives/figures/pooled_analysis/rsfc ]; then
    mkdir -p $project_directory/derivatives/figures/pooled_analysis/rsfc
  fi
  python3 /Users/miriedel/Desktop/GitHub/surflay/make_slice_figures.py \
    -f $project_directory/derivatives/pooled_analysis/rsfc/"$x"_"$y"_"$z"_tstat1_thr001.nii.gz \
    --roi_file $project_directory/derivatives/pooled_analysis/rois/"$x"_"$y"_"$z".nii.gz \
    -x ${x%%.*} \
    -y ${y%%.*} \
    -z ${z%%.*} \
    --outdir $project_directory/derivatives/figures/pooled_analysis/rsfc \
    --colormap hsv

  #macm
  if [ ! -d $project_directory/derivatives/figures/pooled_analysis/macm ]; then
    mkdir -p $project_directory/derivatives/figures/pooled_analysis/macm
  fi
  python3 /Users/miriedel/Desktop/GitHub/surflay/make_slice_figures.py \
    -f $project_directory/derivatives/pooled_analysis/macm/"$x"_"$y"_"$z"/"$x"_"$y"_"$z"_z_corr-FWE-voxel_thresh-0.001.nii.gz \
    --roi_file $project_directory/derivatives/pooled_analysis/rois/"$x"_"$y"_"$z".nii.gz \
    -x ${x%%.*} \
    -y ${y%%.*} \
    -z ${z%%.*} \
    --outdir $project_directory/derivatives/figures/pooled_analysis/macm \
    --colormap hsv

  #conjunction connectivity
  if [ ! -d $project_directory/derivatives/figures/pooled_analysis/connectivity ]; then
    mkdir -p $project_directory/derivatives/figures/pooled_analysis/connectivity
  fi
  python3 /Users/miriedel/Desktop/GitHub/surflay/make_slice_figures.py \
    -f $project_directory/derivatives/pooled_analysis/connectivity/"$x"_"$y"_"$z"/conjunction.nii.gz \
    --roi_file $project_directory/derivatives/pooled_analysis/rois/"$x"_"$y"_"$z".nii.gz \
    -x ${x%%.*} \
    -y ${y%%.*} \
    -z ${z%%.*} \
    --outdir $project_directory/derivatives/figures/pooled_analysis/connectivity \
    --colormap hsv

  done
