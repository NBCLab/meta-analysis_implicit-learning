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

if [ 1 -eq 0 ]; then
  #make figures for rsfc and macm maps
  rois=$(ls $project_directory/derivatives/rois/)

  for roi in $rois; do
    roi_prefix=${roi%.nii.gz*}
    echo $roi_prefix

    x=${roi_prefix%%_*}
    y=${roi_prefix%_*}
    y=${y#*_}
    z=${roi_prefix##*_}

    #rsfc
    if [ ! -d $project_directory/derivatives/figures/rsfc ]; then
      mkdir -p $project_directory/derivatives/figures/rsfc
    fi
    python3 /Users/miriedel/Desktop/GitHub/surflay/make_slice_figures.py \
      -f $project_directory/derivatives/rsfc/"$x"_"$y"_"$z"_tstat1_thr001.nii.gz \
      --roi_file $project_directory/derivatives/rois/"$x"_"$y"_"$z".nii.gz \
      -x ${x%%.*} \
      -y ${y%%.*} \
      -z ${z%%.*} \
      --outdir $project_directory/derivatives/figures/rsfc \
      --colormap hsv

    #macm
    if [ ! -d $project_directory/derivatives/figures/macm ]; then
      mkdir -p $project_directory/derivatives/figures/macm
    fi
    python3 /Users/miriedel/Desktop/GitHub/surflay/make_slice_figures.py \
      -f $project_directory/derivatives/macm/"$x"_"$y"_"$z"/"$x"_"$y"_"$z"_z_corr-FWE-voxel_thresh-0.001.nii.gz \
      --roi_file $project_directory/derivatives/rois/"$x"_"$y"_"$z".nii.gz \
      -x ${x%%.*} \
      -y ${y%%.*} \
      -z ${z%%.*} \
      --outdir $project_directory/derivatives/figures/macm \
      --colormap hsv

    #conjunction connectivity
    if [ ! -d $project_directory/derivatives/figures/connectivity ]; then
      mkdir -p $project_directory/derivatives/figures/connectivity
    fi
    python3 /Users/miriedel/Desktop/GitHub/surflay/make_slice_figures.py \
      -f $project_directory/derivatives/connectivity/"$x"_"$y"_"$z"/conjunction.nii.gz \
      --roi_file $project_directory/derivatives/rois/"$x"_"$y"_"$z".nii.gz \
      -x ${x%%.*} \
      -y ${y%%.*} \
      -z ${z%%.*} \
      --outdir $project_directory/derivatives/figures/connectivity \
      --colormap hsv

    done
fi
