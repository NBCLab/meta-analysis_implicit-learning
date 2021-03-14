import pickle
import nimare
import os
import os.path as op
import nibabel as nib
import numpy as np


#define some functions
def thresh_img(logp_img, z_img, p):
    sig_inds = np.where(logp_img.get_fdata() > -np.log(p))
    z_img_data = z_img.get_fdata()
    z_img_thresh_data = np.zeros(z_img.shape)
    z_img_thresh_data[sig_inds] = z_img_data[sig_inds]
    z_img = nib.Nifti1Image(z_img_thresh_data, z_img.affine)
    return z_img

def run_ale(dset, output_dir, prefix):
    #define the ALE algorithm
    ale = nimare.meta.ALE()

    #calculate ale statistics
    ale.fit(dset)

    #run FWE multiple comparisons correction
    corr = nimare.correct.FWECorrector(method="montecarlo", n_iters=10000, voxel_thresh=0.001)
    cres = corr.transform(ale.results)

    #save output maps and dataset
    os.makedirs(output_dir, exist_ok=True)
    cres.save_maps(output_dir=output_dir, prefix=prefix)
    dset.save(op.join(output_dir, prefix + ".pkl"))

    #the FWE corrected map contains a -logp value for each cluster, so we need to manually theshold the map to account for that
    z_img_logp = nib.load(op.join(output_dir, '{prefix}_logp_level-cluster_corr-FWE_method-montecarlo.nii.gz'.format(prefix=prefix)))
    z_img = nib.load(op.join(output_dir, '{prefix}_z.nii.gz'.format(prefix=prefix)))
    z_img_thresh = thresh_img(z_img_logp, z_img, 0.05)
    nib.save(z_img_thresh, op.join(output_dir, '{prefix}_z_corr-FWE_thresh-05.nii.gz'.format(prefix=prefix)))


def run_subtraction(dset1, dset2, output_dir, prefix1, prefix2):
    #define the ALE algorithm
    #alesubtraction = nimare.meta.ALESubtraction(n_iters=10000)

    #calculate ale statistics
    #sres = alesubtraction.fit(dset1, dset2)

    #save output maps and dataset
    #os.makedirs(output_dir, exist_ok=True)
    #sres.save_maps(output_dir=op.join(output_dir, '{}_gt_{}'.format(prefix1, prefix2), prefix='{}_gt_{}'.format(prefix1, prefix2))

    #get contrast analysis results
    z_img = nib.load(op.join(output_dir, '{}_gt_{}'.format(prefix1, prefix2), '{}_gt_{}_z_desc-group1MinusGroup2.nii.gz'.format(prefix1, prefix2)))

    #threshold each direction at p < 0.05
    #get dset1 info
    z_img_group1 = nib.load(op.join(output_dir, prefix1, '{}_z.nii.gz'.format(prefix1)))
    z_img_group1_data = z_img_group1.get_fdata()
    #get dset2 info
    z_img_group2 = nib.load(op.join(output_dir, prefix2, '{}_z.nii.gz'.format(prefix2)))
    z_img_group2_data = z_img_group2.get_fdata()

    #for first direction
    sig_inds = np.where(z_img.get_fdata() > -np.log(0.05))

    zimg1_sub_zimg2 = z_img_group1_data - z_img_group2_data

    z_img_group1_thresh_data = np.zeros(z_img_group1.shape)
    z_img_group1_thresh_data[sig_inds] = zimg1_sub_zimg2[sig_inds]

    z_img_group1_final = nib.Nifti1Image(z_img_group1_thresh_data, z_img_group1.affine)
    nib.save(z_img_group1_final, op.join(output_dir, '{}_gt_{}'.format(prefix1, prefix2), '{}-{}_z_thresh-05.nii.gz'.format(prefix1, prefix2)))

    #for second direction
    sig_inds = np.where(z_img.get_fdata() < np.log(0.05))

    zimg2_sub_zimg1 = z_img_group2_data - z_img_group1_data

    z_img_group2_thresh_data = np.zeros(z_img_group2.shape)
    z_img_group2_thresh_data[sig_inds] = zimg2_sub_zimg1[sig_inds]

    z_img_group2_final = nib.Nifti1Image(z_img_group2_thresh_data, z_img_group2.affine)
    nib.save(z_img_group2_final, op.join(output_dir, '{}_gt_{}'.format(prefix1, prefix2), '{}-{}_z_thresh-05.nii.gz'.format(prefix2, prefix1)))



#main part of script
project_directory = '/Users/miriedel/Desktop/GitHub/meta-analysis_implicit-learning'

# load the nimare dataset
with open(op.join(project_directory, 'code', 'nimare_dataset.pkl'), 'rb') as fo:
    dset = pickle.load(fo)

# only select those contrasts that should be included
dset_include_idx = dset.get_studies_by_label(labels='Include Contrast? (0,1)')
dset_include = dset.slice(dset_include_idx)

#select the grammatical contrasts
grammatical_idx = dset_include.annotations['id'][dset_include.annotations['Meta-Analysis Group'] == 'Grammatical']
grammatical_dset = dset_include.slice(grammatical_idx)

#define output directory and run ale function defined above
output_dir = op.join(project_directory, 'derivatives', 'grammatical')
#run_ale(grammatical_dset, output_dir, 'grammatical')

#select the ungrammatical contrasts
ungrammatical_idx = dset_include.annotations['id'][dset_include.annotations['Meta-Analysis Group'] == 'Ungrammatical']
ungrammatical_dset = dset_include.slice(ungrammatical_idx)

#define output directory and run ale function defined above
output_dir = op.join(project_directory, 'derivatives', 'ungrammatical')
#run_ale(ungrammatical_dset, output_dir, 'ungrammatical')

#run subtraction analysis
output_dir = op.join(project_directory, 'derivatives')
run_subtraction(grammatical_dset, ungrammatical_dset, output_dir, 'grammatical', 'ungrammatical')
