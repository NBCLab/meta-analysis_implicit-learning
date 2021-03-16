import pickle
import nimare
import os
import os.path as op
import nibabel as nib
import numpy as np
import pandas as pd
from atlasreader import atlasreader


def conjunction(img1_fn, img2_fn, output_dir):
    img1 = nib.load(img1_fn)
    img2 = nib.load(img2_fn)
    inds1 = np.ravel_multi_index(np.nonzero(img1.get_fdata()), img1.shape)
    inds2 = np.ravel_multi_index(np.nonzero(img2.get_fdata()), img2.shape)
    inds12 = np.unravel_index(np.intersect1d(inds1, inds2), img1.shape)
    minstat = np.min(np.vstack((img1.get_fdata()[inds12], img2.get_fdata()[inds12])), axis=0)

    conjunction_img = np.zeros(img1.shape)
    conjunction_img[inds12] = minstat
    os.makedirs(output_dir, exist_ok=True)
    nib.save(nib.Nifti1Image(conjunction_img, img1.affine), op.join(output_dir, 'conjunction.nii.gz'))

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

    #the FWE corrected map contains a -logp value for each cluster, so we need to manually theshold the map to account for that
    z_img_logp = nib.load(op.join(output_dir, '{prefix}_logp_level-cluster_corr-FWE_method-montecarlo.nii.gz'.format(prefix=prefix)))
    z_img = nib.load(op.join(output_dir, '{prefix}_z.nii.gz'.format(prefix=prefix)))
    z_img_thresh = thresh_img(z_img_logp, z_img, 0.05)
    nib.save(z_img_thresh, op.join(output_dir, '{prefix}_z_corr-FWE_thresh-05.nii.gz'.format(prefix=prefix)))

def run_subtraction(dset1, dset2, output_dir, prefix1, prefix2):
    #define the ALE algorithm
    alesubtraction = nimare.meta.ALESubtraction(n_iters=10000)

    #calculate ale statistics
    alesubtraction.fit(dset1, dset2)
    alesubtraction.results.maps['p'] = alesubtraction.results.maps['p_desc-group1MinusGroup2']

    #sres.save_maps(output_dir=op.join(output_dir, '{}_gt_{}'.format(prefix1, prefix2)), prefix='{}_gt_{}'.format(prefix1, prefix2))
    #run FDR multiple comparisons correction
    corr = nimare.correct.FDRCorrector(method="indep")
    cres = corr.transform(alesubtraction.results)

    #save output maps and dataset
    os.makedirs(output_dir, exist_ok=True)
    cres.save_maps(output_dir=op.join(output_dir, '{}_gt_{}'.format(prefix1, prefix2)), prefix='{}_gt_{}'.format(prefix1, prefix2))

    if 1 == 0:
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

def dataset2table(dset, output_dir, prefix):
    dset.metadata.to_csv(op.join(output_dir, '{}.csv'.format(prefix)), sep=',', index=False)

def get_peaks(img_file, output_dir):

    # get cluster + peak information from image
    stat_img = nib.load(img_file)
    atlas=['aal']
    voxel_thresh = np.min(stat_img.get_fdata()[np.nonzero(stat_img.get_fdata())])
    direction = 'pos'
    cluster_extent = 1
    prob_thresh = 5
    min_distance = 15
    out_fn = op.join(output_dir, '{0}_clusterinfo.csv'.format(op.basename(img_file).split('.')[0]))

    _, peaks_frame = atlasreader.get_statmap_info(
        stat_img, atlas=atlas, voxel_thresh=voxel_thresh,
        direction=direction, cluster_extent=cluster_extent,
        prob_thresh=prob_thresh, min_distance=min_distance)

    for i, row in peaks_frame.iterrows():
        tmplabel = row['aal']
        if i == 0:
            if tmplabel.split('_')[-1] in ['L', 'R']:
                hemis = [tmplabel.split('_')[-1]]
                labels = [' '.join(tmplabel.split('_')[:-1])]
            else:
                hemis = ['']
                labels = [' '.join(tmplabel.split('_'))]
        else:
            if tmplabel.split('_')[-1] in ['L', 'R']:
                hemis.append(tmplabel.split('_')[-1])
                labels.append(' '.join(tmplabel.split('_')[:-1]))
            else:
                hemis.append('')
                labels.append(' '.join(tmplabel.split('_')))

    peaks_frame['Hemisphere'] = hemis
    peaks_frame['Label'] = labels
    peaks_frame = peaks_frame.drop(columns=['aal'])
    peaks_frame = peaks_frame.rename(columns={'cluster_id': 'Cluster', 'peak_x': 'x', 'peak_y': 'y', 'peak_z': 'z', 'peak_value': 'Value', 'volume_mm': 'Volume (mm^3)'})

    # write output .csv files
    peaks_frame.to_csv(out_fn,
        index=False, float_format='%5g')


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
#save the dataset
grammatical_dset.save(op.join(output_dir, "grammatical.pkl"))
dataset2table(grammatical_dset, output_dir, 'grammatical')
grammatical_fn = op.join(output_dir, '{}_z_corr-FWE_thresh-05.nii.gz'.format('grammatical'))
get_peaks(grammatical_fn, output_dir)

#select the ungrammatical contrasts
ungrammatical_idx = dset_include.annotations['id'][dset_include.annotations['Meta-Analysis Group'] == 'Ungrammatical']
ungrammatical_dset = dset_include.slice(ungrammatical_idx)

#define output directory and run ale function defined above
output_dir = op.join(project_directory, 'derivatives', 'ungrammatical')
#un_ale(ungrammatical_dset, output_dir, 'ungrammatical')
#save the dataset
ungrammatical_dset.save(op.join(output_dir, "ungrammatical.pkl"))
dataset2table(ungrammatical_dset, output_dir, 'ungrammatical')
ungrammatical_fn = op.join(output_dir, '{}_z_corr-FWE_thresh-05.nii.gz'.format('ungrammatical'))
get_peaks(ungrammatical_fn, output_dir)

output_dir = op.join(project_directory, 'derivatives', 'conjunction')
conjunction_fn = op.join(output_dir, 'conjunction.nii.gz')
conjunction(grammatical_fn, ungrammatical_fn, output_dir)
get_peaks(conjunction_fn, output_dir)
exit()
#run subtraction analysis
output_dir = op.join(project_directory, 'derivatives')
run_subtraction(grammatical_dset, ungrammatical_dset, output_dir, 'grammatical', 'ungrammatical')
