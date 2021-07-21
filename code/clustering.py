
import os
import os.path as op
from glob import glob
import nibabel as nib
import numpy as np
from nilearn import masking
from nilearn import plotting
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import KMeans
from sklearn.metrics import davies_bouldin_score # lower better
from sklearn.metrics import calinski_harabasz_score #lower better
from sklearn.metrics import silhouette_score #higher better
from nilearn.connectome import ConnectivityMeasure
import pandas as pd


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, color_threshold=linkage_matrix[7,2], **kwargs)

project_dir = '/home/data/nbc/misc-projects/meta-analyses/meta-analysis_implicit-learning'

rois = glob(op.join(project_dir, 'derivatives', 'rois-pooled', '*.nii.gz'))
labels_dict = {'-44.0_12.0_20.0': '1) L Frontal Inf Oper',
               '-34.0_22.0_-2.0': '2) L Insula',
               '48.0_26.0_20.0': '3) R Frontal Inf Tri (1)',
               '50.0_26.0_4.0': '4) R Frontal Inf Tri (2)',
               '48.0_4.0_30.0': '5) R Precentral',
               '6.0_26.0_34.0': '6) R Cingulate Mid',
               '2.0_22.0_50.0': '7) L Supp Motor Area',
               '36.0_22.0_-4.0': '8) R Insula',
               '32.0_-66.0_38.0': '9) R Occipital Mid',
               '36.0_-50.0_48.0': '10) R Parietal Inf'}

roi_names = []
for roi in rois:
    tmp_roi_name = op.basename(roi.split('.nii.gz')[0])
    roi_names.append(tmp_roi_name)

labels = [labels_dict[tmp_roi] for tmp_roi in roi_names]
macm_data_array = []
rsfc_data_array = []
macm_data_array_thresh = []
rsfc_data_array_thresh = []

rsfc_mask_list = []
for roi_name in roi_names:
    tmp_rsfc_img = op.join(project_dir, 'derivatives', 'rsfc', '{}_z.nii.gz'.format(roi_name))
    tmpimg_rsfc = nib.load(tmp_rsfc_img)
    rsfc_mask_list.append(masking.compute_background_mask(tmpimg_rsfc))

rsfc_mask = masking.intersect_masks(rsfc_mask_list, threshold=1)

for roi_name in roi_names:
    tmp_macm_img = op.join(project_dir, 'derivatives', 'macm', roi_name, '{}_stat.nii.gz'.format(roi_name))
    tmpimg_macm = nib.load(tmp_macm_img)
    tmp_mask_macm = masking.compute_background_mask(tmpimg_macm)
    macm_data_array.append(masking.apply_mask(tmpimg_macm, tmp_mask_macm))

    tmp_rsfc_img = op.join(project_dir, 'derivatives', 'rsfc', '{}_z.nii.gz'.format(roi_name))
    tmpimg_rsfc = nib.load(tmp_rsfc_img)
    rsfc_data_array.append(masking.apply_mask(tmpimg_rsfc, rsfc_mask))

    tmp_macm_img_thresh = op.join(project_dir, 'derivatives', 'macm', roi_name, '{}_z_corr-FWE-voxel_thresh-0.001.nii.gz'.format(roi_name))
    tmpimg_macm_thresh = nib.load(tmp_macm_img_thresh)
    macm_data_array_thresh.append(masking.apply_mask(tmpimg_macm_thresh, tmp_mask_macm))

    tmp_rsfc_img_thresh = op.join(project_dir, 'derivatives', 'rsfc', '{}_tstat1_thr001.nii.gz'.format(roi_name))
    tmpimg_rsfc_thresh = nib.load(tmp_rsfc_img_thresh)
    rsfc_data_array_thresh.append(masking.apply_mask(tmpimg_rsfc_thresh, rsfc_mask))

print(np.shape(macm_data_array))
#generate correlation matrices
macm_corrmat = ConnectivityMeasure(kind='correlation').fit_transform([np.transpose(np.asarray(macm_data_array))])[0]
rsfc_corrmat = ConnectivityMeasure(kind='correlation').fit_transform([np.transpose(np.asarray(rsfc_data_array))])[0]

average_corrmat = np.average(np.stack((macm_corrmat, rsfc_corrmat), axis=2), axis=2)

plotting.plot_matrix(macm_corrmat, figure=(8, 8), labels=labels)
plt.savefig('macm_correlation_matrix.png')
plt.close()
plotting.plot_matrix(rsfc_corrmat, figure=(8, 8), labels=labels)
plt.savefig('rsfc_correlation_matrix.png')
plt.close()
plotting.plot_matrix(average_corrmat, figure=(8, 8), labels=labels)
plt.savefig('average_correlation_matrix.png')
plt.close()

hierarchical_clustering = AgglomerativeClustering(n_clusters=None, compute_full_tree=True, linkage='ward', distance_threshold=0, affinity='euclidean')
#kmeans_clustering = KMeans(n_clusters=6, random_state=0)

#macm_hierarchical_clustering = hierarchical_clustering.fit(macm_data_array)
#fig = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
#plt.title('MACM Hierarchical Clustering Dendrogram')
#plot_dendrogram(macm_hierarchical_clustering, orientation='right', labels=labels)
#plt.savefig('./macm-hierarchical_input-maps.png')
#plt.close()

#macm_kmeans_clustering = kmeans_clustering.fit(macm_data_array)
#tmp_df = pd.DataFrame()
#tmp_df['cluster'] = macm_kmeans_clustering.labels_
#tmp_df['label'] = labels
#tmp_df.to_csv('macm-kmeans_input-maps.csv')

#rsfc_hierarchical_clustering = hierarchical_clustering.fit(rsfc_data_array)
#fig = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
#plt.title('rsFC Hierarchical Clustering Dendrogram')
#plot_dendrogram(rsfc_hierarchical_clustering, orientation='right', labels=labels)
#plt.savefig('./rsfc-hierarchical_input-maps.png')
#plt.close()

#rsfc_kmeans_clustering = kmeans_clustering.fit(rsfc_data_array)
#tmp_df = pd.DataFrame()
#tmp_df['cluster'] = rsfc_kmeans_clustering.labels_
#tmp_df['label'] = labels
#tmp_df.to_csv('rsfc-kmeans_input-maps.csv')

macm_corrmat_hierarchical_clustering = hierarchical_clustering.fit(macm_corrmat)
fig = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
plt.title('MACM Hierarchical Clustering Dendrogram Correlation Matrix')
plot_dendrogram(macm_corrmat_hierarchical_clustering, orientation='right', labels=labels)
plt.savefig('./macm-hierarchical_input-corrmat.png')
plt.close()

#macm_corrmat_kmeans_clustering = kmeans_clustering.fit(macm_corrmat)
#tmp_df = pd.DataFrame()
#tmp_df['cluster'] = macm_corrmat_kmeans_clustering.labels_
#tmp_df['label'] = labels
#tmp_df.to_csv('macm-kmeans_input-corrmat.csv')

rsfc_corrmat_hierarchical_clustering = hierarchical_clustering.fit(rsfc_corrmat)
fig = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
plt.title('rsFC Hierarchical Clustering Dendrogram Correlation Matrix')
plot_dendrogram(rsfc_corrmat_hierarchical_clustering, orientation='right', labels=labels)
plt.savefig('./rsfc-hierarchical_input-corrmat.png')
plt.close()

#rsfc_corrmat_kmeans_clustering = kmeans_clustering.fit(rsfc_corrmat)
#tmp_df = pd.DataFrame()
#tmp_df['cluster'] = rsfc_corrmat_kmeans_clustering.labels_
#tmp_df['label'] = labels
#tmp_df.to_csv('rsfc-kmeans_input-corrmat.csv')

average_corrmat_hierarchical_clustering = hierarchical_clustering.fit(average_corrmat)
fig = plt.subplots(nrows=1, ncols=1, figsize=(15, 10))
plt.title('Multimodal Hierarchical Clustering Dendrogram Average Correlation Matrix')
plot_dendrogram(average_corrmat_hierarchical_clustering, orientation='right', labels=labels)
plt.savefig('./average-hierarchical_input-corrmat.png')
plt.close()

#average_corrmat_kmeans_clustering = kmeans_clustering.fit(average_corrmat)
#fused_matrix_kmeans_clustering = kmeans_clustering.fit(fused_matrix)
#tmp_df = pd.DataFrame()
#tmp_df['cluster'] = average_corrmat_kmeans_clustering.labels_
#tmp_df['cluster'] = fused_matrix_kmeans_clustering.labels_
#tmp_df['label'] = labels
#tmp_df.to_csv('average-kmeans_input-corrmat.csv')
#tmp_df.to_csv('average-kmeans_input-fused_matrix.csv')

hierarchical_clustering = AgglomerativeClustering(n_clusters=3, compute_full_tree=True, linkage='ward', distance_threshold=None, affinity='euclidean')
rsfc_corrmat_hierarchical_clustering_labels = hierarchical_clustering.fit(rsfc_corrmat).labels_
macm_corrmat_hierarchical_clustering_labels = hierarchical_clustering.fit(macm_corrmat).labels_
average_corrmat_hierarchical_clustering_labels = hierarchical_clustering.fit(average_corrmat).labels_

#make sum images
for a in range(np.max(rsfc_corrmat_hierarchical_clustering_labels)+1):
    cluster_label = '_'.join(np.asarray(labels)[np.where(rsfc_corrmat_hierarchical_clustering_labels == a)[0]])
    cluster_label = cluster_label.replace(' ', '')
    #need to use thresholded images
    tmp_mat = np.asarray(rsfc_data_array_thresh)[np.where(rsfc_corrmat_hierarchical_clustering_labels == a)[0], :]
    print(np.where(rsfc_corrmat_hierarchical_clustering_labels == a)[0])
    print(cluster_label)
    print(np.shape(tmp_mat))
    tmp_mat = np.where(tmp_mat > 0, 1, 0)
    tmp_mat_sum = np.sum(tmp_mat, axis=0)
    print(np.max(tmp_mat_sum))
    tmp_img = np.zeros(tmpimg_rsfc.shape)
    tmp_img[np.where(rsfc_mask.get_fdata())] = tmp_mat_sum
    nib.save(nib.Nifti1Image(tmp_img, rsfc_mask.affine), './rsfc_sum_cluster-{}.nii.gz'.format(cluster_label))
    tmp_mat_sum_thresh = tmp_mat_sum
    tmp_mat_sum_thresh[tmp_mat_sum_thresh < np.max(tmp_mat_sum_thresh) - 2] = 0
    tmp_img = np.zeros(tmpimg_rsfc.shape)
    tmp_img[np.where(rsfc_mask.get_fdata())] = tmp_mat_sum_thresh
    nib.save(nib.Nifti1Image(tmp_img, rsfc_mask.affine), './rsfc_sum_cluster-{}_thresh.nii.gz'.format(cluster_label))

for a in range(np.max(macm_corrmat_hierarchical_clustering_labels)+1):
    cluster_label = '_'.join(np.asarray(labels)[np.where(macm_corrmat_hierarchical_clustering_labels == a)[0]])
    cluster_label = cluster_label.replace(' ', '')
    #need to use thresholded images
    tmp_mat = np.asarray(macm_data_array_thresh)[np.where(macm_corrmat_hierarchical_clustering_labels == a)[0], :]
    tmp_mat_unthresh = np.asarray(macm_data_array)[np.where(macm_corrmat_hierarchical_clustering_labels == a)[0], :]
    print(np.where(macm_corrmat_hierarchical_clustering_labels == a)[0])
    print(cluster_label)
    print(np.shape(tmp_mat))
    tmp_mat = np.where(tmp_mat > 0, 1, 0)
    tmp_mat_sum = np.sum(tmp_mat, axis=0)
    print(np.max(tmp_mat_sum))
    tmp_img = np.zeros(tmpimg_macm.shape)
    tmp_img[np.where(tmp_mask_macm.get_fdata())] = tmp_mat_sum
    nib.save(nib.Nifti1Image(tmp_img, tmp_mask_macm.affine), './macm_sum_cluster-{}.nii.gz'.format(cluster_label))
    tmp_mat_sum_thresh = tmp_mat_sum
    tmp_mat_sum_thresh[tmp_mat_sum_thresh < np.max(tmp_mat_sum_thresh) - 2] = 0
    tmp_img = np.zeros(tmpimg_macm.shape)
    tmp_img[np.where(tmp_mask_macm.get_fdata())] = tmp_mat_sum_thresh
    nib.save(nib.Nifti1Image(tmp_img, tmp_mask_macm.affine), './macm_sum_cluster-{}_thresh.nii.gz'.format(cluster_label))
    tmp_mat_unthresh_mean = np.average(tmp_mat_unthresh, axis=0)
    tmp_img = np.zeros(tmpimg_macm.shape)
    tmp_img[np.where(tmp_mask_macm.get_fdata())] = tmp_mat_unthresh_mean
    nib.save(nib.Nifti1Image(tmp_img, tmp_mask_macm.affine), './macm_average_cluster-{}.nii.gz'.format(cluster_label))

for a in range(np.max(average_corrmat_hierarchical_clustering_labels)+1):
    cluster_label = '_'.join(np.asarray(labels)[np.where(average_corrmat_hierarchical_clustering_labels == a)[0]])
    cluster_label = cluster_label.replace(' ', '')
    #need to use thresholded images
    #tmp_mat_unthresh = np.asarray(macm_data_array)[np.where(average_corrmat_hierarchical_clustering_labels == a)[0], :]
    #tmp_mat_unthresh_mean = np.average(tmp_mat_unthresh, axis=0)
    #tmp_mat_thresh = np.asarray(macm_data_array_thresh)[np.where(average_corrmat_hierarchical_clustering_labels == a)[0], :]
    #tmp_mat_thresh_mean = np.average(tmp_mat_thresh, axis=0)
    tmp_mat_thresh = np.asarray(rsfc_data_array_thresh)[np.where(average_corrmat_hierarchical_clustering_labels == a)[0], :]
    tmp_mat_thresh_mean = np.average(tmp_mat_thresh, axis=0)
    tmp_img = np.zeros(rsfc_mask.shape)
    #tmp_img[np.where(tmp_mask_macm.get_fdata())] = tmp_mat_unthresh_mean
    tmp_img[np.where(rsfc_mask.get_fdata())] = tmp_mat_thresh_mean
    nib.save(nib.Nifti1Image(tmp_img, rsfc_mask.affine), './average_cluster-{}.nii.gz'.format(cluster_label))

sil_macm_corrmat = np.zeros(len(range(2,9,1)))
db_macm_corrmat = np.zeros(len(range(2,9,1)))
ch_macm_corrmat = np.zeros(len(range(2,9,1)))
sil_rsfc_corrmat = np.zeros(len(range(2,9,1)))
db_rsfc_corrmat = np.zeros(len(range(2,9,1)))
ch_rsfc_corrmat = np.zeros(len(range(2,9,1)))
sil_avg_corrmat = np.zeros(len(range(2,9,1)))
db_avg_corrmat = np.zeros(len(range(2,9,1)))
ch_avg_corrmat = np.zeros(len(range(2,9,1)))

for tmp_clust_count, tmp_clust in enumerate(range(2,9,1)):
    tmp_hc = AgglomerativeClustering(n_clusters=tmp_clust, compute_full_tree=True, linkage='ward', distance_threshold=None, affinity='euclidean')
    tmp_macm_corrmat_hc = tmp_hc.fit(macm_corrmat)
    sil_macm_corrmat[tmp_clust_count] = silhouette_score(macm_corrmat, tmp_macm_corrmat_hc.labels_, metric='euclidean')
    db_macm_corrmat[tmp_clust_count] = davies_bouldin_score(macm_corrmat, tmp_macm_corrmat_hc.labels_)
    ch_macm_corrmat[tmp_clust_count] = calinski_harabasz_score(macm_corrmat, tmp_macm_corrmat_hc.labels_)

    tmp_rsfc_corrmat_hc = tmp_hc.fit(rsfc_corrmat)
    sil_rsfc_corrmat[tmp_clust_count] = silhouette_score(rsfc_corrmat, tmp_rsfc_corrmat_hc.labels_, metric='euclidean')
    db_rsfc_corrmat[tmp_clust_count] = davies_bouldin_score(rsfc_corrmat, tmp_rsfc_corrmat_hc.labels_)
    ch_rsfc_corrmat[tmp_clust_count] = calinski_harabasz_score(rsfc_corrmat, tmp_rsfc_corrmat_hc.labels_)

    tmp_avg_corrmat_hc = tmp_hc.fit(average_corrmat)
    sil_avg_corrmat[tmp_clust_count] = silhouette_score(average_corrmat, tmp_avg_corrmat_hc.labels_, metric='euclidean')
    db_avg_corrmat[tmp_clust_count] = davies_bouldin_score(average_corrmat, tmp_avg_corrmat_hc.labels_)
    ch_avg_corrmat[tmp_clust_count] = calinski_harabasz_score(average_corrmat, tmp_avg_corrmat_hc.labels_)

fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))
axs[0,0].plot(range(2,9,1), sil_macm_corrmat)
axs[0,0].set_title('MACM')
axs[0,0].set(ylabel='Silhouette Coefficient')
axs[0,1].plot(range(2,9,1), sil_rsfc_corrmat)
axs[0,1].set_title('rsFC')
axs[0,2].plot(range(2,9,1), sil_avg_corrmat)
axs[0,2].set_title('Average')

axs[1,0].plot(range(2,9,1), db_macm_corrmat)
axs[1,0].set(ylabel='Davies-Bouldin Score')
axs[1,1].plot(range(2,9,1), db_rsfc_corrmat)
axs[1,2].plot(range(2,9,1), db_avg_corrmat)

axs[2,0].plot(range(2,9,1), ch_macm_corrmat)
axs[2,0].set(ylabel='Calinksi-Haragasz Score')
axs[2,1].plot(range(2,9,1), ch_rsfc_corrmat)
axs[2,2].plot(range(2,9,1), ch_avg_corrmat)

plt.savefig('./clustering-metrics.png')
plt.close()
