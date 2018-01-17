from typing import *
import os
import csv
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import numpy_groupies.aggregate_numpy as npg
import scipy.stats
import adolescent_mouse as am
import pandas as pd
from sklearn.preprocessing import LabelEncoder


class PoolL5(luigi.Task):
	"""
	Level 4 pooling of the adolescent dataset
	"""

	def requires(self) -> luigi.Task:
		targets = [
			'Sensory_Neurons',
			'Sympathetic_Neurons',
			'Enteric_Neurons',
			'Mesencephalon_Excitatory',
			'Diencephalon_Excitatory',
			'Hindbrain_Inhibitory',
			'SpinalCord_Inhibitory',
			'Brain_Granule',
			'Brain_CholinergicMonoaminergic',
			'Mesencephalon_Inhibitory',
			'Diencephalon_Inhibitory',
			'Striatum_MSN',
			'Hypothalamus_Peptidergic',
			'Telencephalon_Excitatory',
			'Brain_Neuroblasts',
			'Hindbrain_Excitatory',
			'SpinalCord_Excitatory',
			'Telencephalon_Inhibitory',
			'Olfactory_Inhibitory',
			'Oligos',
			'Astrocytes',
			'Vascular',
			'Immune',
			'PeripheralGlia',
			'Ependymal'
		]
		for target in targets:
			yield am.CurateL4(target=target)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L5_All.loom"))
		
	def run(self) -> None:
		logging = cg.logging(self)
		samples = [x.fn for x in self.input()]
		max_cluster_id = 0
		cluster_ids: List[int] = []
		original_ids: List[int] = []
		samples_per_cell: List[str] = []

		celltypes_summary_file = os.path.join(am.paths().build, "curated_L4", "celltypes_summary_leaforder16-Dec-2017.xlsx")
		celltypes_summary = pd.read_excel(celltypes_summary_file)
		celltypes_dict = {celltypes_summary.columns.values[i]: celltypes_summary.values[:, i] for i in range(celltypes_summary.shape[1])}

		with self.output().temporary_path() as out_file:
			accessions = None  # type: np.ndarray
			for sample in samples:
				with loompy.connect(sample) as ds:
					logging.info(f"Adding {ds.shape[1]} cells from {sample}")
					target = os.path.basename(sample)[3:-5]
					not_excluded = celltypes_dict["OriginalCluster"][celltypes_dict["Bucket"] == target]
					cells = np.where(np.isin(ds.ca.Clusters, not_excluded))[0]
					for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
						cluster_ids += list(view.ca.Clusters + max_cluster_id)
						original_ids += list(view.ca.Clusters)
						samples_per_cell += [sample] * selection.shape[0]
						loompy.create_append(out_file, view.layers, view.ra, view.ca, fill_values="auto")
					max_cluster_id = max(cluster_ids) + 1
					logging.info(f"Found {max_cluster_id} clusters total")
			with loompy.connect(out_file) as ds:
				ds.ca.Clusters = np.array(cluster_ids)
				ds.ca.OriginalClusters = np.array(original_ids)
				ds.ca.Bucket = np.array(samples_per_cell)

				leaf_order = np.zeros(ds.shape[1], dtype='int') - 1
				le = LabelEncoder()
				le.fit(celltypes_dict["ClusterName"])
				new_clusters = np.zeros(ds.shape[1], dtype='int') - 1
				d = {}
				for attr in ["LeafOrder", "Probable_location", "Developmental_compartment", "Region", "Description", "Location_based_on", "Neurotransmitter", "ClusterName", "Taxonomy_group", "Comment", "ClusterName"]:
					d[attr] = np.array([""] * ds.shape[1], dtype=object)

				for ix in range(len(celltypes_dict["Bucket"])):
					bucket = celltypes_dict["Bucket"][ix]
					bucket_name = f"/Users/sten/build_20171205/L4_{bucket}.loom"
					original_cluster = celltypes_dict["OriginalCluster"][ix]
					cells = np.logical_and(ds.ca.Bucket == bucket_name, ds.ca.OriginalClusters == original_cluster)
					leaf_order[cells] = celltypes_dict["LeafOrder"][ix]
					new_clusters[cells] = le.transform([celltypes_dict["ClusterName"][ix]])
					for attr in d.keys():
						d[attr][cells] = celltypes_dict[attr][ix]

				logging.info(f"Found {new_clusters.max() + 1} clusters total")
				ds.ca.Clusters = new_clusters
				ds.ca.LeafOrder = leaf_order
				for key, vals in d.items():
					ds.ca[key] = vals.astype("unicode")

				taxonomy_file = os.path.join(am.paths().build, "curated_L4", "Taxonomy.xlsx")
				taxonomy_table = pd.read_excel(taxonomy_file)
				taxonomy = {taxonomy_table.values[i, 3]: taxonomy_table.values[i, :] for i in range(taxonomy_table.shape[0])}

				tax1 = np.array([""] * ds.shape[1], dtype=object)
				tax2 = np.array([""] * ds.shape[1], dtype=object)
				tax3 = np.array([""] * ds.shape[1], dtype=object)
				tax4 = np.array([""] * ds.shape[1], dtype=object)
				taxs = np.array([""] * ds.shape[1], dtype=object)

				for i in range(ds.shape[1]):
					if ds.ca.Clusters[i] == -1:
						continue
					tax1[i] = taxonomy[d["Taxonomy_group"][i]][0]
					tax2[i] = taxonomy[d["Taxonomy_group"][i]][1]
					tax3[i] = taxonomy[d["Taxonomy_group"][i]][2]
					tax4[i] = taxonomy[d["Taxonomy_group"][i]][3]
					taxs[i] = taxonomy[d["Taxonomy_group"][i]][4]
				ds.ca.TaxonomyRank1 = tax1
				ds.ca.TaxonomyRank2 = tax2
				ds.ca.TaxonomyRank3 = tax3
				ds.ca.TaxonomyRank4 = tax4
				ds.ca.TaxonomySymbol = taxs

				logging.info("Recomputing the list of valid genes")
				nnz = ds.map([np.count_nonzero], axis=0)[0]
				valid_genes = np.logical_and(nnz > 10, nnz < ds.shape[1] * 0.6)
				ds.ra._Valid = valid_genes.astype('int')
				
				logging.info("Learning the manifold")
				ml = cg.ManifoldLearning2(gtsne=True, alpha=1, max_iter=3000)
				(knn, mknn, tsne) = ml.fit(ds)
				ds.col_graphs.KNN = knn
				ds.col_graphs.MKNN = mknn
				ds.ca._X = tsne[:, 0]
				ds.ca._Y = tsne[:, 1]
