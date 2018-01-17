from typing import *
import os
from shutil import copyfile
import csv
#import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import numpy_groupies.aggregate_numpy as npg
import scipy.stats
import adolescent_mouse as am


class CurateL4(luigi.Task):
	"""
	Level 4 manual curation of the adolescent dataset
	"""
	target = luigi.Parameter()  # e.g. L3_Forebrain_Excitatory
	
	def requires(self) -> Iterator[luigi.Task]:
		return am.ClusterL3(target=self.target)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L4_" + self.target + ".loom"))
		
	def run(self) -> None:
		logging = cg.logging(self, True)
		with self.output().temporary_path() as out_file:
			logging.info("Curating clusters in " + self.target)
			with loompy.connect(self.input().fn) as ds:
				n_labels = len(set(ds.ca.Clusters))
				with open(os.path.join(am.paths().build, "curated_L3", "L3_" + self.target + "_agg_cluster.curated.txt")) as f:
					lines = f.readlines()
					curated = [x[:-1].split("\t") for x in lines[1:]]
				old_cluster_ids = [int(x[1]) for x in curated]
				new_cluster_ids = [int(x[2]) for x in curated]
				# Renumber the clusters
				d = dict(zip(old_cluster_ids, new_cluster_ids))
				new_clusters = np.array([d[x] for x in ds.ca.Clusters])

				logging.info("Renumbering clusters according to manual curation")
				for items in curated:
					if len(items) >= 4 and items[3] != "":
						gene = items[3]
						new_cluster = int(items[2])
						# We have a cluster that needs to be split by some gene expression
						logging.info(f"Splitting cluster {items[1]} -> {new_cluster} by {gene}")
						if gene not in ds.ra.Gene:
							raise ValueError(f"Gene {gene} not found!")
						zeros = np.where(np.logical_and(ds[ds.ra.Gene == gene, :][0] == 0, new_clusters == new_cluster))[0]
						new_clusters[zeros] = -1
				outliers = (new_clusters == -1).astype('int')
				cells = np.where(new_clusters >= 0)[0]
				logging.info(f"Keeping {cells.shape[0]} of {ds.shape[1]} cells")
				dsout: loompy.LoomConnection = None
				for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
					loompy.create(out_file, view.layers, view.ra, view.ca)

				with loompy.connect(out_file) as dsout:
					dsout.ca.Clusters = new_clusters[cells]
					dsout.ca.Outliers = outliers[cells]
					d = dict(zip(ds.ca.Clusters[cells], dsout.ca.Clusters))
					for (old, new) in d.items():
						logging.info(f"{old} -> {new}")

					logging.info("Recomputing the list of valid genes")
					nnz = dsout.map([np.count_nonzero], axis=0)[0]
					valid_genes = np.logical_and(nnz > 10, nnz < dsout.shape[1] * 0.6)
					dsout.ra._Valid = valid_genes.astype('int')
					
					logging.info("Learning the manifold")
					ml = cg.ManifoldLearning2(gtsne=True, alpha=1)
					(knn, mknn, tsne) = ml.fit(dsout)
					dsout.col_graphs.KNN = knn
					dsout.col_graphs.MKNN = mknn
					dsout.ca._X = tsne[:, 0]
					dsout.ca._Y = tsne[:, 1]

					labels = new_clusters
					logging.info(f"Curating {n_labels} -> {labels.max() + 1} clusters")
