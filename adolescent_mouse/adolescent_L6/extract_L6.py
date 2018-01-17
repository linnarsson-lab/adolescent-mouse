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


class ExtractL6(luigi.Task):
	"""
	Extract cells by taxon for level 6
	"""
	rank = luigi.IntParameter()
	taxon = luigi.Parameter()

	def requires(self) -> luigi.Task:
		return am.PoolL5()

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L6_R{self.rank}_({self.taxon}).loom"))
		
	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			with loompy.connect(self.input().fn) as ds:
				cells = np.where(ds.ca[f"TaxonomyRank{self.rank}"] == self.taxon)[0]
				if cells.sum() == 0:
					raise ValueError(f"No cells found in taxon {self.taxon}!")
				for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
					loompy.create_append(out_file, view.layers, view.ra, view.ca)
			logging.info("Renumbering the clusters")
			with loompy.connect(out_file) as dsout:
				# Renumber the clusters
				dsout.ca.Clusters = LabelEncoder().fit_transform(dsout.ca.Clusters)

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
