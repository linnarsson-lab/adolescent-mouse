from typing import *
import os
from shutil import copyfile
import numpy as np
#import logging
import luigi
import gc
import cytograph as cg
import loompy
import logging
from scipy import sparse
from scipy.special import polygamma
from sklearn.cluster import AgglomerativeClustering, KMeans, Birch
from sklearn.decomposition import PCA, IncrementalPCA, FastICA
from sklearn.manifold import TSNE
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import BallTree, NearestNeighbors, kneighbors_graph
from sklearn.preprocessing import scale
from sklearn.svm import SVR
from scipy.stats import ks_2samp
import networkx as nx
import hdbscan
from sklearn.cluster import DBSCAN
import adolescent_mouse as am


class ClusterL2C2(luigi.Task):
	"""
	Level 2 clustering
	"""
	major_class = luigi.Parameter()
	tissue = luigi.Parameter()

	def requires(self) -> luigi.Task:
		tissues = cg.PoolSpec().tissues_for_project("Adolescent")
		if self.tissue == "All":
			return [am.ClusterL1C2(tissue=tissue) for tissue in tissues]
		else:
			return [am.ClusterL1C2(tissue=self.tissue)]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L2_" + self.major_class + "_" + self.tissue + ".loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			for clustered in self.input():
				with loompy.connect(clustered.fn, "r") as ds:
					logging.info("Split/pool from " + clustered.fn)
					cells = np.where(ds.ca.Class == self.major_class)[0]
					if self.major_class == "Oligos":
						# Special selection of cells for the oligo class, to balance between tissues
						enough_genes = ds.map((np.count_nonzero,), axis=1)[0] > 1000
						has_pdgfra = ds[ds.ra.Gene == "Pdgfra", :][0] > 0
						has_meg3 = ds[ds.ra.Gene == "Meg3", :][0] > 0
						is_doublet = np.zeros(ds.shape[1], dtype='bool')
						for g in ['Stmn2', 'Aqp4', 'Gja1', 'C1qc', 'Aif1', 'Cldn5', 'Fn1', 'Hbb-bt', 'Hbb-bh1', 'Hbb-bh2', 'Hbb-y', 'Hbb-bs', 'Hba-a1', 'Hba-a2', 'Hba-x']:
							is_doublet = np.logical_or(is_doublet, ds[ds.ra.Gene == g, :][0] > 0)
						ok_cells = enough_genes & (~is_doublet) & (has_pdgfra | ~has_meg3)
						cells = np.intersect1d(cells, np.where(ok_cells)[0])
						if cells.shape[0] > 5000:
							cells = np.random.choice(cells, 5000, False)

						for (_, _, view) in ds.scan(items=cells, axis=1, key="Accession"):
							loompy.create_append(out_file, view.layers, view.ra, view.ca)
					else:
						for (_, _, view) in ds.scan(items=cells, axis=1, key="Accession"):
							loompy.create_append(out_file, view.layers, view.ra, view.ca)

			with loompy.connect(out_file) as ds:
				logging.info(f"Found {ds.shape[1]} valid cells")
				logging.info("Learning the manifold")
				cg.Cytograph2(max_iter=100).fit(ds)
