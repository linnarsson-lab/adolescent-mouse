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


class ClusterL1C2(luigi.Task):
	"""
	Level 1 clustering
	"""
	tissue = luigi.Parameter()

	def requires(self) -> luigi.Task:
		return am.PrepareTissuePool(tissue=self.tissue)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L1_" + self.tissue + ".loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			with loompy.connect(self.input().fn) as ds:
				logging.info("Collecting valid cells")
				for (ix, selection, view) in ds.scan(items=np.where(ds.col_attrs["_Valid"] == 1)[0], axis=1, key="Accession"):
					loompy.create_append(out_file, view.layers, view.ra, view.ca)

			with loompy.connect(out_file) as ds:
				logging.info(f"Found {ds.shape[1]} valid cells")
				logging.info("Learning the manifold")
				cg.Cytograph2(n_factors=100, max_iter=100).fit(ds)
