from typing import *
import os
import csv
#import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import numpy_groupies.aggregate_numpy as npg
import scipy.stats
from sklearn.cluster import DBSCAN
from sklearn.neighbors import BallTree, NearestNeighbors, kneighbors_graph
import tempfile
import adolescent_mouse as am


params = {  # eps_pct and min_pts
    "L2_Neurons_Amygdala": [80, 10],
    "L2_Neurons_Cerebellum": [90, 10],
    "L2_Neurons_Cortex1": [75, 10],
    "L2_Neurons_Cortex2": [75, 10],
    "L2_Neurons_Cortex3": [80, 10],
    "L2_Neurons_DRG": [70, 10],
    "L2_Neurons_Enteric": [60, 10],
    "L2_Neurons_Hippocampus": [90, 10],
    "L2_Neurons_Hypothalamus": [75, 10],
    "L2_Neurons_Medulla": [70, 10],
    "L2_Neurons_MidbrainDorsal": [80, 10],
    "L2_Neurons_MidbrainVentral": [75, 10],
    "L2_Neurons_Olfactory": [75, 10],
    "L2_Neurons_Pons": [75, 10],
    "L2_Neurons_SpinalCord": [90, 10],
    "L2_Neurons_StriatumDorsal": [80, 10],
    "L2_Neurons_StriatumVentral": [75, 10],
    "L2_Neurons_Sympathetic": [60, 10],
    "L2_Neurons_Thalamus": [75, 10],
    "L2_Oligos_All": [95, 500],
    "L2_Astrocytes_All": [70, 40],
    "L2_Ependymal_All": [70, 40],
    "L2_Blood_All": [70, 20],
    "L2_Immune_All": [80, 40],
    "L2_PeripheralGlia_All": [80, 20],
    "L2_Vascular_All": [90, 10]
}


class ClusterL2(luigi.Task):
	"""
	Level 2 clustering of the adolescent dataset
	"""
	major_class = luigi.Parameter()
	tissue = luigi.Parameter(default="All")
	n_genes = luigi.IntParameter(default=1000)
	n_components = luigi.IntParameter(default=30)
	k = luigi.IntParameter(default=5)
	N = luigi.IntParameter(default=5000)
	gtsne = luigi.BoolParameter(default=True)
	alpha = luigi.FloatParameter(default=1)

	def requires(self) -> luigi.Task:
		tissues = cg.PoolSpec().tissues_for_project("Adolescent")
		if self.tissue == "All":
			return [am.ClusterL1(tissue=tissue) for tissue in tissues]
		else:
			return [am.ClusterL1(tissue=self.tissue)]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L2_" + self.major_class + "_" + self.tissue + ".loom"))
		
	def run(self) -> None:
		logging = cg.logging(self)
		dsout = None  # type: loompy.LoomConnection
		accessions = None  # type: np.ndarray
		with self.output().temporary_path() as out_file:
			for clustered in self.input():
				with loompy.connect(clustered.fn, "r") as ds:
					logging.info("Split/pool from " + clustered.fn)

					logging.info("Masking outliers")
					min_pts = 10
					eps_pct = 80
					tsne_pos = np.vstack((ds.col_attrs["_X"], ds.col_attrs["_Y"])).transpose()

					# DBSCAN to find outliers
					nn = NearestNeighbors(n_neighbors=min_pts, algorithm="ball_tree", n_jobs=4)
					nn.fit(tsne_pos)
					knn = nn.kneighbors_graph(mode='distance')
					k_radius = knn.max(axis=1).toarray()
					epsilon = np.percentile(k_radius, eps_pct)

					clusterer = DBSCAN(eps=epsilon, min_samples=min_pts)
					labels = clusterer.fit_predict(tsne_pos)

					# Mask out cells that don't match the class of their local neighbors
					logging.info("Masking cells in bad neighborhoods")
					temp = []
					for ix in range(ds.shape[1]):
						if labels[ix] == -1:
							continue
						if ds.ca.Class[ix] == self.major_class:
							neighbors = ds.col_graphs.KNN.col[np.where(ds.col_graphs.KNN.row == ix)[0]]
							neighborhood = ds.ca.Class[neighbors] == self.major_class
							if neighborhood.sum() / neighborhood.shape[0] > 0.2:
								temp.append(ix)

					cells = np.array(temp)
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

			with loompy.connect(out_file) as dsout:
				logging.info("Learning the manifold")
				if self.major_class == "Oligos":
					ml = cg.ManifoldLearning2(n_genes=self.n_genes, alpha=self.alpha)
				else:
					ml = cg.ManifoldLearning2(n_genes=self.n_genes, gtsne=self.gtsne, alpha=self.alpha)
				(knn, mknn, tsne) = ml.fit(dsout)
				dsout.col_graphs.KNN = knn
				dsout.col_graphs.MKNN = mknn
				dsout.ca._X = tsne[:, 0]
				dsout.ca._Y = tsne[:, 1]

				logging.info("Clustering on the manifold")
				pl = cg.PolishedLouvain()
				labels = pl.fit_predict(dsout)
				dsout.ca.Clusters = labels + 1
				dsout.ca.Outliers = (labels == -1).astype('int')
				logging.info(f"Found {labels.max() + 1} clusters")
