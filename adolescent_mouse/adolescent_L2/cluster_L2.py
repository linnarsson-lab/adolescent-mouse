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
				ds = loompy.connect(clustered.fn)
				logging.info("Split/pool from " + clustered.fn)
				labels = ds.col_attrs["Class"]

				# Mask out cells that do not have the majority label of its cluster
				clusters = ds.col_attrs["Clusters"]

				def mode(x):
					return scipy.stats.mode(x)[0][0]

				majority_labels = npg.aggregate(clusters, labels, func=mode).astype('str')

				temp = []
				for ix in range(ds.shape[1]):
					if labels[ix] == self.major_class and labels[ix] == majority_labels[clusters[ix]]:
						temp.append(ix)
				logging.info("Keeping " + str(len(temp)) + " cells with majority labels")
				if len(temp) == 0:
					continue

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

				# Keep track of the gene order in the first file
				if accessions is None:
					accessions = ds.row_attrs["Accession"]
				
				ordering = np.where(ds.row_attrs["Accession"][None, :] == accessions[:, None])[1]
				for (ix, selection, vals) in ds.batch_scan(cells=cells, axis=1):
					ca = {}
					for key in ds.col_attrs:
						ca[key] = ds.col_attrs[key][selection]
					if dsout is None:
						dsout = loompy.create(out_file, vals[ordering, :], ds.row_attrs, ca)
					else:
						dsout.add_columns(vals[ordering, :], ca)

			#
			# logging.info("Poisson imputation")
			# pi = cg.PoissonImputation(k=self.k, N=self.N, n_genes=self.n_genes, n_components=self.n_components)
			# pi.impute_inplace(dsout)

			logging.info("Learning the manifold")
			ds = loompy.connect(out_file)
			if self.major_class == "Oligos":
				ml = cg.ManifoldLearning2(n_genes=self.n_genes, alpha=self.alpha)
			else:
				ml = cg.ManifoldLearning2(n_genes=self.n_genes, gtsne=self.gtsne, alpha=self.alpha)
			(knn, mknn, tsne) = ml.fit(ds)
			ds.set_edges("KNN", knn.row, knn.col, knn.data, axis=1)
			ds.set_edges("MKNN", mknn.row, mknn.col, mknn.data, axis=1)
			ds.set_attr("_X", tsne[:, 0], axis=1)
			ds.set_attr("_Y", tsne[:, 1], axis=1)

			logging.info("Clustering on the manifold")
			fname = "L2_" + self.major_class + "_" + self.tissue
			(eps_pct, min_pts) = params[fname]
			cls = cg.Clustering(method="mknn_louvain", eps_pct=eps_pct, min_pts=min_pts)
			labels = cls.fit_predict(ds)
			ds.set_attr("Clusters", labels, axis=1)
			logging.info(f"Found {labels.max() + 1} clusters")
			cg.Merger(min_distance=0.2).merge(ds)
			logging.info(f"Merged to {ds.col_attrs['Clusters'].max() + 1} clusters")
			ds.close()
		dsout.close()
