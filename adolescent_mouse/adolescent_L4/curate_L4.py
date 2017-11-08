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
		L3_targets = [
			'Sensory_Neurons',
			'Sympathetic_Neurons',
			'Enteric_Neurons',
			'DiMesencephalon_Excitatory',
			'Hindbrain_Inhibitory',
			'SpinalCord_Inhibitory',
			'Brain_Granule',
			'Brain_CholinergicMonoaminergic',
			'DiMesencephalon_Inhibitory',
			'Striatum_MSN',
			'Hypothalamus_Peptidergic',
			'Forebrain_Excitatory',
			'Brain_Neuroblasts',
			'Hindbrain_Excitatory',
			'SpinalCord_Excitatory',
			'Forebrain_Inhibitory',
			'Olfactory_Inhibitory'
		]
		if self.target in L3_targets:
			return am.FilterL3(target=self.target)
		else:
			return am.FilterL2(major_class=self.target.split("_")[0], tissue="All")

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L4_" + self.target + ".loom"))
		
	def run(self) -> None:
		logging = cg.logging(self, True)
		dsout: loompy.LoomConnection = None
		with self.output().temporary_path() as out_file:
			logging.info("Curating clusters in " + self.target)
			ds = loompy.connect(self.input().fn)
			n_labels = len(set(ds.Clusters))
			curated = np.loadtxt(os.path.join(am.paths().build, "curated", self.target + ".curated.txt"), skiprows=1, usecols=(2,))

			# Renumber the clusters
			d = dict(zip(np.arange(n_labels), curated.astype('int')))
			new_clusters = np.array([d[x] if x in d else -1 for x in ds.Clusters])
			cells = np.where(new_clusters >= 0)[0]
			logging.info(f"Keeping {cells.shape[0]} of {ds.shape[1]} cells")
			for (ix, selection, vals) in ds.batch_scan(cells=cells, axis=1):
				ca = {k: v[selection] for k, v in ds.col_attrs.items()}
				ca["Clusters"] = new_clusters[selection]
				if dsout is None:
					dsout = loompy.create(out_file, vals, ds.row_attrs, ca)
				else:
					dsout.add_columns(vals, ca)

			try:
				# Filter the KNN and MKNN edges
				(a, b, w) = ds.get_edges("KNN", axis=1)
				mask = np.logical_and(np.in1d(a, cells), np.in1d(b, cells))
				a = a[mask]
				b = b[mask]
				w = w[mask]
				d = dict(zip(np.sort(cells), np.arange(cells.shape[0])))
				a = np.array([d[x] for x in a])
				b = np.array([d[x] for x in b])
				dsout.set_edges("KNN", a, b, w, axis=1)

				(a, b, w) = ds.get_edges("MKNN", axis=1)
				mask = np.logical_and(np.in1d(a, cells), np.in1d(b, cells))
				a = a[mask]
				b = b[mask]
				w = w[mask]
				d = dict(zip(np.sort(cells), np.arange(cells.shape[0])))
				a = np.array([d[x] for x in a])
				b = np.array([d[x] for x in b])
				dsout.set_edges("MKNN", a, b, w, axis=1)
			except KeyError:
				logging.info("No edges found in the input file")

			logging.info(f"Curating {n_labels} -> {len(set(new_clusters[cells]))} clusters")

			dsout.close()
			ds.close()
