from typing import *
import os
import csv
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import numpy_groupies.aggregate_numpy as npg
import scipy.stats
import tempfile
import adolescent_mouse as am


class FilterL3(luigi.Task):
	"""
	Level 3 filtering of bad clusters
	"""
	target = luigi.Parameter()

	def requires(self) -> luigi.Task:
		return am.ClusterL3(target=self.target)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L3_" + self.target + ".filtered.loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		logging.info("Filtering bad clusters")
		dsout = None  # type: loompy.LoomConnection
		accessions = None  # type: np.ndarray
		with self.output().temporary_path() as out_file:
			ds = loompy.connect(self.input().fn)
			if not len(set(ds.Clusters)) == ds.Clusters.max() + 1:
				raise ValueError("There are holes in the cluster ID sequence!")
			labels = ds.Clusters
			n_labels = len(set(labels))
			remove = []

			# Remove clusters that lack enriched genes
			logging.info("Checking for clusters with no enriched genes")
			logging.info("Trinarizing")
			trinaries = cg.Trinarizer().fit(ds)
			logging.info("Computing cluster gene enrichment scores")
			(markers, enrichment, qvals) = cg.MarkerSelection(10).fit(ds)
			data = trinaries[markers, :].T
			for ix in range(n_labels):
				total_score = data[ix, ix * 10:(ix + 1) * 10].sum()
				if total_score < 2:
					remove.append(ix)
					logging.info(f"Cluster {ix} score: {total_score:.1f} < 2 (removing).")
				else:
					logging.info(f"Cluster {ix} score: {total_score:.1f}")

			# Remove clusters that express genes of the wrong major class
			logging.info("Checking for clusters with markers of wrong major class")
			nix_genes = {
				"L3_SpinalCord_Inhibitory": [],
				"L3_SpinalCord_Excitatory": [],
				"L3_Olfactory_Inhibitory": [],
				"L3_Enteric_Neurons": ["Sostdc1"],
				"L3_Sensory_Neurons": ["Mpz", "Sostdc1"],
				"L3_Sympathetic_Neurons": ["Mpz", "Sostdc1"],
				"L3_Hypothalamus_Peptidergic": [],
				"L3_Hindbrain_Inhibitory": [],
				"L3_Hindbrain_Excitatory": [],
				"L3_Brain_Neuroblasts": ["Vtn"],
				"L3_Forebrain_Inhibitory": [],
				"L3_Forebrain_Excitatory": ["Mog"],
				"L3_DiMesencephalon_Inhibitory": ["Cd9", "Aqp4"],
				"L3_DiMesencephalon_Excitatory": ["Aqp4"],
				"L3_Brain_Granule": [],
				"L3_Brain_CholinergicMonoaminergic": [],
				"L3_Striatum_MSN": []
			}
			for lbl in range(n_labels):
				# Clusters with markers of other major class
				n_cells_in_cluster = (ds.Clusters == lbl).sum()
				for gene in nix_genes["L3_" + self.target]:
					if gene not in ds.Gene:
						logging.warn("Couldn't use '" + gene + "' to nix clusters")
					gix = np.where(ds.Gene == gene)[0][0]
					if trinaries[gix, lbl] > 0.95:
						logging.info("Nixing cluster {} because {} was detected".format(lbl, gene))
						remove.append(lbl)

			retain = np.sort(np.setdiff1d(np.arange(n_labels), remove))
			temp: List[int] = []
			for i in retain:
				temp += list(np.where(ds.Clusters == i)[0])
			cells = np.sort(np.array(temp))

			# Renumber the clusters
			d = dict(zip(retain, np.arange(len(set(retain)) + 1)))
			new_clusters = np.array([d[x] if x in d else -1 for x in ds.Clusters])
			logging.info(f"Keeping {cells.shape[0]} of {ds.shape[1]} cells")
			for (ix, selection, vals) in ds.batch_scan(cells=cells, axis=1):
				ca = {k: v[selection] for k, v in ds.col_attrs.items()}
				ca["Clusters"] = new_clusters[selection]
				if dsout is None:
					dsout = loompy.create(out_file, vals, ds.row_attrs, ca)
				else:
					dsout.add_columns(vals, ca)

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

			logging.info(f"Filtering {n_labels} -> {len(set(new_clusters[cells]))} clusters")
