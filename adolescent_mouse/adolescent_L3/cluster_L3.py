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


class ClusterL3(luigi.Task):
	"""
	Level 3 clustering of the adolescent dataset
	"""
	target = luigi.Parameter()  # e.g. Hindbrain_Inhibitory
	n_enriched = luigi.Parameter(default=500)  # Number of enriched genes per cluster to use for manifold learning
	
	def requires(self) -> Iterator[luigi.Task]:
		tissues: List[str] = []
		for fname in os.listdir(os.path.join(am.paths().build, "curated_L2")):
			if not fname.startswith("L2"):
				continue
			tissue = fname.split("_")[2]
			with open(os.path.join(am.paths().build, "curated_L2", fname)) as f:
				schedule = [x[:-1].split("\t") for x in f.readlines()]
				for (cluster, n_cells, auto_target, curated_target, comment) in schedule:
					if curated_target == self.target:
						if tissue not in tissues:
							if tissue == "All":
								yield [am.ClusterL2(tissue="All", major_class=self.target), am.AggregateL2(tissue="All", major_class=self.target)]
							else:
								yield [am.ClusterL2(tissue=tissue, major_class="Neurons"), am.AggregateL2(tissue=tissue, major_class="Neurons")]
							tissues.append(tissue)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L3_" + self.target + ".loom"))

	def run(self) -> None:
		logging = cg.logging(self, True)
		dsout: loompy.LoomConnection = None
		accessions: loompy.LoomConnection = None
		with self.output().temporary_path() as out_file:
			logging.info("Gathering cells for " + self.target)
			enriched_markers: List[np.ndarray] = []  # The enrichment vector for each selected cluster
			cells_found = False
			for in_file, agg_file in self.input():
				tissue = os.path.basename(in_file.fn).split("_")[2].split(".")[0]
				ds = loompy.connect(in_file.fn)
				dsagg = loompy.connect(agg_file.fn)
				enrichment = dsagg.layer["enrichment"][:, :]
				labels = ds.col_attrs["Clusters"]
				ordering: np.ndarray = None
				logging.info(tissue)

				# Figure out which cells should be collected
				cells: List[int] = []
				for fname in os.listdir(os.path.join(am.paths().build, "curated_L2")):
					if not fname.startswith("L2"):
						continue
					from_tissue = fname.split("_")[2]
					if from_tissue != tissue:
						continue
					if tissue == "All":
						major_class = fname.split("_")[1]
						if major_class != self.target:
							continue
					logging.info("Gathering cells from " + in_file.fn)
					logging.info("Gathering cells based on " + fname)
					with open(os.path.join(am.paths().build, "curated_L2", fname)) as f:
						schedule = [x[:-1].split("\t") for x in f.readlines()]
						for (cluster_str, n_cells, auto_target, curated_target, comment) in schedule:
							cluster = int(cluster_str)
							if curated_target == self.target:
								if accessions is None:
									accessions = ds.row_attrs["Accession"]
								if ordering is None:
									ordering = np.where(ds.row_attrs["Accession"][None, :] == accessions[:, None])[1]
								cells += list(np.where(labels == cluster)[0])
								enriched_markers.append(np.argsort(-enrichment[:, cluster][ordering]))

				if len(cells) > 0:
					cells = np.sort(np.array(cells))
					cells_found = True
					for (ix, selection, view) in ds.scan(items=cells, axis=1, key="Accession"):
						loompy.create(out_file, view.layers, view.ra, view.ca)

			if not cells_found:
				raise ValueError(f"No cells matched any schedule for {self.target}")

			# Figure out which enriched markers to use
			ix = 0
			temp: List[int] = []
			while len(temp) < self.n_enriched:
				for j in range(len(enriched_markers)):
					if enriched_markers[j][ix] not in temp:
						temp.append(enriched_markers[j][ix])
				ix += 1
			genes = np.sort(np.array(temp))

			logging.info("Learning the manifold")
			with loompy.connect(out_file) as dsout:
				ml = cg.ManifoldLearning2(gtsne=True, alpha=1, genes=genes)
				(knn, mknn, tsne) = ml.fit(dsout)
				dsout.col_graphs.KNN = knn
				dsout.col_graphs.MKNN = mknn
				dsout.ca._X = tsne[:, 0]
				dsout.ca._Y = tsne[:, 1]

				logging.info("Clustering on the manifold")

				special_res = {
					"Astrocytes": 0.6,
					"Sensory_Neurons": 0.35,
					"Brain_Granule": 0.6
				}
				r = 1.0
				if self.target in special_res:
					r = special_res[self.target]

				pl = cg.PolishedLouvain(resolution=r)
				labels = pl.fit_predict(dsout)
				dsout.ca.Clusters = labels + 1
				dsout.ca.Outliers = (labels == -1).astype('int')
				logging.info(f"Found {labels.max() + 1} clusters")
