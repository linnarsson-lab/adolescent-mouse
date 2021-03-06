from typing import *
import os
import csv
#import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import scipy.cluster.hierarchy as hierarchy
import numpy_groupies.aggregate_numpy as npg
import scipy.cluster.hierarchy as hc
import adolescent_mouse as am


class AggregateL2C2(luigi.Task):
	"""
	Aggregate all clusters in a new Loom file
	"""
	major_class = luigi.Parameter()
	tissue = luigi.Parameter()
	n_auto_genes = luigi.IntParameter(default=6)

	def requires(self) -> List[luigi.Task]:
		return am.ClusterL2C2(tissue=self.tissue, major_class=self.major_class)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L2_" + self.major_class + "_" + self.tissue + ".agg.loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			with loompy.connect(self.input().fn) as ds:
				cg.Aggregator().aggregate(ds, out_file)
				with loompy.connect(out_file) as dsagg:
					logging.info("Computing auto-annotation")
					aa = cg.AutoAnnotator(root=am.paths().autoannotation)
					aa.annotate_loom(dsagg)
					aa.save_in_loom(dsagg)

					logging.info("Computing auto-auto-annotation")
					n_clusters = dsagg.shape[1]
					(selected, selectivity, specificity, robustness) = cg.AutoAutoAnnotator(n_genes=self.n_auto_genes).fit(dsagg)
					dsagg.set_attr("MarkerGenes", np.array([" ".join(ds.ra.Gene[selected[:, ix]]) for ix in np.arange(n_clusters)]), axis=1)
					np.set_printoptions(precision=1, suppress=True)
					dsagg.set_attr("MarkerSelectivity", np.array([str(selectivity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
					dsagg.set_attr("MarkerSpecificity", np.array([str(specificity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
					dsagg.set_attr("MarkerRobustness", np.array([str(robustness[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
					dsagg.close()
