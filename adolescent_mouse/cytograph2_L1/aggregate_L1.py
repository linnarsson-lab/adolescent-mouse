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


class AggregateL1C2(luigi.Task):
	"""
	Aggregate all clusters in a new Loom file
	"""
	tissue = luigi.Parameter()
	a = luigi.FloatParameter(default=1)
	b = luigi.FloatParameter(default=5)
	c = luigi.FloatParameter(default=1)
	d = luigi.FloatParameter(default=5)
	n_factors = luigi.IntParameter(default=100)
	k_smoothing = luigi.IntParameter(default=10)
	k = luigi.IntParameter(default=25)
	log = luigi.BoolParameter(default=False)
	normalize = luigi.BoolParameter(default=False)
	accel = luigi.BoolParameter(default=False)

	def requires(self) -> List[luigi.Task]:
		return am.ClusterL1C2(tissue=self.tissue, a=self.a, b=self.b, c=self.c, d=self.d, 
							n_factors=self.n_factors, k_smoothing=self.k_smoothing, k=self.k, 
							log=self.log, normalize=self.normalize, accel=self.accel)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L1_{self.tissue}_nfactors={self.n_factors}_k={self.k}_ksmoothing={self.k_smoothing}_a={self.a}_b={self.b}_c={self.c}_d={self.d}_log={self.log}_normalize={self.normalize}_accel={self.accel}.agg.loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			with loompy.connect(self.input().fn, mode="r+") as ds:
				cg.Aggregator().aggregate(ds, out_file)
				with loompy.connect(out_file) as dsagg:
					logging.info("Computing auto-annotation")
					aa = cg.AutoAnnotator(root=am.paths().autoannotation)
					aa.annotate_loom(dsagg)
					aa.save_in_loom(dsagg)

					logging.info("Computing auto-auto-annotation")
					n_clusters = dsagg.shape[1]
					(selected, selectivity, specificity, robustness) = cg.AutoAutoAnnotator(n_genes=6).fit(dsagg)
					dsagg.set_attr("MarkerGenes", np.array([" ".join(ds.ra.Gene[selected[:, ix]]) for ix in np.arange(n_clusters)]), axis=1)
					np.set_printoptions(precision=1, suppress=True)
					dsagg.set_attr("MarkerSelectivity", np.array([str(selectivity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
					dsagg.set_attr("MarkerSpecificity", np.array([str(specificity[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
					dsagg.set_attr("MarkerRobustness", np.array([str(robustness[:, ix]) for ix in np.arange(n_clusters)]), axis=1)
					dsagg.close()
