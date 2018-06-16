from typing import *
import os
#import logging
import loompy
from scipy import sparse
import numpy as np
import networkx as nx
import cytograph as cg
import luigi
import adolescent_mouse as am


class ExportL1C2(luigi.Task):
	"""
	Luigi Task to export summary files
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
		return [
			am.AggregateL1C2(tissue=self.tissue, a=self.a, b=self.b, c=self.c, d=self.d, 
							n_factors=self.n_factors, k_smoothing=self.k_smoothing, k=self.k, 
							log=self.log, normalize=self.normalize, accel=self.accel),
			am.ClusterL1C2(tissue=self.tissue, a=self.a, b=self.b, c=self.c, d=self.d, 
							n_factors=self.n_factors, k_smoothing=self.k_smoothing, k=self.k, 
							log=self.log, normalize=self.normalize, accel=self.accel)
		]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L1_{self.tissue}_nf={self.n_factors}_a={self.a}_b={self.b}_c={self.c}_d={self.d}_log={self.log}_normalize={self.normalize}_accel={self.accel}"))

	def run(self) -> None:
		logging = cg.logging(self, True)
		logging.info("Exporting cluster data")
		with self.output().temporary_path() as out_dir:
			if not os.path.exists(out_dir):
				os.mkdir(out_dir)
			with loompy.connect(self.input()[0].fn, mode="r+") as dsagg:
				logging.info("Exporting tab files")
				dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_expression.tab"))
				dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_enrichment.tab"), layer="enrichment")
				dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_enrichment_q.tab"), layer="enrichment_q")
				dsagg.export(os.path.join(out_dir, "L1_" + self.tissue + "_trinaries.tab"), layer="trinaries")

				ds = loompy.connect(self.input()[1].fn)

				logging.info("Plotting MKNN graph")
				cg.plot_knn(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.mknn.png"))

				# logging.info("Plotting Louvain resolution")
				# cg.plot_louvain(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.louvain.png"))

				try:
					logging.info("Plotting manifold graph with classes")
					cg.plot_classes(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.classes.png"))
				except Exception:
					pass

				logging.info("Plotting manifold graph with auto-annotation")
				tags = list(dsagg.col_attrs["AutoAnnotation"])
				cg.plot_graph(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.aa.png"), tags)

				logging.info("Plotting manifold graph with auto-auto-annotation")
				tags = list(dsagg.col_attrs["MarkerGenes"])
				cg.plot_graph(ds, os.path.join(out_dir, "L1_" + self.tissue + "_manifold.aaa.png"), tags)

				logging.info("Plotting marker heatmap")
				cg.plot_markerheatmap(ds, dsagg, out_file=os.path.join(out_dir, "L1_" + self.tissue + "_heatmap.pdf"))

				logging.info("Plotting latent factors")
				cg.plot_factors(ds, base_name=os.path.join(out_dir, "L1_" + self.tissue + "_factors"))
