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


class ExportL2C2(luigi.Task):
	"""
	Luigi Task to export summary files
	"""
	major_class = luigi.Parameter()
	tissue = luigi.Parameter()

	def requires(self) -> List[luigi.Task]:
		return [
			am.AggregateL2C2(tissue=self.tissue, major_class=self.major_class),
			am.ClusterL2C2(tissue=self.tissue, major_class=self.major_class)
		]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L2_" + self.major_class + "_" + self.tissue + "_exported"))

	def run(self) -> None:
		logging = cg.logging(self, True)
		with self.output().temporary_path() as out_dir:
			logging.info("Exporting cluster data")
			if not os.path.exists(out_dir):
				os.mkdir(out_dir)

			with loompy.connect(self.input()[0].fn) as dsagg:
				logging.info("Exporting tab files")
				dsagg.export(os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_expression.tab"))
				dsagg.export(os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_enrichment.tab"), layer="enrichment")
				dsagg.export(os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_enrichment_q.tab"), layer="enrichment_q")
				dsagg.export(os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_trinaries.tab"), layer="trinaries")

				logging.info("Plotting manifold graph with auto-annotation")
				tags = list(dsagg.col_attrs["AutoAnnotation"])
				with loompy.connect(self.input()[1].fn) as ds:
					cg.plot_graph(ds, os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_manifold.aa.png"), tags)

					logging.info("Plotting manifold graph with auto-auto-annotation")
					tags = list(dsagg.col_attrs["MarkerGenes"][np.argsort(dsagg.col_attrs["Clusters"])])
					cg.plot_graph(ds, os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_manifold.aaa.png"), tags)

					logging.info("Plotting manifold graph with classes")
					cg.plot_classes(ds, os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_manifold.classes.png"))

					logging.info("Plotting marker heatmap")
					cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=10, out_file=os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_heatmap.pdf"))

					logging.info("Plotting latent factors")
					cg.plot_factors(ds, base_name=os.path.join(out_dir, "L2_" + self.major_class + "_" + self.tissue + "_factors"))
