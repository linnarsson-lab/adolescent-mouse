from typing import *
import os
import logging
import loompy
from scipy import sparse
import numpy as np
import networkx as nx
import cytograph as cg
import luigi
from scipy.spatial.distance import squareform, pdist
import adolescent_mouse as am


class ExportByTaxonL6(luigi.Task):
	"""
	Luigi Task to export summary files by taxon
	"""
	n_markers = luigi.IntParameter(default=10)
	rank = luigi.IntParameter()

	def requires(self) -> List[luigi.Task]:
		return [am.AggregateByTaxonL6(rank=self.rank), am.ExtractByTaxonL6(rank=self.rank)]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L6_R{self.rank}_exported"))

	def run(self) -> None:
		logging = cg.logging(self, True)
		with self.output().temporary_path() as out_dir:
			logging.info("Exporting cluster data")
			if not os.path.exists(out_dir):
				os.mkdir(out_dir)
			with loompy.connect(self.input()[0].fn) as dsagg:
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_expression.tab"))
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_enrichment.tab"), layer="enrichment")
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_enrichment_q.tab"), layer="enrichment_q")
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_trinaries.tab"), layer="trinaries")

				logging.info("Plotting manifold graph with auto-annotation")
				with loompy.connect(self.input()[1].fn) as ds:
					cg.plot_graph(ds, os.path.join(out_dir, f"L6_R{self.rank}_manifold.aa.png"), list(dsagg.ca.AutoAnnotation))

					logging.info("Plotting manifold graph with auto-auto-annotation")
					cg.plot_graph(ds, os.path.join(out_dir, f"L6_R{self.rank}_manifold.aaa.png"), list(dsagg.ca.MarkerGenes))

					logging.info("Plotting manifold graph with taxon names")
					cg.plot_graph(ds, os.path.join(out_dir, f"L6_R{self.rank}_manifold.names.png"), list(dsagg.ca[f"TaxonomyRank{self.rank}"]))

					logging.info("Plotting marker heatmap")
					cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, f"L6_R{self.rank}_heatmap.pdf"))
