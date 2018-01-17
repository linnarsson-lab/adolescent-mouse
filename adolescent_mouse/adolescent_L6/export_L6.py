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
import matplotlib.pyplot as plt


class ExportL6(luigi.Task):
	"""
	Luigi Task to export summary files
	"""
	n_markers = luigi.IntParameter(default=10)
	rank = luigi.IntParameter()
	taxon = luigi.Parameter()

	def requires(self) -> List[luigi.Task]:
		return [am.AggregateL6(rank=self.rank, taxon=self.taxon), am.ExtractL6(rank=self.rank, taxon=self.taxon)]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L6_R{self.rank}_({self.taxon})_exported"))

	def run(self) -> None:
		logging = cg.logging(self, True)
		with self.output().temporary_path() as out_dir:
			logging.info("Exporting cluster data")
			if not os.path.exists(out_dir):
				os.mkdir(out_dir)
			with loompy.connect(self.input()[0].fn) as dsagg:
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_expression.tab"))
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_enrichment.tab"), layer="enrichment")
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_enrichment_q.tab"), layer="enrichment_q")
				dsagg.export(os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_trinaries.tab"), layer="trinaries")

				logging.info("Plotting manifold graph with auto-annotation")
				with loompy.connect(self.input()[1].fn) as ds:
					cg.plot_graph(ds, os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_manifold.aa.png"), list(dsagg.ca.AutoAnnotation))

					logging.info("Plotting manifold graph with auto-auto-annotation")
					cg.plot_graph(ds, os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_manifold.aaa.png"), list(dsagg.ca.MarkerGenes))

					logging.info("Plotting manifold graph with cluster names")
					cg.plot_graph(ds, os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_manifold.names.png"), list(dsagg.ca.ClusterName))

					logging.info("Plotting marker heatmap")
					cg.plot_markerheatmap(ds, dsagg, n_markers_per_cluster=self.n_markers, out_file=os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_heatmap.pdf"))

					size = 200000 / ds.shape[1]
					fig = plt.figure(figsize=(3,3))
					ax = fig.add_axes([0, 0, 1, 1])
					ax.axis('off')
					ix = 0
					if self.rank == 3:
						colors = cg.colorize(np.unique(ds.ca.ClusterName))
						for cluster in np.unique(ds.ca.ClusterName):
							cells = ds.ca.ClusterName == cluster
							plt.scatter(x=ds.ca._X[cells], y=ds.ca._Y[cells], s=size, c=colors[ix, :], marker='.', label=cluster, alpha=0.5, lw=0)
							ix += 1
					else:
						colors = cg.colorize(np.unique(ds.ca.TaxonomyRank4))
						for taxon4 in np.unique(ds.ca.TaxonomyRank4):
							cells = ds.ca.TaxonomyRank4 == taxon4
							plt.scatter(x=ds.ca._X[cells], y=ds.ca._Y[cells], s=size, c=colors[ix, :], marker='.', label=taxon4, alpha=0.5, lw=0)
							ix += 1
					lgnd = ax.legend(fontsize=10, labelspacing=0.2, loc="upper left", bbox_to_anchor=(1, 1),frameon=False)
					for handle in lgnd.legendHandles:
						handle.set_sizes([250])
						handle.set_alpha(1)
					plt.savefig(os.path.join(out_dir, f"L6_R{self.rank}_({self.taxon})_manifold.pretty.png"), dpi=600, transparent=True, bbox_extra_artists=(lgnd,), bbox_inches='tight')
					plt.close()

