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
import matplotlib.pyplot as plt
from matplotlib import cm
import adolescent_mouse as am


class ExportL5(luigi.Task):
	"""
	Luigi Task to export summary files
	"""
	n_markers = luigi.IntParameter(default=10)

	def requires(self) -> List[luigi.Task]:
		return [am.AggregateL5(), am.PoolL5()]

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L5_All_exported"))

	def run(self) -> None:
		logging = cg.logging(self, True)
		with self.output().temporary_path() as out_dir:
			logging.info("Exporting cluster data")
			if not os.path.exists(out_dir):
				os.mkdir(out_dir)
			with loompy.connect(self.input()[0].fn) as dsagg:
				with open(os.path.join(out_dir, "L5_All_taxon_enrichment_0.2.txt"), 'w') as f:
					logging.info("Computing taxon enrichment")
					for rank in [1, 2, 3, 4]:
						taxa = list(set(dsagg.ca[f"TaxonomyRank{rank}"]))
						for taxon in taxa:
							gix = np.where(np.all(dsagg["trinaries"][:, dsagg.ca[f"TaxonomyRank{rank}"] == taxon] > 0.999, axis=1))[0]
							non_group_mean = np.mean(dsagg["trinaries"][gix, :][:, dsagg.ca[f"TaxonomyRank{rank}"] != taxon], axis=1)
							genes = dsagg.ra.Gene[gix[np.argsort(non_group_mean)]][0:20]
							f.write(str(rank) + " " + taxon + "\t" + "\t".join(genes) + "\n")
				with open(os.path.join(out_dir, "L5_All_taxon_enrichment_0.05.txt"), 'w') as f:
					logging.info("Computing taxon enrichment")
					for rank in [1, 2, 3, 4]:
						taxa = list(set(dsagg.ca[f"TaxonomyRank{rank}"]))
						for taxon in taxa:
							gix = np.where(np.all(dsagg["trinaries_0.05"][:, dsagg.ca[f"TaxonomyRank{rank}"] == taxon] > 0.999, axis=1))[0]
							non_group_mean = np.mean(dsagg["trinaries_0.05"][gix, :][:, dsagg.ca[f"TaxonomyRank{rank}"] != taxon], axis=1)
							genes = dsagg.ra.Gene[gix[np.argsort(non_group_mean)]][0:20]
							f.write(str(rank) + " " + taxon + "\t" + "\t".join(genes) + "\n")

				dsagg.export(os.path.join(out_dir, "L5_All_expression.tab"))
				dsagg.export(os.path.join(out_dir, "L5_All_enrichment.tab"), layer="enrichment")
				dsagg.export(os.path.join(out_dir, "L5_All_enrichment_q.tab"), layer="enrichment_q")
				dsagg.export(os.path.join(out_dir, "L5_All_trinaries.tab"), layer="trinaries")

			logging.info("Plotting all cells t-SNE")
			with loompy.connect(os.path.join(out_dir, self.input()[1].fn)) as ds:
				fig = plt.figure(figsize=(3, 3))
				ax = fig.add_axes([0, 0, 1, 1])
				ax.axis('off')
				colors = cg.colorize(np.arange(52))
				ix = 0
				for taxon in np.unique(ds.ca.TaxonomyRank3):
					cells = ds.ca.TaxonomyRank3 == taxon
					plt.scatter(x=ds.ca._X[cells], y=ds.ca._Y[cells], s=10, c=colors[ix, :], marker='.', label=taxon, alpha=0.3, lw=0)
					ix += 1
					lgnd = ax.legend(fontsize=10, labelspacing=0.2, loc="upper left", bbox_to_anchor=(1, 1), frameon=False)
					for handle in lgnd.legendHandles:
						handle.set_sizes([250])
						handle.set_alpha(1)
				plt.savefig(os.path.join(out_dir, "L5_All.png"), dpi=600, transparent=True, bbox_extra_artists=(lgnd,), bbox_inches='tight')
				plt.close()

