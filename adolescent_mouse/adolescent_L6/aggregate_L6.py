from typing import *
import os
import csv
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import scipy.cluster.hierarchy as hierarchy
import numpy_groupies.aggregate_numpy as npg
import scipy.cluster.hierarchy as hc
import adolescent_mouse as am


class AggregateL6(luigi.Task):
	"""
	Aggregate all clusters in a new Loom file
	"""
	n_markers = luigi.IntParameter(default=10)
	n_auto_genes = luigi.IntParameter(default=6)
	rank = luigi.IntParameter()
	taxon = luigi.Parameter()

	def requires(self) -> List[luigi.Task]:
		return am.ExtractL6(rank=self.rank, taxon=self.taxon)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L6_R{self.rank}_({self.taxon}).agg.loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			logging.info("Aggregating loom file")
			ds = loompy.connect(self.input().fn)
			spec = {
				"Age": "tally",
				"Clusters": "first",
				"Class": "mode",
				"_Total": "mean",
				"Sex": "tally",
				"Tissue": "tally",
				"SampleID": "tally",
				"TissuePool": "first",
				"Outliers": "mean",
				"Bucket": "mode",
				"Region": "first",
				"OriginalClusters": "first",
				"Probable_location": "first",
				"Developmental_compartment": "first",
				"Description": "first",
				"Location_based_on": "first",
				"Neurotransmitter": "first",
				"LeafOrder": "first",
				"Comment": "first",
				"ClusterName": "first",
				"TaxonomyRank1": "first",
				"TaxonomyRank2": "first",
				"TaxonomyRank3": "first",
				"TaxonomyRank4": "first",
				"TaxonomySymbol": "first"
			}
			cg.Aggregator(f=[0.2, 0.05]).aggregate(ds, out_file, agg_spec=spec)
			dsagg = loompy.connect(out_file)

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
