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
import scipy.sparse as sparse


def _gene_selection_L5(dsagg: loompy.LoomConnection) -> Tuple[np.ndarray, np.ndarray]:
	trinaries = dsagg["trinaries_0.05"][:, :]
	# Find non-neuronal genes
	neurons = (dsagg.ca.TaxonomyRank1 == "Neurons")
	astrocytes = (dsagg.ca.TaxonomyRank3 == "Astroependymal cells")
	oligos = (dsagg.ca.TaxonomyRank3 == "Oligodendrocytes")
	vascular = (dsagg.ca.TaxonomyRank3 == "Vascular cells")
	immune = (dsagg.ca.TaxonomyRank3 == "Immune cells")
	arcane = (dsagg.ca.TaxonomyRank3 == "Neural crest-like glia")
	nng = (np.mean(trinaries[:, astrocytes], axis=1) == 1) & (np.mean(trinaries[:, neurons], axis=1) >= 0.9)
	nng = nng | (np.mean(trinaries[:, oligos], axis=1) == 1) & (np.mean(trinaries[:, neurons], axis=1) >= 0.9)
	nng = nng | (np.mean(trinaries[:, vascular], axis=1) == 1) & (np.mean(trinaries[:, neurons], axis=1) >= 0.9)
	nng = nng | (np.mean(trinaries[:, immune], axis=1) == 1) & (np.mean(trinaries[:, neurons], axis=1) >= 0.9)
	nng = nng | (np.mean(trinaries[:, arcane], axis=1) == 1) & (np.mean(trinaries[:, neurons], axis=1) >= 0.9)

	# Block specific troublesome genes
	blocked = np.isin(dsagg.ra.Gene, ['Xist','Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d','Trf','Mobp','Mfge8','Mbp','Junb','Hbb-bs','H2-DMb2','Fos','Plp1','Mog','Egr1','Jun'])
	# Remove housekeeping genes
	blocked = blocked | ((trinaries >= 0.99).mean(axis=1) > 0.8)

	return (nng, blocked)


class AggregateL5(luigi.Task):
	"""
	Aggregate all clusters in a new Loom file
	"""
	n_markers = luigi.IntParameter(default=10)
	n_auto_genes = luigi.IntParameter(default=6)

	def requires(self) -> List[luigi.Task]:
		return am.PoolL5()

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L5_All.agg.loom"))

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
				"LeafOrder": "first",
				"Probable_location": "first",
				"Developmental_compartment": "first",
				"Description": "first",
				"Location_based_on": "first",
				"Neurotransmitter": "first",
				"ClusterName": "first",
				"Comment": "first",
				"ClusterName": "first",
				"TaxonomyRank1": "first",
				"TaxonomyRank2": "first",
				"TaxonomyRank3": "first",
				"TaxonomyRank4": "first",
				"TaxonomySymbol": "first"
			}
			cg.Aggregator(f=[0.2, 0.05]).aggregate(ds, out_file, agg_spec=spec)

			with loompy.connect(out_file) as dsagg:
				logging.info("Finding non-neuronal, housekeeping, and troublemaking genes")
				(nng, blocked) = _gene_selection_L5(dsagg)

				logging.info("Manifold learning on the aggregate file")
				normalizer = cg.Normalizer(False)
				normalizer.fit(dsagg)
				pca = cg.PCAProjection(np.arange(dsagg.shape[1] * 10), max_n_components=50)
				pca.fit(dsagg, normalizer)
				transformed = pca.transform(dsagg, normalizer)
				k = 40
				bnn = cg.BalancedKNN(k=k, maxl=2*k)
				bnn.fit(transformed)
				knn = bnn.kneighbors(mode='connectivity')[1][:, 1:]
				n_cells = knn.shape[0]
				a = np.tile(np.arange(n_cells), k)
				b = np.reshape(knn.T, (n_cells * k,))
				w = np.repeat(1 / np.power(np.arange(1, k + 1), 1.8), n_cells)
				knn = sparse.coo_matrix((w, (a, b)), shape=(n_cells, n_cells))
				threshold = w > 0.025
				mknn = sparse.coo_matrix((w[threshold], (a[threshold], b[threshold])), shape=(n_cells, n_cells))
				mknn = mknn.minimum(mknn.transpose()).tocoo()
				tsne = cg.TSNE(perplexity=5).layout(transformed)
				dsagg.col_graphs.KNN = knn
				dsagg.col_graphs.MKNN = mknn
				dsagg.ca._X = tsne[:, 0]
				dsagg.ca._Y = tsne[:, 1]
				
				logging.info("Manifold learning on all cells")
				init = np.zeros((ds.shape[1], 2))
				for lbl in np.unique(ds.ca.Clusters):
					init[ds.ca.Clusters == lbl, :] = tsne[lbl, :] + np.random.normal(size=((ds.ca.Clusters == lbl).sum(), 2))
				ml = cg.ManifoldLearning2(gtsne=True, alpha=1, max_iter=3000)
				(knn, mknn, tsne) = ml.fit(ds, initial_pos=init, nng=nng, blocked_genes=blocked)
				ds.col_graphs.KNN = knn
				ds.col_graphs.MKNN = mknn
				ds.ca._X = tsne[:, 0]
				ds.ca._Y = tsne[:, 1]

				logging.info("Computing auto-annotation")
				aa = cg.AutoAnnotator(root="../auto-annotation/Adolescent/")
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


