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


class AggregateL2(luigi.Task):
	"""
	Aggregate all clusters in a new Loom file
	"""
	major_class = luigi.Parameter()
	tissue = luigi.Parameter(default="All")
	n_markers = luigi.IntParameter(default=10)
	n_auto_genes = luigi.IntParameter(default=6)

	def requires(self) -> List[luigi.Task]:
		return am.ClusterL2(tissue=self.tissue, major_class=self.major_class)

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, "L2_" + self.major_class + "_" + self.tissue + ".agg.loom"))

	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			logging.info("Aggregating loom file")
			with loompy.connect(self.input().fn) as ds:
				cg.Aggregator(self.n_markers).aggregate(ds, out_file)
				with loompy.connect(out_file) as dsagg:
					for ix, score in enumerate(dsagg.col_attrs["ClusterScore"]):
						logging.info(f"Cluster {ix} score {score:.1f}")

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

					tissue = self.tissue
					labels = ds.col_attrs["Clusters"]

					# Figure out which cells should be collected
					cells: List[int] = []
					# clusters_seen: List[int] = []  # Clusters for which there was some schedule
					clusters_seen: Dict[int, str] = {}
					schedule = pooling_schedule_L3[self.tissue]

					# Where to send clusters when no rules match
					_default_schedule: str = None
					for aa_tag, sendto in schedule:
						if aa_tag == "*":
							_default_schedule = sendto

					# For each cluster in the tissue
					bucket_list = []
					for ix, agg_aa in enumerate(dsagg.ca.AutoAnnotation):
						# For each rule in the schedule
						for aa_tag, sendto in schedule:
							if aa_tag in agg_aa.split(","):
								if ix in clusters_seen:
									logging.info(f"{tissue}/{ix}/{agg_aa}: {aa_tag} -> {sendto} (overruled by '{clusters_seen[ix]}')")
								else:
									clusters_seen[ix] = f"{aa_tag} -> {sendto}"
									logging.info(f"{tissue}/{ix}/{agg_aa}: {aa_tag} -> {sendto}")
									bucket_list.append(sendto)
						if ix not in clusters_seen:
							if _default_schedule is None:
								logging.info(f"{tissue}/{ix}/{agg_aa}: No matching rule")
								bucket_list.append("Excluded")
							else:
								clusters_seen[ix] = f"{aa_tag} -> {_default_schedule}"
								logging.info(f"{tissue}/{ix}/{agg_aa}: {aa_tag} -> {_default_schedule}")
								bucket_list.append(_default_schedule)
					dsagg.ca.Bucket = np.array(bucket_list)


pooling_schedule_L3 = {
	"Cortex1": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"Cortex2": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"Cortex3": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"Hippocampus": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@NBL", "Brain_Neuroblasts"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
	],
	"StriatumDorsal": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"StriatumVentral": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"Amygdala": [
		("@GABAGLUT1", "Amygdala_Other"),
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "Forebrain_Inhibitory"),
		("DG-GC", "Brain_Granule"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"Olfactory": [
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@NIPC", "Brain_Neuroblasts"),
		("@VGLUT1", "Forebrain_Excitatory"),
		("@VGLUT2", "Forebrain_Excitatory"),
		("@VGLUT3", "Forebrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts"),
		("*", "Olfactory_Inhibitory")
	],
	"Hypothalamus": [
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@OXT", "Hypothalamus_Peptidergic"),
		("@AVP", "Hypothalamus_Peptidergic"),
		("@GNRH", "Hypothalamus_Peptidergic"),
		("@AGRP", "Hypothalamus_Peptidergic"),
		("@HCRT", "Hypothalamus_Peptidergic"),
		("@PMCH", "Hypothalamus_Peptidergic"),
		("@POMC", "Hypothalamus_Peptidergic"),
		("@TRH", "Hypothalamus_Peptidergic"),
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@VGLUT1", "DiMesencephalon_Excitatory"),
		("@VGLUT2", "DiMesencephalon_Excitatory"),
		("@VGLUT3", "DiMesencephalon_Excitatory"),
		("@GABA", "DiMesencephalon_Inhibitory")
	],
	"MidbrainDorsal": [
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@SER", "Brain_CholinergicMonoaminergic"),
		("@DA", "Brain_CholinergicMonoaminergic"),
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "DiMesencephalon_Inhibitory"),
		("@VGLUT1", "DiMesencephalon_Excitatory"),
		("@VGLUT2", "DiMesencephalon_Excitatory"),
		("@VGLUT3", "DiMesencephalon_Excitatory")
	],
	"MidbrainVentral": [
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@SER", "Brain_CholinergicMonoaminergic"),
		("@DA", "Brain_CholinergicMonoaminergic"),
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "DiMesencephalon_Inhibitory"),
		("@VGLUT1", "DiMesencephalon_Excitatory"),
		("@VGLUT2", "DiMesencephalon_Excitatory"),
		("@VGLUT3", "DiMesencephalon_Excitatory")
	],
	"Thalamus": [
		("MSN-D1", "Striatum_MSN"),
		("MSN-D2", "Striatum_MSN"),
		("@SER", "Brain_CholinergicMonoaminergic"),
		("@DA", "Brain_CholinergicMonoaminergic"),
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GABA", "DiMesencephalon_Inhibitory"),
		("@VGLUT1", "DiMesencephalon_Excitatory"),
		("@VGLUT2", "DiMesencephalon_Excitatory"),
		("@VGLUT3", "DiMesencephalon_Excitatory")
	],
	"Cerebellum": [
		("@NIPC", "Brain_Neuroblasts"),
		("CB-PC", "Hindbrain_Inhibitory"),
		("CB-GC", "Brain_Granule"),
		("@GLY", "Hindbrain_Inhibitory"),
		("@GABA", "Hindbrain_Inhibitory"),
		("@VGLUT1", "Hindbrain_Excitatory"),
		("@VGLUT2", "Hindbrain_Excitatory"),
		("@VGLUT3", "Hindbrain_Excitatory"),
		("@NBL", "Brain_Neuroblasts")
	],
	"Pons": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@SER", "Brain_CholinergicMonoaminergic"),
		("@DA", "Brain_CholinergicMonoaminergic"),
		("@NOR", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@GLY", "Hindbrain_Inhibitory"),
		("@VGLUT1", "Hindbrain_Excitatory"),
		("@VGLUT2", "Hindbrain_Excitatory"),
		("@VGLUT3", "Hindbrain_Excitatory")
	],
	"Medulla": [
		("@CHOL", "Brain_CholinergicMonoaminergic"),
		("@SER", "Brain_CholinergicMonoaminergic"),
		("@DA", "Brain_CholinergicMonoaminergic"),
		("@NOR", "Brain_CholinergicMonoaminergic"),
		("@ADR", "Brain_CholinergicMonoaminergic"),
		("@NIPC", "Brain_Neuroblasts"),
		("@VGLUT1", "Hindbrain_Excitatory"),
		("@VGLUT2", "Hindbrain_Excitatory"),
		("@VGLUT3", "Hindbrain_Excitatory"),
		("@GLY", "Hindbrain_Inhibitory")
	],
	"SpinalCord": [
		("@NIPC", "Brain_Neuroblasts"),
		("@VGLUT1", "SpinalCord_Excitatory"),
		("@VGLUT2", "SpinalCord_Excitatory"),
		("@VGLUT3", "SpinalCord_Excitatory"),
		("@GLY", "SpinalCord_Inhibitory"),
		("@GABA", "SpinalCord_Inhibitory"),
		("PSN", "Exclude"),
		("*", "SpinalCord_Excitatory"),
	],
	"DRG": [
		("*", "Sensory_Neurons")
	],
	"Sympathetic": [
		("*", "Sympathetic_Neurons")
	],
	"Enteric": [
		("*", "Enteric_Neurons")
	],
}