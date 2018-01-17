from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am


class Level123(luigi.WrapperTask):
	"""
	Luigi Task to run all analyses
	"""

	def requires(self) -> Iterator[luigi.Task]:
		targets = [
			'Sensory_Neurons',
			'Sympathetic_Neurons',
			'Enteric_Neurons',
			'Mesencephalon_Excitatory',
			'Diencephalon_Excitatory',
			'Hindbrain_Inhibitory',
			'SpinalCord_Inhibitory',
			'Brain_Granule',
			'Brain_CholinergicMonoaminergic',
			'Mesencephalon_Inhibitory',
			'Diencephalon_Inhibitory',
			'Striatum_MSN',
			'Hypothalamus_Peptidergic',
			'Telencephalon_Excitatory',
			'Brain_Neuroblasts',
			'Hindbrain_Excitatory',
			'SpinalCord_Excitatory',
			'Telencephalon_Inhibitory',
			'Olfactory_Inhibitory'
		]
		tissues = cg.PoolSpec().tissues_for_project("Adolescent")
		classes = ["Oligos", "Ependymal", "Astrocytes", "Vascular", "Immune", "PeripheralGlia"]
		for tissue in tissues:
			yield am.ExportL1(tissue=tissue)
			yield am.ExportL2(tissue=tissue, major_class="Neurons")

		for cls in classes:
			yield am.ExportL2(tissue="All", major_class=cls)
			yield am.ExportL3(target=cls)

		for target in targets:
			yield am.ExportL3(target=target)
