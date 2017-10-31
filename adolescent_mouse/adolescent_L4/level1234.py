from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am


class Level1234(luigi.WrapperTask):
	"""
	Luigi Task to run all Level 1-4 analyses
	"""

	def requires(self) -> Iterator[luigi.Task]:
		targets = [
			'Sensory_Neurons',
			'Sympathetic_Neurons',
			'Enteric_Neurons',
			'DiMesencephalon_Excitatory',
			'Hindbrain_Inhibitory',
			'SpinalCord_Inhibitory',
			'Brain_Granule',
			'Brain_CholinergicMonoaminergic',
			'DiMesencephalon_Inhibitory',
			'Striatum_MSN',
			'Hypothalamus_Peptidergic',
			'Forebrain_Excitatory',
			'Brain_Neuroblasts',
			'Hindbrain_Excitatory',
			'SpinalCord_Excitatory',
			'Forebrain_Inhibitory',
			'Olfactory_Inhibitory'
		]
		tissues = cg.PoolSpec().tissues_for_project("Adolescent")
		classes = ["Oligos", "Ependymal", "Astrocytes", "Vascular", "Immune", "PeripheralGlia"]
		for tissue in tissues:
			yield am.ExportL1(tissue=tissue)
			yield am.ExportL2(tissue=tissue, major_class="Neurons")

		for cls in classes:
			yield am.ExportL2(tissue="All", major_class=cls)

		for target in targets:
			yield am.ExportL3(target=target)

		yield am.ExportL4()
