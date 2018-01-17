from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am


class Level4(luigi.WrapperTask):
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
		classes = ["Oligos", "Ependymal", "Astrocytes", "Vascular", "Immune", "PeripheralGlia"]

		for cls in classes:
			yield am.ExportL4(target=cls)

		for target in targets:
			yield am.ExportL4(target=target)
