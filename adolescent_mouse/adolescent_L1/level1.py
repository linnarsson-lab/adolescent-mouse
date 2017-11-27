from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am


class Level1(luigi.WrapperTask):
	"""
	Luigi Task to run all analyses
	"""

	def requires(self) -> Iterator[luigi.Task]:
		tissues = cg.PoolSpec().tissues_for_project("Adolescent")
		for tissue in tissues:
			yield am.ExportL1(tissue=tissue)
