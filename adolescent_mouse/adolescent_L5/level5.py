from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am


class Level5(luigi.WrapperTask):
	"""
	Luigi Task to run level 5
	"""

	def requires(self) -> Iterator[luigi.Task]:
		yield am.ExportL5()
