from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am


class OptimizeC2(luigi.WrapperTask):
	"""
	Luigi Task to run all analyses
	"""
	tissue = luigi.Parameter()

	def requires(self) -> Iterator[luigi.Task]:
		for (a, b) in [(0.3, 0.3), (0.3, 5), (0.3, 10), (1, 10)]:
			for (c, d) in [(0.3, 0.3), (0.3, 5), (0.3, 10), (1, 10)]:
				for accel in [True, False]:
					for log in [True, False]:
						normalize = log
						for k_smoothing in [10]:
							for k in [25, 100]:
								for n_factors in [100]:
									yield am.ExportL1C2(tissue=self.tissue, a=a, b=b, c=c, d=d, 
										n_factors=n_factors, k_smoothing=k_smoothing, k=k, 
										log=log, normalize=normalize, accel=accel)
