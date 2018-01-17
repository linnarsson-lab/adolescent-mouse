from typing import *
import os
import csv
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import numpy_groupies.aggregate_numpy as npg
import scipy.stats
import adolescent_mouse as am
import pandas as pd
from shutil import copyfile
from sklearn.preprocessing import LabelEncoder


class ExtractByTaxonL6(luigi.Task):
	"""
	Extract cells by taxon for level 6
	"""
	rank = luigi.IntParameter()

	def requires(self) -> luigi.Task:
		return am.PoolL5()

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"L6_R{self.rank}.loom"))
		
	def run(self) -> None:
		logging = cg.logging(self)
		with self.output().temporary_path() as out_file:
			copyfile(self.input().fn, out_file)
			with loompy.connect(out_file) as ds:
				labels = ds.ca[f"TaxonomyRank{self.rank}"]
				le = LabelEncoder()
				new_clusters = le.fit_transform(labels)
				ds.ca.Clusters = new_clusters

