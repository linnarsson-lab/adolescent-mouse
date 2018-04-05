from typing import *
import os
import logging
import pickle
import loompy
import numpy as np
import cytograph as cg
import luigi
import adolescent_mouse as am
import pandas as pd


class Level6(luigi.WrapperTask):
	"""
	Luigi Task to run level 6
	"""

	def requires(self) -> Iterator[luigi.Task]:
		taxonomy_file = os.path.join(am.paths().build, "curated_L4", "Taxonomy.xlsx")
		taxonomy_table = pd.read_excel(taxonomy_file)
		taxonomy = {taxonomy_table.columns.values[i]: taxonomy_table.values[:, i] for i in range(taxonomy_table.shape[1])}
		for taxon in list(set(taxonomy["TaxonomyRank1"])):
			yield am.ExportL6(rank=1, taxon=taxon)
		for taxon in list(set(taxonomy["TaxonomyRank2"])):
			yield am.ExportL6(rank=2, taxon=taxon)
		for taxon in list(set(taxonomy["TaxonomyRank3"])):
			yield am.ExportL6(rank=3, taxon=taxon)
#		for taxon in list(set(taxonomy["TaxonomyRank4"])):
#			yield am.ExportL6(rank=4, taxon=taxon)
		for rank in [1, 2, 3, 4]:
			yield am.ExportByTaxonL6(rank=rank)
		yield am.ExportL5()
