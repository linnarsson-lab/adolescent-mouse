from typing import *
import os
import csv
import logging
import pickle
import loompy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cytograph as cg
import luigi
import numpy_groupies.aggregate_numpy as npg
import scipy.stats
import adolescent_mouse as am
import pandas as pd
from sklearn.preprocessing import LabelEncoder


class OligosL6(luigi.Task):
	"""
	Extract and plot oligos cell types
	"""

	def requires(self) -> luigi.Task:
		return am.PoolL5()

	def output(self) -> luigi.Target:
		return luigi.LocalTarget(os.path.join(am.paths().build, f"F_Oligos"))
		
	def run(self) -> None:
		logging = cg.logging(self, True)
		with self.output().temporary_path() as out_dir:
			logging.info("Exporting oligo cell types")
			if not os.path.exists(out_dir):
				os.mkdir(out_dir)
			with loompy.connect(self.input().fn) as ds:
				celltypes = [
					"COP1",
					"COP2",
					"NFOL2",
					"NFOL1",
					"OPC"
				]
				selected = np.array([], dtype='int')
				for ct in celltypes:
					print(ct)
					cells = np.where(ds.ca.ClusterName == ct)[0]
					if cells.shape[0] > 820:
						cells = np.random.choice(cells, size=820, replace=False)
					selected = np.union1d(selected, cells)

				ngfile = os.path.join(out_dir, "F_Oligos.loom")
				for (_, _, view) in ds.scan(items=selected, axis=1):
					loompy.create_append(ngfile, view.layers, view.ra, view.ca)
			
			with loompy.connect(ngfile) as ds:
				logging.info("Learning the manifold")
				ml = cg.ManifoldLearning2(gtsne=False, alpha=1)
				(knn, mknn, tsne) = ml.fit(ds)
				ds.col_graphs.KNN = knn
				ds.col_graphs.MKNN = mknn
				ds.ca._X = tsne[:, 0]
				ds.ca._Y = tsne[:, 1]

				fig = plt.figure(figsize=(3, 3))
				ax = fig.add_axes([0, 0, 1, 1])
				lc = LineCollection(zip(tsne[mknn.row], tsne[mknn.col]), linewidths=0.25, zorder=0, color='grey', alpha=0.1)
				ax.add_collection(lc)
				ax.axis('off')
				colors = cg.colorize(np.unique(ds.ca.ClusterName))
				ix = 0
				for ct in np.unique(ds.ca.ClusterName):
					cells = (ds.ca.ClusterName == ct)
					plt.scatter(x=ds.ca._X[cells], y=ds.ca._Y[cells], s=40, c=colors[ix, :], marker='.', label=ct, alpha=0.5, lw=0)
					ix += 1
					lgnd = ax.legend(fontsize=10, labelspacing=0.2, loc="upper left", bbox_to_anchor=(1, 1), frameon=False)
					for handle in lgnd.legendHandles:
						handle.set_sizes([250])
						handle.set_alpha(1)
				plt.savefig(os.path.join(out_dir, "Fig_Oligos_Types.png"), dpi=600, transparent=True, bbox_extra_artists=(lgnd,), bbox_inches='tight')
				plt.close()

				fig = plt.figure(figsize=(3, 3))
				ax = fig.add_axes([0, 0, 1, 1])
				ax.axis('off')
				plt.scatter(x=ds.ca._X, y=ds.ca._Y, s=40, c=cg.colors75[(ds[ds.ra.Gene == "Cdk1", :][0] != 0).astype('int')], marker='.', label=ct, alpha=0.5, lw=0)
				plt.savefig(os.path.join(out_dir, "Fig_Oligos_Cdk1.png"), dpi=600, transparent=True, bbox_inches='tight')
