from setuptools import setup, find_packages

__version__ = "0.0.0"
exec(open('adolescent_mouse/_version.py').read())

setup(
	name="adolescent-mouse",
	version=__version__,
	packages=find_packages(),
	install_requires=[
		'loompy',
		'numpy',
		'scikit-learn',
		'scipy',
		'networkx',
		'cytograph',
		'luigi',
		'polo',
		'numpy-groupies'
	],
	
	# metadata
	author="Linnarsson Lab",
	author_email="sten.linnarsson@ki.se",
	description="Pipeline for single-cell RNA-seq analysis",
	license="MIT",
	url="https://github.com/linnarsson-lab/cytograph",
)
