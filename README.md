# adolescent-mouse

Analysis pipeline for the adolescent mouse nervous system project

## Installation

1. Install [luigi](https://luigi.readthedocs.io/en/stable/), [cytograph](https://github.com/linnarsson-lab/cytograph) and [loompy](https://github.com/linnarsson-lab/loompy)

2. Clone the repository to your computer:

```
git clone https://github.com/linnarsson-lab/adolescent-mouse.git
```

3. Install in development mode using `pip`:

```
cd adolescent-mouse
pip install -e .
```

(note the trailing dot)

## Setting up your environment

In order to run the pipeline, you will need to set up your samples folder:

* A folder containing the raw samples (loom files named like `10X23_1.loom`)
* Inside it, a sub-folder named `classified` and containing pre-classified clusters (loom files named like `L0_Cortex1.loom`)
* Also inside it, a sub-folder named `metadata` and containing a file named `metadata.xlsx` with sample metadata

For example, you could have something like this:

```
loom_samples/
    10X23_1.loom
    10X23_2.loom
    10X53_1.loom
    ...
    classified/
        L0_Cortex1.loom
        L0_Hippocampus.loom
        ...
        classifier.pickle
    metadata/
        metadata.xlsx
```

**Note:** On monod, a samples folder is available at `/data/proj/chromium/`. Please use this directly instead of making a copy.

**Note:** Instead of sample files directly in the top-level folder, sample files can also be in subfolders named after the sample ID, or in a sub-folder of such a subfolder, named `velocyto`. Thus, these are all equivalent:

```
10X43_2.loom
10X43_2/10X43_2.loom
10X43_2/velocyto/10X43_2.loom
```

If it exists, the sample in `velocyto` will take precedence.

Furthermore, you need a file in the current directory (typically, `adolescent-mouse`) named `pooling_specification.tab`, which gives a list of all the samples and their pool names. This file has four columns: *SampleID* (like `10X04_1`), *Pool* (like `Hippocampus`), *TimepointPool* (always `none`), *QC* (`OK` or `FAILED`), and *Project* (`Adolescent`). 

Samples with `QC == FAILED` will be ignored for all analyses.

## Running the pipeline (example: level 2 astrocytes)

1. Create a folder to hold the output of the build, e.g. `~/build_20171027`.

2. Run `luigi`. For example:

```
luigi --local-scheduler --module adolescent_mouse ExportL2 --major-class Astrocytes --tissue All --paths-samples /data/proj/chromium/ --paths-build ~/build_20171027/ 
```

This will run the `ExportL2` task, and all of its dependencies as needed. The command-line arguments are as follows:

Argument|Effect
----|----
`--local-scheduler` |run the pipeline directly, not in client/server mode.
`--module adolescent_mouse ExportL2`| run the `ExportL2` task in the Python module `adolescent_mouse`.
`--major-class Astrocytes`| an argument to `ExportL2` indicating which cell class we want to run
`--tissue All` |an argument to `ExportL2` indicating that we want astrocytes from all tissues
`--paths-samples /data/proj/chromium/loom_samples/` | a configuration of the `paths` object, setting the sample path
`--paths-build ~/build_20171027/` | a configuration of the `paths` object, setting the build path

**Note**: You can add the `--workers 20` to use 20 cores (useful when running on the server). This will significantly speed up the pipeline, because luigi will automatically figure out which tasks can be run in parallel.


## Running everything

The following command runs all levels of the analysis and exports the results at every level:

```
luigi --local-scheduler --module adolescent_mouse Level1234Adolescent --paths-samples /data/proj/chromium/l --paths-build ~/build_20171027/ 
```

## Rerunning parts of the pipeline

To rerun only a subset of the pipeline, simply delete the files (or folders) you want regenerated, and rerun the pipeline. Luigi will figure out what needs to be done. 

Note, however, that luigi will start with the module you ask it to run, and work its way backward. 

For example, if you delete the `L1_Hippocampus.loom` file, but leave `L2_Neurons_Hippocampus.loom`, and ask luigi to run `ExportL2(tissue=Hippocampus,major_class=Neurons)`, then the `L1_Hippocampus.loom` file will not be generated because the result you asked for already exists. You can then either delete `L2_Neurons_Hippocampus.loom` or run `ExportL1(tissue=Hippocampus)` to force luigi to regenerate the file.

## Running level 0 analysis

Level 0 consists in pooling raw samples by tissue according to the `pooling_specification.tab` file, and classifying major cell classes using the `classifier.pickle` pre-trained classifier. Optionally, you can train a new classifier.

Task(args)|Purpose
----|-----
`TrainClassifier`|Train a new classifier, and create a `classifier.pickle` file
`PrepareTissuePool(tissue)`|Prepare a pool of the samples corresponding to `tissue`

You will normally not need to run level 0 alone, but it will be triggered by running a level 1 analysis.

## Running level 1 analysis

Level 1 performs manifold learning, clustering and annotation by tissue.

Task(args)|Purpose|Output|Depends on
----|-----|----|----
`ExportL1(tissue)`| Export the results | `L1_{tissue}_exported` | `ClusterL1`, `AggregateL1`
`AggregateL1(tissue)`| Aggregate by cluster, computing enrichment, trinarization and auto-annotation | `L1_{tissue}.agg.loom` | `ClusterL1`
`ClusterL1(tissue)`| Manifold learning and clustering | `L1_{tissue}.loom` | `PrepareTissuePool`

Example:

```
luigi --local-scheduler --module adolescent_mouse ExportL1 --tissue Hippocampus --paths-samples /data/proj/chromium/ --paths-build ~/build_20171027/ 
```


## Running level 2 analysis

Level 2 performs manifold learning, clustering and annotation by major class and (for neurons only) tissue. The major classes (as determined by the level 0 classifier) are: `Astrocytes`, `Ependymal`, `Vascular`, `Immune`, `Oligos`, `PeripheralGlia` and `Neurons`.

Task(args)|Purpose|Output|Depends on
----|-----|----|----
`ExportL2(major_class, tissue)`| Export the results | `L2_{major_class}_{tissue}_exported` | `FilterL2`, `AggregateL2`
`FilterL2(major_class, tissue)`| Filter to remove bad clusters | `L2_{major_class}_{tissue}.filtered.loom` | `ClusterL2`
`AggregateL2(major_class, tissue)`| Aggregate by cluster, computing enrichment, trinarization and auto-annotation | `L2_{major_class}_{tissue}.agg.loom` | `FilterL2`
`ClusterL2(major_class, tissue)`| Manifold learning and clustering | `L2_{major_class}_{tissue}.loom` | `ClusterL1(tissue)`

**Note:** For classes other than `Neurons`, use `All` as tissue name.

Example:

```
luigi --local-scheduler --module adolescent_mouse ExportL2 --major-class Neurons --tissue Hippocampus --paths-samples /data/proj/chromium/ --paths-build ~/build_20171027/ 
```


## Running level 3 analysis

Level 3 performs manifold learning, clustering and annotation (of neurons only) by regional and neurotransmitter target categories. The targets are: 
```    
	SpinalCord_Inhibitory
	SpinalCord_Excitatory
	Olfactory_Inhibitory
	Enteric_Neurons
	Sensory_Neurons
	Sympathetic_Neurons	
	Hypothalamus_Peptidergic
	Hindbrain_Inhibitory
	Hindbrain_Excitatory
	Brain_Neuroblasts
	Forebrain_Inhibitory
	Forebrain_Excitatory
	DiMesencephalon_Inhibitory
	DiMesencephalon_Excitatory
	Brain_Granule
	Brain_CholinergicMonoaminergic
	Striatum_MSN
```
Task(args)|Purpose|Output|Depends on
----|-----|----|----
`ExportL3(target)`| Export the results | `L3_{target}_exported` | `ClusterL3`, `AggregateL3`
`AggregateL3(target)`| Aggregate by cluster, computing enrichment, trinarization and auto-annotation | `L3_{target}.agg.loom` | `ClusterL3`
`ClusterL3(target)`| Manifold learning and clustering | `L3_{target}.loom` | `ClusterL2`


Example:

```
luigi --local-scheduler --module adolescent_mouse ExportL3 --target Sympathetic_Neurons --paths-samples /data/proj/chromium/ --paths-build ~/build_20171027/ 
```




## Running level 4 analysis

Level 4 pools all cells into a single file, and performs aggregation and visualization. 

Task(args)|Purpose|Output|Depends on
----|-----|----|----
`ExportL4`| Export the results | `L4_All_exported` | `PoolL4`, `AggregateL4`
`AggregateL4`| Aggregate by cluster, computing enrichment, trinarization and auto-annotation | `L4_All.agg.loom` | `PoolL4`
`PoolL4`| Pool all level 3 (neurons) and level 2 (non-neurons) cells | `L4_All.loom` | `ClusterL3, FilterL2`


Example:

```
luigi --local-scheduler --module adolescent_mouse ExportL4 --paths-samples /data/proj/chromium/ --paths-build ~/build_20171027/ 
```



