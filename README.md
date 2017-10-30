# adolescent-mouse

Analysis pipeline for the adolescent mouse nervous system project

## Installation

1. Install [cytograph](https://github.com/linnarsson-lab/cytograph) and [loompy](https://github.com/linnarsson-lab/loompy)

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

**Note:** On monod, a samples folder is available at `/data/proj/chromium/loom_samples`. Please use this directly instead of making a copy.

Furthermore, you need a file in the current directory (typically, `adolescent-mouse`) named `pooling_specification.tab`, which gives a list of all the samples and their pool names. This file has four columns: *SampleID* (like `10X04_1`), *Pool* (like `Hippocampus`), *TimepointPool* (always `none`), *QC* (`OK` or `FAILED`), and *Project* (`Adolescent`). Samples with `QC == FAILED` will be ignored.




