# MethylationToActivity
## a deep-learning framework that reveals promoter activity landscapes from DNA methylomes in individual tumors

MethylationToActivity (M2A) is a machine learning framework using convolutional neural networks (CNN) to infer histone modification (HM) enrichment from whole genome bisulfite sequencing (WGBS). To date, both H3K27ac and H3K4me3 enrichment prediction from WGBS is supported, from a tab-delimited text file format of M-values. Optionally, we also support transfer-learning where a user may have matching H3K27ac or H3K4me3 data with appropriate controls in addition to WGBS data.


### M2A is comprised 5 parts, including transfer-learning:
|Process                | Description                                                                                                        |
|-----------------------|--------------------------------------------------------------------------------------------------------------------|
|1_ResponseVariable     | Generate histone enrichment for each unique promote region (transfer-learning only)                                |
|2_MethylationFeatures  | Process WGBS features for model input                                                                              |
|3_CombineInput         | Scale and recombine features, and for transfer learning, calculated HM values                                      |
|4_TransferLearning     | Train fully-connected layers of a particular model for increased performance in your domain of interest (Optional) |
|5_RunModel             | Using pre-generated input, get HM predictions for each unique promoter region                                      |

### Prerequisites

Python 3.6.5 or greater:
1) [pyBigWig] v0.3.13
2) [numpy] v1.17.1
3) [pandas] v0.25.1
4) [pandarallel] v1.4.2
5) [scikit-learn] 0.20.2
6) [h5py] v2.9.0
7) [keras] v2.2.4
8) [tensorflow] v1.10.1
9) [scipy] v1.3.1
10) [matplotlib] v3.3.0
11) [cwltool] v1.0
12) [psutil] v5.6.1

[Python]: https://www.python.org/
[cwltool]: https://github.com/common-workflow-language/cwltool
[numpy]: https://numpy.org/
[pandas]: https://pandas.pydata.org/
[h5py]: https://www.h5py.org/
[keras]: https://keras.io/
[scikit-learn]: https://scikit-learn.org/
[psutil]: https://pypi.org/project/psutil/
[scipy]: https://www.scipy.org/
[tensorflow]: https://www.tensorflow.org/
[matplotlib]: https://matplotlib.org/
[pyBigWig]: https://github.com/deeptools/pyBigWig
[pandarallel]: https://github.com/nalepae/pandarallel

### Obtain M2A
Clone M2A from GitHub: 
```
git clone https://github.com/chenlab-sj/M2A.git
```

### Inputs
M2A requires five inputs, defined in a YAML file as [CWL inputs]. E.g., `inputs.yml`:
```
chipBigwig:
  class: File
  path: sample.bw
inputBigwig:
  class: File
  path: input.bw
curated:
  class: File
  path: sites.txt
promoterDefinitions:
  class: File
  path: promoters.txt
model: 
  class: File
  path: model.h5
```
[CWL inputs]: https://www.commonwl.org/user_guide/02-1st-example/index.html
[Gencode]: https://www.gencodegenes.org/

#### Input description
| Name | Description |
|-------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| Sample HM bigwig file (only if using M2A with Transfer)           | HM ChIP-seq experiment bigwig track.                                              |
| Sample HM control (Input) bigwig (only if using M2A with Transfer)| ChIP-seq Experiment control (Input) bigwig track.                                 |
| WGBS data file                                                    | M-values by chromosome and position (non-standard format, see below).             |
| Promoter region definition file (*provided, or user defined*)     | File describing promoter regions to be predicted. (non-standard format, see below)|
| Model weights (*provided, or user defined from transfer*)         | hdf5 model weights for either H3K27ac prediction OR H3K4me3 prediction            |

##### Promoter region definition file 
A tab delimited file containing the unique promoter-regions for either:
* **hg19-based data: 2_Promoter_Definitions_hg19.txt, or**
* **GRCh38-based data: 2_Promoter_Definitions_GRCh38.txt**

|Column      | Description                                   |
|------------|-----------------------------------------------|
|EnsmblID_T  | Ensemble transcript ID (unique)               |
|EnsmblID_G  | Ensemble gene ID (not unique)                 |      
|Gene        | human readable gene name (abbrev, not unique) |
|Strand      | +, -                                          |
|Chr         | chr1, chr2, ... chr22, etc.                   |
|Start       | Beginning of transcript definition            | 
|End         | End of transcript definition                  |
|RStart      | TSS - 1000bp                                  |
|REnd        | TSS + 1000bp                                  |

##### WGBS data file
A bed-like file of genomic positions with corresponding M-values, tab delimited:

|Column      | Description                                                            |
|------------|------------------------------------------------------------------------|
|chrom       | chromosome ID, e.g. 1,2,3 ...22                                        |
|pos         | position of 5' cytosine of a CpG on the positive strand                |      
|mval        | caluclated mvalue of a given CpG, typically M-value=log2(Beta/1-Beta)  |


### Run M2A with transfer learning 

M2A uses [CWL] to describe its workflow. To run an example workflow, update `sample_data/input_data/inputs.yml` with the path to a promoter definitions file.
Then run the following.

```
$ mkdir results
$ cwltool --outdir results cwl/m2a.cwl sample_data/input_data/inputs.yml
```
### Run M2A without transfer learning 

M2A without transfer learning enabled is contained in the CWL workflow `cwl/m2a_without_transfer_learning.cwl`. It requires the same inputs as the with transfer learning pipeline, with the exception of the bigwig files.

[CWL]: https://www.commonwl.org/

## Docker

M2A provides a [Dockerfile] that builds an image with all the included dependencies. To use this image, install [Docker] for your platform. This Docker image is used by the CWL workflow and contains the prerequisites.

[Dockerfile]: https://docs.docker.com/engine/reference/builder/
[Docker]: https://www.docker.com/

### Build Docker image

In the M2A project directory, build the Docker image.

```
$ docker build --tag stjude/m2a:0.0.1 .
```

## Evaluate test data results

Today, the M2A pipeline does not produce an interactive visualization. If M2A with Transfer was run, the easiest measurment of training prediction accuracy would be caluclating the Pearson's R<sup>2</sup>, or root mean square error (RMSE) between the measured and M2A predicted values. Furthermore, comparisons of sample-sample consistency with the same/similar cancer-type (as determiend by Pearson's R<sup>2</sup>) is a good start for a contextual understanding of the predictions produced by M2A.

## St. Jude Cloud

To run M2A in St. Jude Cloud, please follow the directions at https://www.stjude.cloud/docs/guides/tools/methylation-to-activity/

## Availability

Copyright 2019 St. Jude Children's Research Hospital

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## Seeking help

For questions and bug reports, please open an issue on the [GitHub] project page.

[GitHub]: https://github.com/stjude/M2A/issues

## Publication analyses

All scripts describing the experiments and analyses in the M2A publication, including previous (unsupported) versions of M2A, can be found in the M2A_analyses directory.

## Citing M2A

(In submission) Justin Williams, Beisi Xu, Daniel Putnam, Andrew Thrasher, and Xiang Chen. MethylationToActivity: a deep-learning framework that reveals promoter activity landscapes from DNA methylomes in individual tumors.
