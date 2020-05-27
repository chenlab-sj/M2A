# MethylationToActivity: a deep-learning framework that reveals promoter activity landscapes from DNA methylomes in individual tumors

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
```
Linux Operating System
Python 3.6.5 or greater

the following non-standard python packages, available through pip3 install [package]:
1) pyBigWig v0.3.13
2) argparse v1.1
3) numpy v1.17.1
4) pandas v0.25.1
5) pandarallel v1.4.2
6) sklearn 0.20.2+
7) h5py v2.9.0
8) keras v2.2.4
9) tensorflow v1.10.1
10) scipy v1.3.1
for more details, please visit individual scripts
```

## 1_ResponseVariable

Calculate the response variable (histone enrichment) for each promoter region, for transfer learning.

#### Mandatory Arguments
|Argument            | Description                                  |
|--------------------|----------------------------------------------|
|ChIP_Path           | Sample HM bigwig file path                   |
|Input_Path          | Sample HM control (Input) bigwig file path   |
|PromoterDefinitions | Path to Promoter Definitions file            |

#### Optional Arguments
|Argument            | Description                                  |
|--------------------|----------------------------------------------|
|--outFileName       | File name to use for output                  |
|--outDirectory      | Directory to write output to                 |


#### Commandline Usage
```
python 1_getResponseVariable.py [ChIP_Path] [Input_Path] [PromoterDefinitions]
```

## 2_MethylationFeatures

Extract methylation features on a per window, promoter (region) basis.

#### Mandatory Arguments
|Argument            | Description                                  |
|--------------------|----------------------------------------------|
|MethylFilePath      | Path to methylation bed-like file            |
|PromoterDefinitions | Path to Promoter Definitions file            |

#### Optional Arguments
|Argument            | Description                                  |
|--------------------|----------------------------------------------|
|--nbWorkers         | No. of threads to use, default 2             |
|--outFileName       | File name to use for output                  |
|--outDirectory      | Directory to write output to                 |

#### Commandline Usage
```
python 2_getMethylation.py [MethylFilePath] [PromoterDefinitions] 
```

## 3_CombineInput

Combine/interleave all features into pseudo image arrays for input to a CNN model.

#### Mandatory Arguments
|Argument            | Description                                                                     |
|--------------------|---------------------------------------------------------------------------------|
|MethylationFilePath | Full path to methylation features file (output from step 2_MethylationFeatures) |

#### Optional Arguments
|Argument               | Description                                                              |
|-----------------------|--------------------------------------------------------------------------|
|--ResponseVariablePath | Full path to response variable file, only required for transfer learning |
|--outFileName          | File name to use for output                                              |
|--outDirectory         | Directory to write output to                                             |

#### Commandline Usage
```
python 3_combineInput.py [MethylationFilePath]
```

## 4_TransferLearning

Perform transfer learning on a single sample, updating a current model.

#### Mandatory Arguments
|Argument            | Description                                                                     |
|--------------------|---------------------------------------------------------------------------------|
|FeatureFilePath     | Full path to methylation features .h5 file (output from step 3_CombineInput)    |
|ModelFilePath       | Full path to model .h5 file (Provided H3K27ac model, or H3K4me3 model)          |

#### Optional Arguments
|Argument               | Description                                                              |
|-----------------------|--------------------------------------------------------------------------|
|--outFileName          | File name to use for output                                              |
|--outDirectory         | Directory to write output to                                             |

#### Commandline Usage
```
python 4_getTransferModel.py [FeatureFilePath] [ModelFilePath]
```

## 5_RunModel

Generate HM predictions from DNA methylation feature files.

#### Mandatory Arguments
|Argument            | Description                                                                                                    |
|--------------------|----------------------------------------------------------------------------------------------------------------|
|FeatureFilePath     | Full path to methylation features .h5 file (output from step 3_CombineInput)                                   |
|ModelFilePath       | Full path to model .h5 file (Provided H3K27ac model, or H3K4me3 model), alternatively a transfer-learned model |

#### Optional Arguments
|Argument               | Description                                                              |
|-----------------------|--------------------------------------------------------------------------|
|--outFileName          | File name to use for output                                              |
|--outDirectory         | Directory to write output to                                             |

#### Commandline Usage
```
python 5_getPredictions.py [FeatureFilePath] [ModelFilePath]
```

## File definitions of required input 
#### PromoterDefinitions
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

#### methylation bed-like file
A bed like file of genomic positions with corresponding M-values, tab delimited:

|Column      | Description                                                            |
|------------|------------------------------------------------------------------------|
|chrom       | chromosome ID, e.g. 1,2,3 ...22                                        |
|pos         | position of 5' cytosine of a CpG on the positive strand                |      
|mval        | caluclated mvalue of a given CpG, typically M-value=log2(Beta/1-Beta)  |

## Authors

* **Justin Williams, Beisi Xu, Daniel Putnam, and Xiang Chen** 
