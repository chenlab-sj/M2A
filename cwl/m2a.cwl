class: Workflow
cwlVersion: v1.0
id: m2a
label: m2a
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: curated
    type: File
    'sbg:x': -580
    'sbg:y': 244
  - id: model
    type: File
    'sbg:x': -154
    'sbg:y': -396
  - id: chipBigwig
    type: File
    'sbg:x': -597.9090576171875
    'sbg:y': 89.5
  - id: inputBigwig
    type: File
    'sbg:x': -656.958251953125
    'sbg:y': -25.5
  - id: promoterDefinitions
    type: File
    'sbg:x': -639.671875
    'sbg:y': -296
outputs:
  - id: predictionsTxt
    outputSource:
      - getpredictions/output
    type: File
    'sbg:x': 920
    'sbg:y': 24
  - id: featuresTxt
    outputSource:
      - get_methylation/output
    type: File?
    'sbg:x': 905.9365234375
    'sbg:y': 386.6665954589844
  - id: featuresModel
    outputSource:
      - combine/outputModel
    type: File
    'sbg:x': 858.1072998046875
    'sbg:y': 232.57901000976562
steps:
  - id: get_methylation
    in:
      - id: curated_sites
        source: curated
      - id: promoterDefinitions
        source: promoterDefinitions
    out:
      - id: output
    run: ./getmethylation.cwl
    label: getMethylation
    'sbg:x': 13.206361770629883
    'sbg:y': 207.0237579345703
  - id: getpredictions
    in:
      - id: inputFeatures
        source: combine/outputModel
      - id: model
        source: gettransfermodel/updatedModel
    out:
      - id: output
    run: ./getpredictions.cwl
    label: getPredictions
    'sbg:x': 702
    'sbg:y': -18
  - id: gettransfermodel
    in:
      - id: FeatureFilePath
        source: combine/outputModel
      - id: ModelFilePath
        source: model
    out:
      - id: updatedModel
    run: ./gettransfermodel.cwl
    label: getTransferModel
    'sbg:x': 501
    'sbg:y': -49
  - id: getresponsevariable
    in:
      - id: chipBigwig
        source: chipBigwig
      - id: inputBigwig
        source: inputBigwig
      - id: promoterDefinitions
        source: promoterDefinitions
    out:
      - id: responseVariable
    run: ./getresponsevariable.cwl
    label: getResponseVariable
    'sbg:x': 83.14013671875
    'sbg:y': -154.5
  - id: combine
    in:
      - id: curatedFeatures
        source: get_methylation/output
      - id: responseVariable
        source: getresponsevariable/responseVariable
    out:
      - id: outputModel
    run: ./combine.cwl
    label: combine
    'sbg:x': 323.0421142578125
    'sbg:y': 120.5
requirements: []
