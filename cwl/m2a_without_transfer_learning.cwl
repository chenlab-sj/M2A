class: Workflow
cwlVersion: v1.0
id: m2a_without_transfer_learning
label: m2a_without_transfer_learning
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: curated_sites
    type: File
    'sbg:x': -772.3968505859375
    'sbg:y': -54.5
  - id: model
    type: File
    'sbg:x': -61.3968505859375
    'sbg:y': -331.5
  - id: promoterDefinitions
    type: File
    'sbg:x': -751
    'sbg:y': -255
outputs:
  - id: predictionsTxt
    outputSource:
      - getpredictions/output
    type: File
    'sbg:x': 546.6031494140625
    'sbg:y': -159.5
  - id: featuresTxt
    outputSource:
      - get_methylation/output
    type: File?
    'sbg:x': 540
    'sbg:y': 211
  - id: featuresModel
    outputSource:
      - combine/outputModel
    type: File
    'sbg:x': 495.6031494140625
    'sbg:y': 54.5
steps:
  - id: get_methylation
    in:
      - id: curated_sites
        source: curated_sites
      - id: promoterDefinitions
        source: promoterDefinitions
    out:
      - id: output
    run: ./getmethylation.cwl
    label: getMethylation
    'sbg:x': -360
    'sbg:y': -150
  - id: combine
    in:
      - id: curatedFeatures
        source: get_methylation/output
    out:
      - id: outputModel
    run: ./combine.cwl
    label: combine
    'sbg:x': -83
    'sbg:y': -202
  - id: getpredictions
    in:
      - id: inputFeatures
        source: combine/outputModel
      - id: model
        source: model
    out:
      - id: output
    run: ./getpredictions.cwl
    label: getPredictions
    'sbg:x': 243
    'sbg:y': -186
requirements: []
