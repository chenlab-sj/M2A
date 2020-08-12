class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: getpredictions
baseCommand:
  - 4_getPredictions.py
inputs:
  - id: inputFeatures
    type: File
    inputBinding:
      position: 0
  - id: model
    type: File
    inputBinding:
      position: 0
outputs:
  - id: output
    type: File
    outputBinding:
      glob: output/Predictions*.txt
label: getPredictions
requirements:
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
