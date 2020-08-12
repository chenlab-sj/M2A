class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: gettransfermodel
baseCommand:
  - 5_getTransferModel.py
inputs:
  - id: FeatureFilePath
    type: File
    inputBinding:
      position: 0
  - id: ModelFilePath
    type: File
    inputBinding:
      position: 0
outputs:
  - id: updatedModel
    type: File
    outputBinding:
      glob: output/NewModel*.h5
label: getTransferModel
requirements:
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
