class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: combine
baseCommand:
  - 3_Combine.py
inputs:
  - id: curatedFeatures
    type: File
    inputBinding:
      position: 0
  - id: responseVariable
    type: File?
    inputBinding:
      position: 0
      prefix: '--ResponseVariablePath'
  - id: outFileName
    type: string?
    inputBinding:
      position: 0
      prefix: '--outFileName'
  - id: outDirectory
    type: string?
    inputBinding:
      position: 0
      prefix: '--outDirectory'
outputs:
  - id: outputModel
    type: File
    outputBinding:
      glob: output/*.h5
label: combine
requirements:
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
