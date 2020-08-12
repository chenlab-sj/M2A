class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: getresponsevariable
baseCommand:
  - 1_getResponseVariable.py
inputs:
  - id: chipBigwig
    type: File
    inputBinding:
      position: 0
  - id: inputBigwig
    type: File
    inputBinding:
      position: 1
  - id: promoterDefinitions
    type: File
    inputBinding:
      position: 3
outputs:
  - id: responseVariable
    type: File
    outputBinding:
      glob: output/*.txt
label: getResponseVariable
requirements:
  - class: ResourceRequirement
    ramMin: 2000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
