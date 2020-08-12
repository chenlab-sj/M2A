class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: concatenatemethylation
baseCommand:
  - combine.sh
inputs:
  - id: methylation
    type: 'File[]'
    inputBinding:
      position: 0
      shellQuote: false
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: Features_Curated.txt
label: concatenateMethylation
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
  - class: InlineJavascriptRequirement
