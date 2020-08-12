class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: parsegencodegff
baseCommand:
  - 01_parseGencodeGFF.py
inputs:
  - id: gencodeGFF
    type: File
    inputBinding:
      position: 0
outputs:
  - id: parsedGencode
    type: File
    outputBinding:
      glob: output/1_Gencode_parsedGFF.txt
label: parseGencodeGFF
requirements:
  - class: ResourceRequirement
    ramMin: 2000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
