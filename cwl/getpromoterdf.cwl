class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: getpromoterdf
baseCommand:
  - 02_getPromoterDF.py
inputs:
  - id: parsedGencode
    type: File
    inputBinding:
      position: 0
outputs:
  - id: output
    type: File
    outputBinding:
      glob: output/Promoter_Definitions.txt
label: getPromoterDF
requirements:
  - class: ResourceRequirement
    ramMin: 2000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
