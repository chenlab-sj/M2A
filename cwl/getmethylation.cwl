class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: get_methylation
baseCommand:
  - 2_getMethylationV2.py
inputs:
  - id: curated_sites
    type: File
    inputBinding:
      position: 0
  - id: promoterDefinitions
    type: File
    inputBinding:
      position: 2
outputs:
  - id: output
    type: File?
    outputBinding:
      glob: output/*Features*
label: getMethylation
requirements:
  - class: ResourceRequirement
    ramMin: 8000
    coresMin: 0
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
