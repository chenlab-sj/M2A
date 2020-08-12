class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
id: getchromosomes
baseCommand: []
inputs:
  - id: curated
    type: File
outputs:
  - id: output
    type: 'string[]'
    outputBinding:
      loadContents: true
      glob: stdout.txt
      outputEval: '$(self[0].contents.split("\n").filter(Boolean))'
label: getChromosomes
arguments:
  - position: 0
    prefix: ''
    shellQuote: false
    valueFrom: tail -n +2 $(inputs.curated.path) | cut -f 1 | sort | uniq
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 1000
    coresMin: 1
  - class: DockerRequirement
    dockerPull: 'stjude/m2a:0.0.1'
  - class: InlineJavascriptRequirement
stdout: stdout.txt
