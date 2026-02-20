nextflow.enable.dsl=2

process hello {
  publishDir 'output', mode: 'copy'

  output:
  path 'result.txt'

  script:
  """
  echo "Welcome to Nextflow!" > result.txt
  """
}

workflow {
  // Run the hello process
  hello()
}
