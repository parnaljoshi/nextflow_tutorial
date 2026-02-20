nextflow.enable.dsl=2

params.query_fasta = "query.fasta"
params.db_query    = "reviewed:true%20AND%20organism_id:559292"   // reviewed yeast proteins
params.outdir      = "results"

process DOWNLOAD_UNIPROT_FASTA {

  output:
  path "uniprot_yeast.fasta"

  script:
  """
  curl -L \
    "https://rest.uniprot.org/uniprotkb/stream?query=${params.db_query}&format=fasta&compressed=false" \
    -o uniprot_yeast.fasta
  """
}

process MAKE_BLAST_DB {

  input:
  path db_fasta

  output:
  path "uniprot_yeast.*"

  script:
  """
  makeblastdb -in ${db_fasta} -dbtype prot -out uniprot_yeast
  """
}

process BLASTP {

  publishDir "${params.outdir}", mode: 'copy'

  input:
  path query_fasta
  path db_files

  output:
  path "blast.tsv"

  script:
  """
  blastp \
    -query ${query_fasta} \
    -db uniprot_yeast \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    > blast.tsv
  """
}

workflow {
  q = file(params.query_fasta)
  db_fa   = DOWNLOAD_UNIPROT_FASTA()
  db_idx  = MAKE_BLAST_DB(db_fa)
  BLASTP(q, db_idx)
}
