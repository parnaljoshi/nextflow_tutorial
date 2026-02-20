nextflow.enable.dsl=2

params.query_fasta = "data/query.fasta"
params.db_query    = "reviewed:true%20AND%20organism_id:559292"   // reviewed yeast proteins
params.outdir      = "results"
params.chunkSize = 5

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
  tuple val(chunk_id), path(query_chunk)
  path db_files

  output:
  path "blast_${chunk_id}.tsv"

  script:
  """
  blastp \
    -query ${query_chunk} \
    -db uniprot_yeast \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    > blast_${chunk_id}.tsv
  """
}

process MERGE_BLAST {

  publishDir "${params.outdir}", mode: 'copy'

  input:
  path blast_parts

  output:
  path "blast.tsv"

  script:
  """
  cat *.tsv > blast.tsv
  """
}

workflow {

  q = file(params.query_fasta)
  db_fa   = DOWNLOAD_UNIPROT_FASTA()
  db_idx  = MAKE_BLAST_DB(db_fa)

  ch_chunks = Channel
    .fromPath(params.query_fasta)
    .splitFasta(by: params.chunkSize, file: true)
    .map { f -> tuple( "chunk_${f.baseName}", f ) }
  
  ch_blast_parts = BLASTP(ch_chunks, db_idx)

  MERGE_BLAST(ch_blast_parts.collect())
}
