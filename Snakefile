rule download_reads:
  input:
    "raw_data_accessions.sh"
  run:
    shell("sh {input}")


