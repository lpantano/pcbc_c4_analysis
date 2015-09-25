clean_table = function(what, mat, meta_clean){
  i = meta_clean %>% filter(Diffname_short %in% what)
  m = mat[,i$UID]

  d = data.frame(row.names = colnames(m), cov=colnames(m))
  dds = DESeqDataSetFromMatrix(m,colData = d,design = ~ 1)
  dds = estimateSizeFactors(dds)
  mat_norm = counts(dds, normalized=T)
  cutoff = round(ncol(mat_norm)*0.9)
  clean = as.data.frame((mat_norm[rowSums(mat_norm>5)>cutoff,]))
  list(mat=clean, meta=i)
}

clean_and_sync_mirna_tables = function(metadata, mat, reads, collapse=FALSE){
  metadata_filtered <-
    metadata %>%
    filter(Diffname_short != "") %>%
    filter(UID %in% colnames(mat)) %>%
    filter(Cell_Type == "PSC") %>%
    filter(C4_Karyotype_Result != "abnormal")

  REMOVED_UID <- setdiff(colnames(mat), metadata_filtered$UID)
  metadata <- metadata[metadata_filtered$UID,]

  for (nc in setdiff(colnames(metadata), "UID") )
    metadata[,nc] = as.factor(substr(gsub("[^[:alnum:]]", "", metadata[,nc]), 1, 40))

  mat <- mat[, metadata$UID]

  reads = reads[rownames(metadata),]
  keep = ( (reads %>% mutate(samples=rownames(reads)) %>% filter(mirna > 500000 & norm > 0.2))[,"samples"] )
  metadata$size = reads[as.character(row.names(metadata)),"norm"]
  metadata$size_cat = reads[as.character(row.names(metadata)),"breaks"]
  metadata$mirna = log2(reads[as.character(row.names(metadata)),"mirna"])
  # join meso-15 and meso-30
  stages =  gsub("-", "",metadata$Diffname_short)
  metadata$Diffname_short = as.character(stages)

  meta_clean = metadata[colnames(mat), ]
  meta_clean[meta_clean == 'N/A'] = NA
  meta_clean[is.na(meta_clean)] = 'unk'
  meta_clean[meta_clean == ''] = 'unk'

  meta_clean = meta_clean[keep,]
  mat = mat[,keep]

  keep_meta = meta_clean$Diffname_short %in% c("DE", "SC", "ECTO", "MESO5", "MESO15", "MESO30", "EB")
  meta_clean = meta_clean[keep_meta,]
  mat = mat[,keep_meta]

  PROCESSED_COUNTS <- lapply(
                             unique(meta_clean$Diffname_short),
                             function(state){
                               cols = row.names(meta_clean[meta_clean$Diffname_short==state,])
                               PROCESSED_COUNTS = getGeneFilteredGeneExprMatrix(mat[,cols],
                                                                                MIN_GENE_CPM=10,
                                                                                MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.9,
                                                                                verbose=TRUE)
                             })
  names(PROCESSED_COUNTS) = unique(meta_clean$Diffname_short)
  sapply(PROCESSED_COUNTS, function(x) print(x$plotHist))
  if (collapse){
    keep = unique(unlist(sapply(PROCESSED_COUNTS, function(x) x$filteredExprMatrix$genes[,1])))
    PROCESSED_COUNTS = list( filteredExprMatrix=DGEList(mat[keep, ], genes=row.names(mat[keep,])) )
  }

  return(list(count=PROCESSED_COUNTS, metadata=meta_clean))
}
