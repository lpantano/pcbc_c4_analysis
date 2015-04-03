library(org.Hs.eg.db)
library(AnnotationDbi)
# miRDB_v5. Need file from (better filtered by hsa first): http://mirdb.org/miRDB/download.html
nm = read.table("~/repos/pcbc_c4_analysis/data/targets/hsa_miRDB_v5.0_prediction_result.txt")
symbol = AnnotationDbi::select(org.Hs.eg.db, as.character(unlist(nm[,2])), "ENSEMBL", keytype="REFSEQ")
write.table(symbol, "~/repos/pcbc_c4_analysis/data/targets/hsa_refseq_mirdb_ensembl",row.names=F,quote=F,col.names=F,sep="\t")

# mirtarbase. CSV of whole file from web: http://mirtarbase.mbc.nctu.edu.tw/php/download.php
mtb = read.csv("~/repos/pcbc_c4_analysis/data/targets/hsa_MTI.csv")
symbol = AnnotationDbi::select(org.Hs.eg.db, unique(as.character(unlist(mtb[,4]))), "ENSEMBL", keytype="SYMBOL")
write.table(symbol, "~/repos/pcbc_c4_analysis/data/targets/hsa_symbol_mirtarbase_ensembl",row.names=F,quote=F,col.names=F,sep="\t")


