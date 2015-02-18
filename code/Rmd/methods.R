
get_metadata = function(project){
    require(synapseClient)
    metadataTable <- synGet(project)

    colsToUse <- c("UID", "C4_Cell_Line_ID", "bad_lines", "pass_qc", "SampleType", "Cell_Type", "Cell_Line_Type",
                   "Cell_Type_of_Origin", "Cell_Line_of_Origin", "Tissue_of_Origin", "Reprogramming_Vector_Type",
                   "Reprogramming_Gene_Combination", "C4_Karyotype_Result", "High_Confidence_Donor_ID")
    q <- sprintf("SELECT * FROM %s", metadataTable@properties$id)

    metadata <- synTableQuery(q)@values

    ## Remove fail qc and bad lines
    metadataFinal <- filter(metadata, !bad_lines, pass_qc)

    rownames(metadataFinal) <- metadataFinal$UID
    metadataFinal
}
