#' Obtain site names for RES results
#'
#' Assigns each RES from a Sprint output table a site name denoting position and type of change
#'
#' @param RES.table data frame of RES table from Sprint obtained from read.delim() function
#'
#' @export
#'
RES.name <- function(RES.table){
  site.name = c()

  strand_name = rep("dot", nrow(RES.table))
  strand_name[which(as.character(RES.table$Strand) == "+")] = "plus"
  strand_name[which(as.character(RES.table$Strand) == "-")] = "minus"

  RES.name = c()
 for(i in 1:nrow(RES.table)){
    select.table = RES.table[i,]

    RES.name[i] = paste(select.table$X.Chrom, select.table$Start.0base., select.table$End.1base.,
                         select.table$Type, strand_name[i])
 }


  return(RES.name)
}

#' Obtain editing levels for results
#' Reformats results from Sprint file
#'
#' @param RES.table data frame of RES table from Sprint obtained from read.delim() function
#'
#' @export
#'
AD_DP <- function(RES.table){
  ADDP = as.character(RES.table$AD.DP)
  write.table(x = as.data.frame(ADDP), row.names = F, col.names = F, quote = F, file = "tempADDP.txt")
  ADDP = read.table("tempADDP.txt", sep = ":")
  system("rm tempADDP.txt")
  AD = ADDP$V1
  DP = ADDP$V2
  Editing_Level = AD/DP
  AD_DP = data.frame(RES.table, AD, DP, Editing_Level)
}

#' A function to intersect user region data with annotation data
#' Taken from annotatr package
#' Annotate genomic regions to selected genomic annotations while preserving the data associated with the genomic regions.
#'
#' @param regions The GRanges object to annotate
#' @param annotations The annotations to overlap with
#' @param minoverlap A scalar, positive integer, indicating the minimum required overlap of regions with annotations.
#' @param ignore.strand Logical indicating whether strandedness should be respected in findOverlaps(). Default FALSE.
#' @param quiet Print progress messages (FALSE) or not (TRUE).
#'
#'' @export
#'
annotate_regions <- function (regions, annotations = hg38_annotations, minoverlap = 1L, ignore.strand = TRUE,
          quiet = FALSE)
{
  if (class(regions)[1] != "GRanges") {
    stop("Error in annotate_regions(...): regions object is not GRanges.")
  }
  if (class(annotations)[1] != "GRanges") {
    stop("Error in annotate_regions(...): annotations object is not GRanges. Use build_annotations(...) to construct the annotations before calling annotate_regions(...).")
  }
  if (!quiet) {
    message("Annotating...")
  }
  intersections = GenomicRanges::findOverlaps(regions, annotations,
                                              minoverlap = minoverlap, ignore.strand = ignore.strand)
  if (length(intersections) > 0) {
    gr = regions[S4Vectors::queryHits(intersections)]
    GenomicRanges::mcols(gr)$annot = annotations[S4Vectors::subjectHits(intersections)]
    return(gr)
  }
  else {
    stop("No annotations intersect the regions.")
  }
}


#' A function for annotating the Sprint results with genomic regions and repeats
#'
#' @param input.file.path The input file path
#' @param annotation.granges.path The path to output initial annotations
#' @param complete.annotated.path The path to output the complete annotations
#'
#' @export
#'
run_ADARR <- function(input.file.path, annotation.granges.path, complete.annotated.path){

  input.sprint.read = read.delim(input.file.path)

  input.sprint = input.sprint.read[grep(x=as.character(input.sprint.read$X.Chrom),pattern = "chr"),]
  print(paste(nrow(input.sprint)," sites present"))

  res_name = RES.name(input.sprint)
  input.sprint = data.frame(input.sprint, res_name)
  input.sprint = AD_DP(input.sprint)
  input_granges = makeGRangesFromDataFrame(input.sprint,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("X.Chrom"),
                                           start.field="Start.0base.",
                                           end.field=c("End.1base."),
                                           strand.field="Strand",
                                           starts.in.df.are.0based=FALSE)


  input_annotated = annotate_regions(
    regions = input_granges,
    annotations = hg38_annotations,
    ignore.strand = FALSE,
    quiet = FALSE)


  input_annotated_df = as.data.frame(input_annotated, row.names = NULL)

  save(list="input_annotated", file = annotation.granges.path)

  ### reformat the metadata for hg38 repeats
  hg38_repeat_list = as.character(hg38_repeat_meta$repeat.)
  names(hg38_repeat_list) = hg38_repeat_meta$id

  print("Assigning annotations to RES")

  for(i in 1:nrow(input.sprint)){
    start.time = proc.time()
    select.res = as.character(input.sprint$res_name[i])
    subset.i = subset(input_annotated_df, as.character(res_name) == select.res)

    annot_identified = unique(as.character(subset.i$annot.id))
    genes_identified = unique(as.character(subset.i$annot.symbol))
    annot_identified = annot_identified[which(!is.na(annot_identified))]
    genes_identified = genes_identified[which(!is.na(genes_identified))]

    identified_annotations = paste(annot_identified,collapse = "|")
    identified_genes = paste(genes_identified, collapse = "|")

    if(length(annot_identified) == 0){
      identified_annotations = NA
    }

    if(length(genes_identified) == 0){
      identified_genes = NA
    }

    vector.select = annot_identified
    which.repeat = grep(x=vector.select, pattern = "repeat")

    if(length(which.repeat) > 0){
      repeat_assignments = paste(hg38_repeat_list[vector.select[which.repeat]], collapse = ";")
      #print(repeat_assignments)
    }else{
      repeat_assignments = NA
    }

    assign_annot = data.frame(identified_genes, identified_annotations, repeat_assignments)

    complete_annotated = cbind(input.sprint[i,],assign_annot)

    if(i == 1){

      write.table(x = complete_annotated, quote = F, sep = "\t", row.names = F,
                  file = complete.annotated.path)
    }else{
      write.table(x = complete_annotated, quote = F, sep = "\t", row.names = F,
                  file = complete.annotated.path, col.names = F, append = T)
    }

    total.time = proc.time() - start.time
    if(i %% 100 == 0){
      print("RES processed: ")
      print(i)
      print("Minutes left:")
      print(total.time[3]*(nrow(input.sprint)-i)/60)

    }

  }



  run_ADARR = NULL

  return(run_ADARR)

}



#' A function for annotating the Sprint results with genomic regions and repeats (allows for partition)
#'
#' @param input.file.path The input file path
#' @param annotation.granges.path The path to output initial annotations
#' @param complete.annotated.path The path to output the complete annotations
#' @param partition the number of RES to partition by
#'
#' @export
#'
run_ADARR_partition <- function(input.file.path, annotation.granges.path, complete.annotated.path, partition = 10000){

  ### reading initial files
  input.sprint.read = read.delim(input.file.path)

  input.sprint = input.sprint.read[grep(x=as.character(input.sprint.read$X.Chrom),pattern = "chr"),]
  print(paste(nrow(input.sprint)," sites present"))

  res_name = RES.name(input.sprint)
  input.sprint = data.frame(input.sprint, res_name)
  input.sprint = AD_DP(input.sprint)
  input_granges = makeGRangesFromDataFrame(input.sprint,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("X.Chrom"),
                                           start.field="Start.0base.",
                                           end.field=c("End.1base."),
                                           strand.field="Strand",
                                           starts.in.df.are.0based=FALSE)


  input_annotated = annotate_regions(
    regions = input_granges,
    annotations = hg38_annotations,
    ignore.strand = FALSE,
    quiet = FALSE)


  #input_annotated_df = as.data.frame(input_annotated, row.names = NULL)

  save(list="input_annotated", file = annotation.granges.path)

  rm(input_annotated)
  rm(input_granges)

  all.input.sprint = input.sprint

  rm(input.sprint)

  all.partitions = ceiling(nrow(all.input.sprint)/partition)

  initial.partitions = round(nrow(all.input.sprint)/partition)

  dir.create("partitions")
  all.indexes = list()
  for(m in 1:initial.partitions){
    all.indexes[[m]] = 1:partition + partition*(m-1)
    input.sprint = all.input.sprint[all.indexes[[m]],]
    save(list="input.sprint", file = paste("partitions/partition", m, ".RData", sep = ""))
    rm(input.sprint)

  }

  if(all.partitions > initial.partitions){
    m = all.partitions
    all.indexes[[all.partitions]] = max(unlist(all.indexes)):nrow(all.input.sprint)
    input.sprint = all.input.sprint[all.indexes[[all.partitions]],]
    save(list="input.sprint", file = paste("partitions/partition", m, ".RData", sep = ""))
    rm(input.sprint)
  }
  total_counts = nrow(all.input.sprint)
  rm(all.input.sprint)


  total_done = 0
  for(j in 1:all.partitions){

  load(paste("partitions/partition", j, ".RData", sep = ""))
  input_granges = makeGRangesFromDataFrame(input.sprint,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("X.Chrom"),
                                           start.field="Start.0base.",
                                           end.field=c("End.1base."),
                                           strand.field="Strand",
                                           starts.in.df.are.0based=FALSE)

  input_annotated = annotate_regions(
      regions = input_granges,
      annotations = hg38_annotations,
      ignore.strand = FALSE,
      quiet = FALSE)


  input_annotated_df = as.data.frame(input_annotated, row.names = NULL)
  ### reformat the metadata for hg38 repeats
  hg38_repeat_list = as.character(hg38_repeat_meta$repeat.)
  names(hg38_repeat_list) = hg38_repeat_meta$id


  print("Assigning annotations to RES")

  for(i in 1:nrow(input.sprint)){
    start.time = proc.time()
    select.res = as.character(input.sprint$res_name[i])
    subset.i = subset(input_annotated_df, as.character(res_name) == select.res)

    annot_identified = unique(as.character(subset.i$annot.id))
    genes_identified = unique(as.character(subset.i$annot.symbol))
    annot_identified = annot_identified[which(!is.na(annot_identified))]
    genes_identified = genes_identified[which(!is.na(genes_identified))]

    identified_annotations = paste(annot_identified,collapse = "|")
    identified_genes = paste(genes_identified, collapse = "|")

    if(length(annot_identified) == 0){
      identified_annotations = NA
    }

    if(length(genes_identified) == 0){
      identified_genes = NA
    }

    vector.select = annot_identified
    which.repeat = grep(x=vector.select, pattern = "repeat")

    if(length(which.repeat) > 0){
      repeat_assignments = paste(hg38_repeat_list[vector.select[which.repeat]], collapse = ";")
      #print(repeat_assignments)
    }else{
      repeat_assignments = NA
    }

    assign_annot = data.frame(identified_genes, identified_annotations, repeat_assignments)

    complete_annotated = cbind(input.sprint[i,],assign_annot)

    if(i == 1 & j == 1){

      write.table(x = complete_annotated, quote = F, sep = "\t", row.names = F,
                  file = complete.annotated.path)
    }else{
      write.table(x = complete_annotated, quote = F, sep = "\t", row.names = F,
                  file = complete.annotated.path, col.names = F, append = T)
    }

    total_done = total_done + 1


    total.time = proc.time() - start.time

    if(i == 1 & j == 1){
      time_span = total.time[3]
    }

    time_span = (total.time[3]+time_span)/2
    names(time_span) = NULL

    if(i %% 100 == 0){
      print("RES processed: ")
      print(total_done)
      print("Minutes left:")



      print(time_span*(total_counts - total_done)/60)

    }

  }

  rm(input.sprint)

}

  run_ADARR_partition = NULL

  return(run_ADARR_partition)

}



#' A function for annotating the Sprint results with genomic regions and repeats
#' This requires inputting a data frame rather than a file path
#'
#' @param input.sprint The input file path
#' @param annotation.granges.path The path to output initial annotations
#' @param complete.annotated.path The path to output the complete annotations
#'
#' @export
#'
run_ADARR_from_df <- function(input.sprint, annotation.granges.path, complete.annotated.path){


  input_granges = makeGRangesFromDataFrame(input.sprint,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("X.Chrom"),
                                           start.field="Start.0base.",
                                           end.field=c("End.1base."),
                                           strand.field="Strand",
                                           starts.in.df.are.0based=FALSE)


  input_annotated = annotate_regions(
    regions = input_granges,
    annotations = hg38_annotations,
    ignore.strand = FALSE,
    quiet = FALSE)


  input_annotated_df = as.data.frame(input_annotated, row.names = NULL)

  save(list="input_annotated", file = annotation.granges.path)

  ### reformat the metadata for hg38 repeats
  hg38_repeat_list = as.character(hg38_repeat_meta$repeat.)
  names(hg38_repeat_list) = hg38_repeat_meta$id

  print("Assigning annotations to RES")

  for(i in 1:nrow(input.sprint)){
    start.time = proc.time()
    select.res = as.character(input.sprint$res_name[i])
    subset.i = subset(input_annotated_df, as.character(res_name) == select.res)

    annot_identified = unique(as.character(subset.i$annot.id))
    genes_identified = unique(as.character(subset.i$annot.symbol))
    annot_identified = annot_identified[which(!is.na(annot_identified))]
    genes_identified = genes_identified[which(!is.na(genes_identified))]

    identified_annotations = paste(annot_identified,collapse = "|")
    identified_genes = paste(genes_identified, collapse = "|")

    if(length(annot_identified) == 0){
      identified_annotations = NA
    }

    if(length(genes_identified) == 0){
      identified_genes = NA
    }

    vector.select = annot_identified
    which.repeat = grep(x=vector.select, pattern = "repeat")

    if(length(which.repeat) > 0){
      repeat_assignments = paste(hg38_repeat_list[vector.select[which.repeat]], collapse = ";")
      #print(repeat_assignments)
    }else{
      repeat_assignments = NA
    }

    assign_annot = data.frame(identified_genes, identified_annotations, repeat_assignments)

    complete_annotated = cbind(input.sprint[i,],assign_annot)

    if(i == 1){

      write.table(x = complete_annotated, quote = F, sep = "\t", row.names = F,
                  file = complete.annotated.path)
    }else{
      write.table(x = complete_annotated, quote = F, sep = "\t", row.names = F,
                  file = complete.annotated.path, col.names = F, append = T)
    }

    total.time = proc.time() - start.time
    if(i %% 100 == 0){
      print("RES processed: ")
      print(i)
      print("Minutes left:")
      print(total.time[3]*(nrow(input.sprint)-i)/60)

    }

  }



  run_ADARR_from_df = NULL

  return(run_ADARR_from_df)

}


#' A function for annotating the Sprint results with genomic regions and repeats (allows for partition)
#'
#' @param input.file.path The input file path
#' @param partition.path Directory to export partitions to
#' @param partition the number of RES to partition by
#'
#' @export
#'
partition_files <- function(input.file.path, partition.path = "partitions", partition = 10000){

  ### reading initial files
  input.sprint.read = read.delim(input.file.path)

  input.sprint = input.sprint.read[grep(x=as.character(input.sprint.read$X.Chrom),pattern = "chr"),]
  print(paste(nrow(input.sprint)," sites present"))

  res_name = RES.name(input.sprint)
  input.sprint = data.frame(input.sprint, res_name)
  input.sprint = AD_DP(input.sprint)

  all.input.sprint = input.sprint

  rm(input.sprint)

  all.partitions = ceiling(nrow(all.input.sprint)/partition)

  initial.partitions = floor(nrow(all.input.sprint)/partition)

  #dir.create("partitions")
  all.indexes = list()
  for(m in 1:initial.partitions){
    all.indexes[[m]] = 1:partition + partition*(m-1)
    input.sprint = all.input.sprint[all.indexes[[m]],]
    save(list="input.sprint", file = paste(partition.path,"/partition", m, ".RData", sep = ""))
    print("file saved")
    print(m)
    print(nrow(input.sprint))
    rm(input.sprint)

  }

  if(all.partitions > initial.partitions){
    m=all.partitions
    all.indexes[[m]] = (max(unlist(all.indexes))+1):nrow(all.input.sprint)
    input.sprint = all.input.sprint[all.indexes[[m]],]
    save(list="input.sprint", file = paste(partition.path,"/partition", m, ".RData", sep = ""))
    print("file saved")
    print(m)
    print(nrow(input.sprint))
    rm(input.sprint)
  }
  total_counts = nrow(all.input.sprint)
  rm(all.input.sprint)

  partition_files = NULL
  return(partition_files)


}




#' A function for filtering annotated Sprint results by editing ratio and/or number of supporting reads
#'
#' @param input.sprint The Sprint data frame to input
#' @param editing_ratio The minimum editing ratio to have
#' @param supporting_reads The minimum number of supporting reads needed
#'
#' @export
#'
filter_annotated_results <- function(input.sprint, editing_ratio = 0, supporting_reads = 1){
filter.sprint1 = subset(input.sprint, Supporting_reads >= supporting_reads)
rm(input.sprint)
filter.sprint2 = subset(filter.sprint1, Editing_Level >= editing_ratio)
rm(filter.sprint1)
filter_annotated_results = filter.sprint2
rm(filter.sprint2)
return(filter_annotated_results)
}


#' A function for filtering annotated Sprint results by type of transitions
#'
#' @param input.sprint The Sprint data frame to input
#' @param transitions_list The list of transitions you want
#'
#' @export
#'
filter_by_transitions <- function(input.sprint, transitions_list = c("AG", "TC")){
  filter_by_transitions = input.sprint[which(as.character(input.sprint$Type) %in% transitions_list),]

  return(filter_by_transitions)
}



#' A function for obtaining table of identified genes
#'
#' @param input.sprint The Sprint data frame to input
#'
#' @export
#'
table_identified_genes <- function(input.sprint){
  gene_list = as.character(input.sprint$identified_genes)

  gene_table = as.data.frame(table(gene_list))

  colnames(gene_table) = c("gene", "RES")

  table_identified_genes = gene_table

  return(table_identified_genes)
}



#' A function for obtaining a gene matrix from a list of tables of identified genes
#'
#' @param list_table_identified_genes A list of tables from table_identified_genes() with names of tables being sample names
#'
#' @export
#'
getGeneMatrix <- function(list_table_identified_genes){

  len.list = length(list_table_identified_genes)
  names.list = names(list_table_identified_genes)


  all_genes <- c()

  for(j in 1:len.list){
    selected.table = list_table_identified_genes[[j]]
    rownames(selected.table) = selected.table$gene

    all_genes = c(all_genes, rownames(selected.table))
  }

  all_genes = unique(all_genes)

  ### The next step is the creation of a gene matrix with counts for each time a site is on a gene
  gene_RES_matrix = matrix(0,nrow = length(all_genes), ncol = len.list)

  rownames(gene_RES_matrix) = all_genes
  colnames(gene_RES_matrix) = names.list

  for(j in 1:len.list){
    selected.table = all_gene_res_list[[j]]
    rownames(selected.table) = selected.table$gene
    gene_RES_matrix[as.character(selected.table$gene),j] = selected.table$RES

  }

  getGeneMatrix = gene_RES_matrix
  return(getGeneMatrix)
}


