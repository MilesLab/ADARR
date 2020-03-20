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
