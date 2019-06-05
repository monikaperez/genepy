## Guillaume Kugener
## for BroadInsitute
## in 2017

######################
# The function below fills in the gaps in the segmented level data
#
#	In this release, the following changes were made:
#
#	* Number of new cell lines: `r length(new_cell_lines)`
#	* Number of cell lines moving from Sanger WES/Broad SNP -> Broad WES CN: `r length(replaced_cell_lines)`
#	* Total cell lines with CN this release: `r length(unique(combined_new_prioritized_dataset$DepMap_ID))`
#
#	## Interpolate segments
#
#	In this section we perform to operations on the segment level data
#
#	1. Fill in the gaps: there may be gaps between segments in the copy number data, leading to the possibility genes mapping to these gaps and being NAed further downstream.
#	2. Extend the ends: there are genes that map outside the targeting regions (in WES). To address these cases, we can extend the ends of the segments so that these genes are not NAed.
# 
# @args:
#   - segments: 
#       data.frame with DepMap_ID (what samples are separated on), 
#       Chromosome, Start, End, Segment_Mean, Num_Probes
#
# @returns: a segments data.frame of the same size with gaps in ranges filled
######################
#
library(org.Hs.eg.db) # This is using hg38 
source('~/gkugener/RScripts/load_libraries_and_annotations.R')


interpolate_gaps_in_segmented <- function(segment_file) {
  segments_as_granges <- GenomicRanges::makeGRangesFromDataFrame(segment_file, keep.extra.columns = T)
  segments_as_granges_list <- split(segments_as_granges, segments_as_granges$DepMap_ID)
  
  # Determine which ones have too many gaps
  segments_gaps_filled <- NULL
  count <- 0
  disjoin_different_cls <- c()
  for (cl in names(segments_as_granges_list)) {
    if (count %% 100 == 0) {
      print(count/length(names(segments_as_granges_list))* 100)
    }
    count <- count + 1
    
    if (length(segments_as_granges_list[[cl]] %>% disjoin()) != length(segments_as_granges_list[[cl]])) {
      disjoin_different_cls %<>% c(., cl)
      source_using_current <- unique(segments_as_granges_list[[cl]]$Source)
      
      if (source_using_current != 'Broad SNP') {
        print(paste0('Disjoin different for: ', cl, ' source: ', source_using_current))
      }
      
      ttt <- segments_as_granges_list[[cl]] %>%
        as.data.frame() %>%
        dplyr::select(-Source, -width, -DepMap_ID)
      
      ddd <- segments_as_granges_list[[cl]] %>% 
        disjoin() %>%
        as.data.frame()
      
      get_missing_in_disjoined <- ddd %>%
        left_join(., ttt, by=c('seqnames', 'start', 'end', 'strand')) %>%
        dplyr::rename(SM_start_end=Segment_Mean, Num_Probes_start_end=Num_Probes) %>%
        left_join(., ttt %>% dplyr::select(-end), by=c('seqnames', 'start', 'strand')) %>%
        dplyr::rename(SM_start=Segment_Mean, Num_Probes_start=Num_Probes) %>%
        left_join(., ttt %>% dplyr::select(-start), by=c('seqnames','end', 'strand')) %>%
        dplyr::rename(SM_end=Segment_Mean, Num_Probes_end=Num_Probes) %>%
        filter(!(is.na(SM_start_end) & !is.na(SM_start) & !is.na(SM_end))) %>%
        mutate(SM_final=case_when(!is.na(SM_start_end) ~ SM_start_end,
          !is.na(SM_start) ~ SM_start,
          !is.na(SM_end) ~ SM_end,
          TRUE ~ 1)) %>%
        mutate(Num_Probes_final=case_when(!is.na(Num_Probes_start_end) ~ Num_Probes_start_end,
          !is.na(Num_Probes_start) ~ Num_Probes_start,
          !is.na(Num_Probes_end) ~ Num_Probes_end,
          TRUE ~ as.integer(1))) %>%
        mutate(DepMap_ID=cl, Source=source_using_current) %>%
        dplyr::select(seqnames, start, end, Num_Probes=Num_Probes_final, Segment_Mean=SM_final, DepMap_ID, Source) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
      
      segments_as_granges_list[[cl]] <- get_missing_in_disjoined
    }
    # We want to extend the ranges to fill in any gaps. We use GenomicRanges gaps function to identify
    # the gaps in the segmented data that are present per cell line. If there is a gap between
    # two segments, we extend the ends to meet halfway between the gap. This is what the code snippet below
    # is performing
    segments_gaps_filled %<>% rbind(segments_as_granges_list[[cl]] %>%
      GenomicRanges::as.data.frame() %>%
      dplyr::select(-strand, -width) %>%
      left_join(., 
        gaps(segments_as_granges_list[[cl]]) %>% 
          GenomicRanges::as.data.frame() %>%
          dplyr::select(-strand, -width) %>%
          dplyr::mutate(right_start=start, right_end=end) %>%
          dplyr::mutate(end=start-1) %>%
          dplyr::select(-start),
        by=c('seqnames', 'end')) %>%
      left_join(., 
          gaps(segments_as_granges_list[[cl]]) %>% 
          GenomicRanges::as.data.frame() %>%
          dplyr::select(-strand, -width) %>%
          dplyr::mutate(left_start=start, left_end=end) %>%
          dplyr::mutate(start=end+1) %>%
          dplyr::select(-end),
        by=c('seqnames', 'start')) %>%
      # Edit the range values if needed
      dplyr::mutate(
        start=case_when(
          !is.na(left_end) ~ (left_start + floor((left_end-left_start)/2)+1),
          TRUE ~ start
        ),
        end=case_when(
          !is.na(right_end) ~ (right_start + floor((right_end-right_start)/2)),
          TRUE ~ end
      )) %>% dplyr::select(DepMap_ID, seqnames, start, end, Num_Probes, Segment_Mean, Source))
  } 
  return(list(segs=segments_gaps_filled, disjoin_different_cls=disjoin_different_cls))
}

# We use this method to extend the ends of the segments
extend_ends_of_segments <- function(segments) {
  # Starts to 1
  segments %<>%
    group_by(DepMap_ID, seqnames) %>%
    dplyr::mutate(m=min(start)) %>%
    dplyr::mutate(start=ifelse(start==m, 1, start)) %>%
    ungroup() %>%
    dplyr::select(-m)
  
  # Ends to max of chromosome value
  chromosome_maxes <- read_tsv(hg38_cyto_band_reference, col_names = F) %>%
    group_by(X1) %>% dplyr::summarise(max=max(X3)) %$% setNames(max, X1)
  
  segments %<>%
    group_by(DepMap_ID, seqnames) %>%
    dplyr::mutate(m=max(end)) %>%
    dplyr::mutate(end=ifelse(end==m, chromosome_maxes[paste0('chr', seqnames)], end)) %>%
    ungroup() %>%
    dplyr::select(-m)
  return(segments)
}

filter.coordinates <- function(loc, coordinate.end, max.spread=2000000) {
  ## Only keep coordinates on propper chromosomes
  loc <- as.numeric(loc[names(loc) %in% c(as.character(1:22), "X", "Y")])
  if (length(loc)==0) {
    return(NA)
  } else if (coordinate.end=="start") {
    ## start coordinate
    if (diff(range(loc)) > max.spread) {
      # print(loc)
    }
    return(ifelse(diff(range(loc))>max.spread, NA, min(abs(loc))))
  } else {
    ## end coordinate
    return(ifelse(diff(range(loc))>max.spread, NA, max(abs(loc))))
  }
}

## Only keep proper chromosomes
filter.chromosomes <- function(chrom) {
  chrom <- intersect(chrom, c(as.character(1:22), "X", "Y"))
  ifelse(length(chrom)>0, chrom, NA)
}