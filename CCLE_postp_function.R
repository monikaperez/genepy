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
source('gkugener/RScripts/load_libraries_and_annotations.R')


process_segments <- function(segment_file) {


  new_copy_number <- readr::read_tsv(segment_file, col_types = cols(
    DepMap_ID = col_character(), 
    CONTIG=col_character(), 
    START=col_double(), END=col_double(), 
    NUM_POINTS_COPY_RATIO=col_double(), MEAN_LOG2_COPY_RATIO=col_double(), 
    CALL=col_character()
  )) %>% dplyr::select(DepMap_ID, Chromosome=CONTIG, Start=START, End=END, Num_Probes=NUM_POINTS_COPY_RATIO, Segment_Mean=MEAN_LOG2_COPY_RATIO) %>%
    dplyr::mutate(Source=case_when(
      grepl('^(dm|ccle2)_', DepMap_ID) ~ 'Broad WES',
      grepl('^chordoma_', DepMap_ID) ~ 'Chordoma WES'),
      'Other WES')

  # Untransform data if it is log2 transformed (which it will be if from the FireCloud pipeline)
  if (min(new_copy_number$Segment_Mean) < 0) {
    new_copy_number %<>% mutate(Segment_Mean=2^Segment_Mean)
  }
  return(new_copy_number)
}

filter_for_CCLE <- function(new_copy_number){

  # We shouldn't have anything labelled other, so check this
  if (nrow(new_copy_number %>% filter(Source=='Other WES')) > 0) {
    print('ERROR. THERE IS A SAMPLE NOT FROM BROAD OR THE CHORDOMA FOUNDATION')
  }

  # Make sure there are no duplicates
  counts_unique_by_source <- new_copy_number %>%
    mutate(DepMap_ID=stringr::str_extract(string=Sample, pattern='ACH\\-[0-9]+')) %>%
    distinct(Source, DepMap_ID, Sample) %>%
    group_by(Source, DepMap_ID) %>%
    dplyr::summarise(count=n()) %>%
    filter(count > 1)

  if (nrow(counts_unique_by_source) != 0) {
    print('ERROR. THERE IS A DUPLICATE SAMPLE IN THE SET FOR A SINGLE SOURCE')
  }
  # If it passes above, then can remove prefix from DepMap IDs
  new_copy_number %<>% dplyr::mutate(DepMap_ID=stringr::str_extract(pattern='ACH\\-[0-9]+', string = DepMap_ID))
  return(new_copy_number)
}


interpolate_gaps_in_segmented <- function(segments) {
  segments_as_granges <- GenomicRanges::makeGRangesFromDataFrame(segments, keep.extra.columns = T)
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
extend_ends_of_segments <- function(segments, hg38_cyto_band_reference='data/hg38_cytoband.gz') {
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




# This chunk is kept in here for reference, but does not need to be run for the release
# This was used to combine the copy number calls from Sanger and the Broad from WES
# For the Broad samples, we want the calls that use the ICE/AGILENT PON for chr1-22,X
# For Sanger samples we want the calls that use the SANGER specific AGILENT PON for chr1-22 and then the AGILENT PONT for X
# This is what this chunk accomplished and then saved to a tsv for upload to taiga

# Format for upload to taiga
combine_broad_sanger <- function(Sanger_filename= 'Downloads/Sanger.called.seg',
  Broad_filename = '~/Downloads/DM19Q2_COMPLETE.called.seg', outfile="'~/Downloads/all_wes_data_19q2_v2.tsv'"){
  agilent_ice_pon_based_calls <- readr::read_tsv(filename, col_types = cols(
    Sample = col_character(), 
    CONTIG=col_character(), 
    START=col_double(), END=col_double(), 
    NUM_POINTS_COPY_RATIO=col_double(), MEAN_LOG2_COPY_RATIO=col_double(), 
    CALL=col_character()
  )) %>% dplyr::select(Sample, Chromosome=CONTIG, Start=START, End=END, Num_Probes=NUM_POINTS_COPY_RATIO, Segment_Mean=MEAN_LOG2_COPY_RATIO) %>%
    dplyr::mutate(Source=ifelse(grepl('sanger', Sample), 'Sanger WES', ifelse(grepl('chordoma', Sample), 'Chordoma WES', 'Broad WES'))) %>%
    # We only keep the X calls from here 
    filter((Source %in% c('Chordoma WES', 'Broad WES')) | (Source == 'Sanger WES' & gsub('chr', '', Chromosome) == 'X'))

  sanger_pon_based_calls <- readr::read_tsv(, col_types = cols(
    Sample = col_character(), 
    CONTIG=col_character(), 
    START=col_double(), END=col_double(), 
    NUM_POINTS_COPY_RATIO=col_double(), MEAN_LOG2_COPY_RATIO=col_double(), 
    CALL=col_character()
  )) %>% dplyr::select(Sample, Chromosome=CONTIG, Start=START, End=END, Num_Probes=NUM_POINTS_COPY_RATIO, Segment_Mean=MEAN_LOG2_COPY_RATIO) %>%
    dplyr::mutate(Source='Sanger WES') %>%
    # We remove the sex chromosomes from here as this was a mixed panel
    filter(gsub('chr', '', Chromosome) %in% seq(1,22))

  # Combine two datasets above, validate that only one cell line per source, change the DepMap_ID column to have only the DepMap_ID
  combined_copynumber_calls <- rbind(sanger_pon_based_calls, agilent_ice_pon_based_calls)

  counts_unique_by_source <- combined_copynumber_calls %>%
    mutate(DepMap_ID=stringr::str_extract(string=Sample, pattern='ACH\\-[0-9]+')) %>%
    distinct(Source, DepMap_ID, Sample) %>%
    group_by(Source, DepMap_ID) %>%
    dplyr::summarise(count=n()) %>%
    filter(count > 1)

  if (nrow(counts_unique_by_source) != 0) {
    print('ERROR. THERE IS A DUPLICATE SAMPLE IN THE SET FOR A SINGLE SOURCE')
  }

  combined_copynumber_calls %<>% mutate(DepMap_ID=stringr::str_extract(pattern='ACH\\-[0-9]+', string=Sample)) %>%
    dplyr::select(-Sample) %>%
    mutate(Chromosome=gsub('chr', '', Chromosome)) %>%
    arrange(DepMap_ID, Chromosome)

  write.table(combined_copynumber_calls, file = outfile, sep = '\t', quote = F, row.names = F) # What we upload to taiga
}


process_counts <- function(filepath){
  # Process counts
  counts_genes <- read_tsv(
    file = filepath
    # col_types = cols(.default = "c")
  )

  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  counts_samples_to_add_one_off <- c('~/Downloads/dm_ACH-000309.rsem.genes.results', '~/Downloads/ibm_ACH-001852.rsem.genes.results')

  if (length(counts_samples_to_add_one_off) > 0) {
    for (f in counts_samples_to_add_one_off) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      
      # Check that these samples are not already present in the dataset
      if (s_id %in% colnames(counts_genes)) {
        print(paste0(s_id, ' is already in the dataset, skipping...'))
        next
      }
      
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, `transcript_id(s)`, expected_count) %>%
        set_colnames(c('gene_id', 'transcript_id(s)', s_id))
      
      counts_genes %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id(s)'))
    }
  }
}

process_tpm <- function(filepath){
  # TPM (genes)
  tpm_genes <- read_tsv(
    file = filepath
    # col_types = cols(.default = "c")
  )

  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  samples_to_add_one_off <- c('~/Downloads/dm_ACH-000309.rsem.genes.results', '~/Downloads/ibm_ACH-001852.rsem.genes.results')

  if (length(samples_to_add_one_off) > 0) {
    for (f in samples_to_add_one_off) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      
      # Check that these samples are not already present in the dataset
      if (s_id %in% colnames(counts_genes)) {
        print(paste0(s_id, ' is already in the dataset, skipping...'))
        next
      }
      
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, `transcript_id(s)`, TPM) %>%
        set_colnames(c('gene_id', 'transcript_id(s)', s_id))
      
      tpm_genes %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id(s)'))
    }
  }
}

process_transcripts <- function(filepath){
  # Transcripts (genes)
  transcripts <- read_tsv(
    file = ,
    # col_types = cols(.default = "c")
  )

  # loop to add one offs (in the case that a subset of the samples didn't run in time as happened in 19Q2)
  transcripts_samples_to_add_one_off <- c('~/Downloads/dm_ACH-000309.rsem.isoforms.results', '~/Downloads/ibm_ACH-001852.rsem.isoforms.results')

  if (length(transcripts_samples_to_add_one_off) > 0) {
    for (f in transcripts_samples_to_add_one_off) {
      s_id <- stringr::str_extract(pattern = '(dm|ibm|ccle2)_ACH\\-[0-9]+', string = f)
      data_to_add <- read_tsv(f) %>%
        dplyr::select(gene_id, transcript_id, TPM) %>%
        set_colnames(c('gene_id', 'transcript_id', s_id))
      
      transcripts %<>% left_join(., data_to_add, by=c('gene_id', 'transcript_id'))
    }
  }
}

compare_releases <- function(){
  # Compare the previous release to the current release (using correlations)
  tentative_new_release_tpm <- tpm_genes %>% 
    mutate(gene=stringr::str_extract(string=gene_id, pattern='ENSG[0-9]+')) %>%
    filter(!grepl('PAR_Y', gene_id)) %>% # Not sure what these are
    filter(!grepl('ERCC', gene_id)) %>% # These are the ERCC spike ins
    dplyr::select(c('gene', colnames(.)[grepl('ACH\\-[0-9]+', colnames(.))])) %>%
    column_to_rownames(var='gene') %>%
    t() %>%
    # If we have two samples with the same arx span id (we shouldn't if the sample set is designed correctly) then it will complain below
    set_rownames(stringr::str_extract(string=row.names(.), pattern = 'ACH\\-[0-9]+'))

  # Now do correlations (of log2+1 TPM)
  overlap_genes <- intersect(colnames(previous_release_tpm), colnames(tentative_new_release_tpm))
  overlap_cell_lines <- intersect(row.names(previous_release_tpm), row.names(tentative_new_release_tpm))

  # Check to see if any cell lines from the previous release are not present in this dataset (this should not be the case unless there is a known processing error, so this list should be empty)
  row.names(previous_release_tpm) %>% setdiff(row.names(tentative_new_release_tpm))

  # Correlations of samples (could also just look at the set of most variable )
  # Intersect the top 2000 most variable in both
  tpm_19q1_most_variable <- apply(previous_release_tpm, 2, sd) %>% .[order(-.)] %>% names() %>% .[1:2000]
  tpm_19q2_most_variable <- apply(log2(tentative_new_release_tpm+1), 2, sd) %>% .[order(-.)] %>% names() %>% .[1:2000]

  # 95% overlap
  most_variable_for_correlations <- intersect(tpm_19q1_most_variable, tpm_19q2_most_variable)
  length(most_variable_for_correlations)/2000

  correlation_rnaseq_data_releases <- cor(
    t(previous_release_tpm[overlap_cell_lines, most_variable_for_correlations]),
    t(log2(tentative_new_release_tpm[overlap_cell_lines, most_variable_for_correlations]+1)),
  )
}

# TODO: process the exons in future releaes...

rename_function <- function(columns) {
  columns_new <- ifelse(columns %in% c('Name', 'Description', 'gene_id', 'transcript_id', "transcript_id(s)"), columns, ifelse(
    grepl('ACH\\-[0-9]+$', columns), 
    stringr::str_extract(string=columns, pattern='ACH\\-[0-9]+'), ccle.to.arxspan(columns, ignore.problems = T)
  ))
  
  return(columns_new)
}

