#' ---
#' title: "ADC2001 GCP retrospective study - circRNA analysis"
#' author: "Marieke Vromman"
#' output: 
#'    html_document:
#'       toc: TRUE
#'       toc_float: TRUE
#'       theme: paper
#'       df_print: paged
#'       highlight: tango
#'       code_folding: "show"
#'       number_sections: TRUE
#' ---

#' # Setup
#' set working dir, load lib, resolve conflicts
#+ message = FALSE

library(tidyverse)

library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("desc", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("slice", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer('intersect', 'dplyr')
conflict_prefer("count", "dplyr")

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


#' # read data

all_circ = read_tsv("./CircRNA/20210806_all_circ_retro.txt", col_types = 'ccccccccccccccccccccccccccccccccc')

all_circ = all_circ %>%
  separate(Sample_ID, into = c('main_diagnosis', "rest_diagnosis", "rest_diagnosis_2"), sep = " ", remove = F) %>%
  mutate(count = as.numeric(count))

all_circ = all_circ %>% select(chr, start, end, count, strand, tool, preprocessing, SMAPLE_RNA, SMAPLE_donor, circ_id, circ_id_strand, Subject_ID, Sample_ID, main_diagnosis, NAFLD_Diagnosis_AMC, NASH_Diagnosis_AMC, NASH_Atrisk_AMC)


samples = read_tsv("./CircRNA/compiled_data.tsv") %>%
  separate(Sample_ID, into = c('main_diagnosis', "rest_diagnosis", "rest_diagnosis_2"), sep = " ", remove = F)


n_reads = read_tsv("./CircRNA/nr_reads.txt") %>% 
  mutate(fastq = "unprocessed") %>%
  bind_rows(read_tsv("./CircRNA/nr_reads_dedup.txt") %>% 
              mutate(fastq = "dedup")) %>%
  full_join(samples) 

n_reads
  
fc_circ = all_circ %>% filter(tool == "find_circ")

#' nr of reads

ggplot(n_reads %>% filter(fastq == "unprocessed"), 
       aes(Sample_ID, nr_reads, fill = main_diagnosis)) +
  geom_bar(stat = "identity") +
  geom_bar(data = n_reads %>% filter(fastq == "dedup"), fill = "grey", alpha = 0.8, stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


#' # nr of circRNAs

ggplot(all_circ %>% filter(tool == "find_circ")) + 
  geom_bar(aes(Sample_ID, fill = main_diagnosis)) +
  geom_bar(data = all_circ %>% filter(tool == "CIRCexplorer2"), 
           mapping = aes(Sample_ID), fill = "grey", alpha = 0.9) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


#' nr of circRNAs, corrected for nr of reads

fc_circ %>% group_by(SMAPLE_RNA) %>%
  tally() %>%
  left_join(n_reads %>% filter(fastq == "dedup")) %>%
  mutate(CPM = 1000000 * n / nr_reads) %>%
  ggplot() +
  geom_bar(aes(Sample_ID, CPM, fill = main_diagnosis), stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylab("nr of circRNAs per 1M reads")


#' # circRNA counts
#' find_circ
all_circ %>% filter(tool == "find_circ") %>%
  ggplot(aes(x = Sample_ID, y = count)) + 
  geom_point(aes(color = main_diagnosis)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  geom_boxplot(outlier.alpha = 0) + # add boxplots without outliers
  scale_y_continuous(trans='log10')
  
#' CE2
all_circ %>% filter(tool == "CIRCexplorer2") %>%
  ggplot(aes(x = Sample_ID, y = count)) + 
    geom_point(aes(color = main_diagnosis)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    geom_boxplot(outlier.alpha = 0) + # add boxplots without outliers
    scale_y_continuous(trans='log10')



#' # circRNAs that are present in multiple samples

fc_circ %>% 
  group_by(circ_id_strand) %>%
  tally() %>%
  ggplot(aes(n)) +
  geom_bar() +
  scale_y_continuous(trans='log10') +
  xlab("# samples in which a circRNA is found") +
  ylab('# circRNA')


fc_circ %>% 
  group_by(circ_id_strand) %>%
  tally() %>%
  filter(n > 111)



#' # diff exp

#' NAs are fitlered out for diff exp
samples %>% count(NASH_Atrisk_AMC)
samples %>% count(NASH_Atrisk_AMC)
samples %>% count(NASH_Atrisk_AMC)


library("DESeq2")
library(ggrepel)



#' # diff exp - spike normalisation
#' code based on code from Jasper A.
counttable = read_tsv("./CircRNA/htseq_counts_noCC.tsv")

counttable = as_tibble(counttable) %>% 
  mutate_all(funs(replace(., is.na(.), 0))) %>% #replace NA by 0 
  select(-ensembl_transcript_id)

#' get sum of sequin counts
seq_counts = filter(counttable, str_detect(ensembl_gene_id,"R1|R2")) %>% summarize_if(is.numeric, sum)

seq_counts


#'   geom_text_repel()

#' # diff exp - lin + circ - spike normalized

#' ## NASH_Diagnosis_AMC
#' make count matrix
count_table = fc_circ %>% filter (!(is.na(NASH_Diagnosis_AMC))) %>% arrange(SMAPLE_RNA) %>%
  select(circ_id, SMAPLE_RNA, count, NASH_Diagnosis_AMC) %>%
  #group_by(circ_id) %>% filter(any(count>4)) %>% ungroup() %>%
  select(circ_id, SMAPLE_RNA, count) %>%
  spread(SMAPLE_RNA, count, fill = 0)

count_table
counttable
#' add lin (but only for same samples)
count_table = count_table %>%
  bind_rows(counttable %>% rename(circ_id = ensembl_gene_id) %>% 
              select(count_table %>% colnames()))

row_names_ = count_table$circ_id
count_table = count_table %>%
  select(-circ_id)
count_matrix = as.matrix(count_table)
row.names(count_matrix) = row_names_

#count_matrix = count_matrix+1

col_data = samples %>%
  filter(!(is.na(NASH_Diagnosis_AMC))) %>% 
  arrange(SMAPLE_donor) %>%
  select(SMAPLE_donor, NASH_Diagnosis_AMC)
row_names_ = col_data$SMAPLE_donor
col_data = col_data %>%
  select(-SMAPLE_donor) %>%
  as.data.frame()
col_data$NASH_Diagnosis_AMC = factor(col_data$NASH_Diagnosis_AMC)
row.names(col_data) = row_names_


#' construct DESEQDataSet Object
dds = DESeqDataSetFromMatrix(countData = count_matrix, 
                             colData = col_data, 
                             design = ~ NASH_Diagnosis_AMC)

dds

#' run DESeq
seq_counts_select <- dplyr::select(seq_counts, colnames(count_matrix))
sizeFactors(dds) = estimateSizeFactorsForMatrix(seq_counts_select)
dds = DESeq(dds)
res = results(dds)

#' run DESeq
dds = estimateSizeFactors(dds, type = 'poscounts')
dds = DESeq(dds)
res = results(dds)

summary(res)

tmp = as.data.frame(res) %>%
  rownames_to_column("circ_id")

res_df = tmp[complete.cases(tmp), ]

res_df = res_df %>% left_join(read_tsv("./circRNA/circRNAs_host_gene.txt") %>% 
                                mutate(circ_id = paste(chr, ":", start, "-", end, sep = "")) %>%
                                select(circ_id, gene_id, gene_name))

res_df = res_df %>%
  mutate(diffexpressed = "NO") %>%
  mutate(diffexpressed = ifelse(log2FoldChange > 0.6 & padj < 0.05, "UP", diffexpressed)) %>%
  mutate(diffexpressed = ifelse(log2FoldChange < -0.6 & padj < 0.05, "DOWN", diffexpressed)) %>%
  mutate(label_on = NA) %>%
  mutate(label_on = ifelse(diffexpressed == "NO", NA, gene_name))

write_tsv(res_df %>% select(-label_on), "./circRNA/circ_lin_NASH_Diagnosis_AMC_cl_spikes.tsv")

ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label = label_on)) + 
  geom_point() +
  theme_minimal() +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_color_manual(values=c("black", "darkolivegreen3")) +
  geom_text_repel()



#' # visualize RNA amount 
#' by (sum non-spike counts) / (sum spike counts)

tmp = counttable %>% 
  filter(str_detect(ensembl_gene_id,"R1|R2")) %>%
  pivot_longer(cols = -ensembl_gene_id, names_to = "SMAPLE_donor", values_to = "count_nr") %>%
  group_by(SMAPLE_donor) %>%
  summarise(spike_count = sum(count_nr)) %>%
  full_join(counttable %>% 
              pivot_longer(cols = -ensembl_gene_id, names_to = "SMAPLE_donor", values_to = "count_nr") %>%
              filter(substr(ensembl_gene_id, 1, 4) == "ENSG" ) %>%
              group_by(SMAPLE_donor) %>%
              summarise(non_spike_count = sum(count_nr))) %>%
  mutate(RNA_amount = non_spike_count / spike_count) 

tmp = tmp %>% left_join(samples %>% select(SMAPLE_RNA, NAFLD_Diagnosis_AMC), by = c('SMAPLE_donor' = "SMAPLE_RNA")) %>%
  mutate(NAFLD_Diagnosis_AMC = as.character(NAFLD_Diagnosis_AMC))

tmp %>%
  ggplot(aes(SMAPLE_donor, RNA_amount, fill = NAFLD_Diagnosis_AMC)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tmp %>%
  ggplot(aes(SMAPLE_donor, non_spike_count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

tmp %>%
  ggplot(aes(SMAPLE_donor, spike_count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

t.test(tmp %>% filter(NAFLD_Diagnosis_AMC == "0") %>% pull(RNA_amount), 
       y = tmp %>% filter(NAFLD_Diagnosis_AMC == "1") %>% pull(RNA_amount)) 

wilcox.test(tmp %>% filter(NAFLD_Diagnosis_AMC == "0") %>% pull(RNA_amount), 
            tmp %>% filter(NAFLD_Diagnosis_AMC == "1") %>% pull(RNA_amount))

tmp %>% 
  ggplot(aes(NAFLD_Diagnosis_AMC, RNA_amount)) +
  geom_boxplot() +
  scale_y_log10()
