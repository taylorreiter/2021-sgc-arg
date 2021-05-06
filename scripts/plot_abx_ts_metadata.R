library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(purrr)

hmp_metadata <- read_csv("https://ibdmdb.org/tunnel/products/HMP2/Metadata/hmp2_metadata.csv")
h4017 <- hmp_metadata %>%
  select(participant_id = "Participant ID", data_id = "External ID", data_type, week_num, diagnosis, antibiotics = "Antibiotics") %>%
  filter(data_type == "metagenomics") %>%
  filter(data_id %in% c('HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
                        'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
                        'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD'))

ggplot(h4017, aes(x = week_num, y = antibiotics)) +
  geom_point() +
  theme_minimal()

# add in ARG info from groot

arg <- list.files("outputs/groot", full.names = T) %>%
  purrr::set_names() %>% 
  map_dfr(read_tsv, col_names = c("arg", "read_count", "gene_length", "coverage"), .id = "source") %>%
  mutate(source = gsub("outputs/groot/", "", source)) %>%
  filter(source != "all_arg90_report.txt") %>%
  mutate(source =gsub("_arg90_report.txt", "", source)) 

arg <- left_join(arg, h4017, by = c("source" = "data_id"))
ggplot(arg %>%
         group_by(week_num, antibiotics) %>%
         tally(), aes(x = week_num, y = n, color = antibiotics)) +
  geom_count() +
  theme_minimal()

# proportion of ARGs

proportion <- list.files("outputs/groot", pattern = "csv$", full.names = T) %>%
  map_dfr(read_csv, col_names = c("source", "proportion")) %>%
  left_join(h4017, by = c("source" = "data_id"))

pdf("../outputs/figures/prop_reads_arg.pdf", height = 2.5, width = 4.25)
ggplot(proportion, aes(x = week_num, y = proportion, fill = antibiotics)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "week number", y = "proportion of reads")
dev.off()


# inspect gather results ----------------------------------------------------------


gather <- list.files("outputs/gather", pattern = ".csv$", full.names = T) %>%
  purrr::set_names() %>% 
  map_dfr(read_csv, .id = "source") %>%
  mutate(source = gsub("outputs/gather/", "", source)) %>%
  mutate(source =gsub("_cfxA4_AY769933.fna.contigs.csv", "", source)) %>%
  separate(col = "source", into = c("source", "radius"), sep = "_") %>%
  left_join(h4017, by = c("source" = "data_id"))


gather_tmp <- gather %>%
  filter(radius == "r10")
table(gather_tmp$name, gather_tmp$source)
