library(tidyverse)
library(data.table)


bed <- fread("../debug/epe10.5_Cas_NGG_K562.UMI.ODN.trimmed.filtered.multi.realigned.bed",header=F)
bed <- bed %>%  mutate(V4 = str_remove(V4,"/2$"))


info <- fread("../debug/epe10.5_Cas_NGG_K562.UMI.ODN.trimmed.filtered.multi.realigned.info",header=F)
info <- info %>% mutate(AS = as.numeric(str_remove(V5,"AS:i:")),
                        FragLength = abs(V4)) %>% 
  select(AS,FragLength)


# bind bed and info -----------------------------------------------------------------

bed <- bed %>% bind_cols(info)


# Filter the bed file ---------------------------------------------------------------
# keep alignment with higher AS

bed <- bed %>% 
  group_by(V4) %>%
  slice_max(n = 1,order_by = AS,with_ties = T)


bed %>%ungroup %>%  count(V1,V2) %>% arrange(desc(n))
