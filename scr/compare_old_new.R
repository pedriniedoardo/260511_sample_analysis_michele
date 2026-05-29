# aim ---------------------------------------------------------------------
# compare the new and the old dataset

data_new <- read.csv("data/PCAdata2.0.csv",row.names = "ms_id") %>%
  # filter only the MS patient
  filter(pt == 1) %>%
  dplyr::select(age, NfL, GFAP, CHIT,OPN) %>%
  # remove the two NA from OPN
  filter(!is.na(OPN)) %>%
  # now we do not have zeros anymore, therefore use the log instead of the log1p
  # mutate(NfL = log1p(NfL),
  #        GFAP = log1p(GFAP),
  #        CHIT = log1p(CHIT),
  #        OPN = log1p(OPN)) %>%
  rownames_to_column("sample_id") %>%
  select(-age) %>%
  pivot_longer(values_to = "con",names_to = "marker",-sample_id)

data_old <- read.csv("data/dataset_short_all.csv", row.names = "MS_ID") %>%
  # filter only the MS patient
  filter(CTR0_MS1 == 1) %>%
  dplyr::select(age = Age_at_prelievo, NfL, GFAP, CHIT,OPN) %>%
  # remove the two zeros from OPN
  filter(OPN>0) %>%
  # mutate(NfL = log1p(NfL),
  #        GFAP = log1p(GFAP),
  #        CHIT = log1p(CHIT),
  #        OPN = log1p(OPN)) %>%
  rownames_to_column("sample_id") %>%
  select(-age) %>%
  pivot_longer(values_to = "con",names_to = "marker",-sample_id)

# compare the values
test <- left_join(data_old,data_new,by = c("sample_id","marker"),suffix = c(".old",".new")) %>%
  mutate(delta = con.old - con.new)

test %>%
  filter(delta != 0) %>%
  arrange(marker,desc(delta)) %>%
  print(n = 40)

test %>%  
  ggplot(aes(x=con.old,y=con.new)) +
  geom_point(shape = 1) +
  geom_abline(slope = 1,intercept = 0,col="red") +
  facet_wrap(~marker,scales="free") +
  theme_bw()+theme(strip.background = element_blank())


# michele has updated the values ------------------------------------------

data_new <- read.csv("data/biomarkers_lod2.csv",row.names = "MS_ID") %>%
  dplyr::select(NfL, GFAP, CHIT,OPN) %>%
  # remove the two NA from OPN
  filter(!is.na(OPN)) %>%
  # now we do not have zeros anymore, therefore use the log instead of the log1p
  # mutate(NfL = log1p(NfL),
  #        GFAP = log1p(GFAP),
  #        CHIT = log1p(CHIT),
  #        OPN = log1p(OPN)) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(values_to = "con",names_to = "marker",-sample_id)

data_old <- read.csv("data/dataset_short_all.csv", row.names = "MS_ID") %>%
  # filter only the MS patient
  filter(CTR0_MS1 == 1) %>%
  dplyr::select(age = Age_at_prelievo, NfL, GFAP, CHIT,OPN) %>%
  # remove the two zeros from OPN
  filter(OPN>0) %>%
  # mutate(NfL = log1p(NfL),
  #        GFAP = log1p(GFAP),
  #        CHIT = log1p(CHIT),
  #        OPN = log1p(OPN)) %>%
  rownames_to_column("sample_id") %>%
  select(-age) %>%
  pivot_longer(values_to = "con",names_to = "marker",-sample_id)

# compare the values
test <- left_join(data_old,data_new,by = c("sample_id","marker"),suffix = c(".old",".new")) %>%
  mutate(delta = con.old - con.new)

test %>%
  filter(delta != 0) %>%
  arrange(marker,desc(delta)) %>%
  print(n = 40)

test %>%  
  ggplot(aes(x=con.old,y=con.new)) +
  geom_point(shape = 1) +
  geom_abline(slope = 1,intercept = 0,col="red") +
  facet_wrap(~marker,scales="free") +
  theme_bw()+theme(strip.background = element_blank())

# create the new dataset for the updated analysis

# read the new data anly to keep the metadata
data_new_meta <- read.csv("data/PCAdata2.0.csv",row.names = "ms_id") %>%
  # filter only the MS patient
  filter(pt == 1) %>%
  dplyr::select(age, pira_group, phenot) %>%
  rownames_to_column("sample_id")

# pull also the other useful covariates from the old analysis to be added to the new dataset
data_old_meta <- read.csv("data/dataset_short_all.csv", row.names = "MS_ID") %>%
  # filter only the MS patient
  filter(CTR0_MS1 == 1) %>%
  select(c("Sex_F1_M2","future_PIRA", "previous_PIRA","RELAPSE_act", "MRI_act", "Inflammatory_act")) %>%
  rownames_to_column("sample_id")

data_new_full <- read.csv("data/biomarkers_lod2.csv",row.names = "MS_ID") %>%
  dplyr::select(NfL, GFAP, CHIT,OPN) %>%
  filter(!is.na(OPN)) %>%
  rownames_to_column("sample_id") %>%
  left_join(data_new_meta,by="sample_id") %>%
  left_join(data_old_meta,by="sample_id")

# save the final dataset
data_new_full %>%
  write_csv("data/data_new_full.csv")
