rm(list=ls())
gc()
library(tidyverse)

# read in sign-up survey data
data <- read.csv("./all_apps_wide_2026-02-26.csv")

# select session: 
session_code = "hw7cmmay"
data2 <- data[data$session.code == session_code, ]

# select eligible participants (who signed up, gave consent, passed checks) 
# who confirmed their participation
data2 %>%
  filter(return.1.player.confirm == 1) %>%
  select(participant.label, participant.role) -> data3

# this needs to be included in the experiment repository (on arrival, roles are assigned based on ID)
write.csv(data3, file = "participants.csv")

# to invite participants via Prolific
writeLines(data3$participant.label, "participants_ids.txt")



