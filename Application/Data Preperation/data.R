library(rstan)
library(frailtypack)

data(readmission)
data <- readmission

# Count the number of patients and number of patients that have passed away in the original data
length(unique(data$id))
length(unique(data$id[which(data$death == 1)]))

table(aggregate(data$event, by = list(data$id), FUN = sum)[,2])
table(aggregate(data$event[which(data$chemo == "Treated")], by = list(data$id[which(data$chemo == "Treated")]), FUN = sum)[,2])
table(aggregate(data$event[which(data$chemo == "NonTreated")], by = list(data$id[which(data$chemo == "NonTreated")]), FUN = sum)[,2])

# Delete patients that have multiple deaths recorded in the data
death_total <- aggregate(data$death, by = list(data$id), FUN = sum)
delete_id_1 <- death_total[which(death_total[,2] > 1),][,1]

# Check the extreme cases of patients with many hospitalizations
data[which(data$id %in% data$id[which(data$enum > 5)]),]
delete_id_2 <- c(274, 350)

# Check whether the variables are changing over time
# Charlson Comorbidity Index (CCI) is changing over time for 96 patients, we do not consider this variable in the analysis
baseline_vars <- c("sex", "dukes", "charlson", "chemo")
var_flags <- aggregate(data[baseline_vars], list(id = data$id),
                       function(x) length(unique(x)) > 1)
names(var_flags)[-1] <- paste0("varies_", baseline_vars)

# Check whether all hospitalizations are recorded 
has_missing_hospitalization <- sapply(tapply(data$enum, data$id, function(x) {
  length(x) != (max(x) - min(x) + 1)
}), identity)
delete_id_3 <- which(has_missing_hospitalization)

# Delete patients
data <- data[-which(data$id %in% c(delete_id_1, delete_id_2, delete_id_3)),]

# Reassign IDs, such that there are no gaps of numbers, this is convenient for stan
old_ids <- sort(unique(data$id))
id_map <- setNames(seq_along(old_ids), old_ids)
data$id <- as.integer(id_map[as.character(data$id)])

# Split the data in a recurrent event dataset and a terminal event dataset
data_recurrent <- data[-which(data$time == 1 & data$death == 1),]
data_terminal <- do.call(rbind, by(data, data$id, function(data) {
  # if any death==1, keep that row
  if (any(data$death == 1)) {
    data[data$death == 1, ][1, ]   # take the first such row if multiple
  } else {
    data[which.max(data$t.stop), ]
  }
}))

# Count the number of patients and number of patients that have passed away in the original data
length(unique(data$id))
length(unique(data$id[which(data$death == 1)]))

table(aggregate(data$event, by = list(data$id), FUN = sum)[,2])
table(aggregate(data$event[which(data$chemo == "Treated")], by = list(data$id[which(data$chemo == "Treated")]), FUN = sum)[,2])
table(aggregate(data$event[which(data$chemo == "NonTreated")], by = list(data$id[which(data$chemo == "NonTreated")]), FUN = sum)[,2])

# chemo = 1 means chemotherapy and chemo = 0 means no chemotherapy
# sex = 1 means female and sex = 0 means male
# create the covariate matrices for the recurrent event data and the terminal event data
data_terminal <- do.call(rbind, by(data, data$id, function(data) {
  # if any death==1, keep that row
  if (any(data$death == 1)) {
    data[data$death == 1, ][1, ]   # take the first such row if multiple
  } else {
    data[which.max(data$t.stop), ]
  }
}))

# duke = 1 corresponds to early-stage 
# chemo = 1 corresponds to chemotherapy 
# sex = 1 corresponds to female 
table(data_terminal$sex)
table(data_terminal$chemo)
table(data_terminal$dukes == "A-B")

# Prepare data in a list for stan 
# Risk factors in matrices
X_recurrent_multiple <- cbind(as.numeric(data_recurrent$sex) - 1, as.numeric(data_recurrent$chemo) - 1, as.numeric(data_recurrent$dukes == "A-B"))
X_terminal <- cbind(as.numeric(data_terminal$sex) - 1, as.numeric(data_terminal$chemo) - 1, as.numeric(data_terminal$dukes == "A-B"))
X_recurrent <- X_terminal

# Total number of subjects, events and variables
N_subject <- length(unique(data_recurrent$id))
N_recurrent <- nrow(data_recurrent)
N_recurrent_covariate <- ncol(X_recurrent)
N_terminal_covariate <- ncol(X_terminal)

# Patient id
id_recurrent <- data_recurrent$id

# Observed gap times for recurrent event data and time to death 
y_recurrent <- data_recurrent$time/365.25
y_terminal <- data_terminal$t.stop/365.25

# Right-censoring for recurrent event data and time to death
censoring_recurrent <- data_recurrent$event
censoring_terminal <- data_terminal$death

# Times for which we compute the number of events and length of grid
renewal_times <- seq(from = 0, to = 10, by = 0.1)
N_times <- length(renewal_times)

# List for stan data
stan_data <- list(
  N_subject = N_subject,
  N_recurrent = N_recurrent,
  id_recurrent = id_recurrent,
  N_recurrent_covariate = N_recurrent_covariate,
  N_terminal_covariate = N_terminal_covariate,
  X_recurrent = X_recurrent,
  X_terminal = X_terminal,
  X_recurrent_multiple = X_recurrent_multiple,
  y_recurrent = y_recurrent,
  y_terminal = y_terminal,
  censoring_recurrent = censoring_recurrent,
  censoring_terminal = censoring_terminal,
  renewal_times = renewal_times,
  N_times = N_times
)

