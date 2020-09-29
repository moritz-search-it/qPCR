library(readr)
library(tidyverse)
datapoints1 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm4_e.coli/bla/datapoints1.txt", na = "0", skip=1)
datapoints2 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm4_e.coli/bla/datapoints2.txt", na = "0", skip=1)
datapoints3 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm4_e.coli/bla/datapoints3.txt", na = "0", skip=1)


select(datapoints1, c("Name","Cp","Concentration","Status"))

all_datapoints <- datapoints1 %>% 
  full_join(datapoints2) %>% 
  full_join(datapoints3) %>% 
  select( , c("Name","Cp","Concentration","Status"))

#----------------------------------------------------------------------------------------------------------------------------------------


##Perform data transformation, calculating mean of replicates in Cp and concentration

datapoints<- all_datapoints %>% separate(Name,into = c("name","time"),sep = "t")
datapoints$Status<- replace_na(datapoints$Status,"0")
View(datapoints)
datapoints_end <- drop_na(datapoints) #drop 5 values incl. 24a_AB_2 in e.faecium // drop 3 values incl. 48_++_3 in e.coli
View(datapoints_end)


datapoints_end <- datapoints_end%>%
  separate(time, into=c("Time","condition","replicate"), sep= "_") %>% select(-("name"))

samples_bla <- sort(unique(datapoints$time))
samples_bla

grouped_data <- datapoints_end %>% group_by(Time,condition,replicate)

grouped_data

#for(i in 1:length(grouped_data$Time)){
  #grouped_data$Concentration[i] <- grouped_data$Concentration[i]/100
#}

grouped_data
for(i in 1:length(grouped_data$Time)){
  if(grouped_data$Status[i] == "E")
    grouped_data$Concentration[i] <- 0
  else if(grouped_data$Status[i] == ">,")
   grouped_data$Concentration[i] <-  0
  else(grouped_data$Concentration <- grouped_data$Concentration)
}
grouped_data

#mean calculations
mean_data <- grouped_data %>% summarise(cp=mean(Cp),concentration=mean(Concentration))

amplificated <- mean_data %>% filter(str_detect(Time,'7a|24a')) %>% group_by(Time)
bla_amplificated <- amplificated
bla_non_amplificated <- anti_join(mean_data,amplificated)
ungroup(bla_amplificated)
ungroup(bla_non_amplificated)

bla_non_amplificated <- bla_non_amplificated %>% ungroup(Time) %>% mutate(Time = parse_number(as.character(Time)))

bla_non_amplificated <- bla_non_amplificated %>% ungroup(condition) %>% mutate(condition= as.factor(condition))

bla_non_amplificated$condition <- relevel(bla_non_amplificated$condition, "ST")
bla_non_amplificated$condition <- relevel(bla_non_amplificated$condition, "AB")
bla_non_amplificated$condition <- relevel(bla_non_amplificated$condition, "inoc")

for (i in 1:length(bla_amplificated$Time)) {
  if(bla_amplificated$Time[i] == "24a")
    bla_amplificated$Time[i] <- "24"
  
  if(bla_amplificated$Time[i] == "7a")
    bla_amplificated$Time[i] <- "7"
  else( bla_amplificated$Time)
  
}



bla_amplificated <- bla_amplificated %>% ungroup(Time) %>% mutate(Time = parse_number(as.character(Time)))

bla_amplificated <- bla_amplificated %>% ungroup(condition) %>% mutate(condition= as.factor(condition))

bla_amplificated$condition <- relevel(bla_amplificated$condition, "ST")
bla_amplificated$condition <- relevel(bla_amplificated$condition, "AB")
bla_amplificated$condition <- relevel(bla_amplificated$condition, "inoc")
#################################
##-------------für UN ----------------------------------------------------------
##################################
datapoints1 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm4_e.coli/UN/datapoints1.txt", na = "0", skip=1)
datapoints2 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm4_e.coli/UN/datapoints2.txt", na = "0", skip=1)
datapoints3 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm4_e.coli/UN/datapoints3.txt", na = "0", skip=1)


select(datapoints1, c("Name","Cp","Concentration","Status"))

all_datapoints <- datapoints1 %>% 
  full_join(datapoints2) %>% 
  full_join(datapoints3) %>% 
  select( , c("Name","Cp","Concentration","Status"))

#----------------------------------------------------------------------------------------------------------------------------------------


##Perform data transformation, calculating mean of replicates in Cp and concentration

datapoints<- all_datapoints %>% separate(Name,into = c("name","time"),sep = "t")
datapoints$Status<- replace_na(datapoints$Status,"0")
View(datapoints)
datapoints_end <- drop_na(datapoints) #drop 3 values all of 0_inoc3 in e.faecium// no drops in e.coli
View(datapoints_end)


datapoints_end <- datapoints_end%>%
  separate(time, into=c("Time","condition","replicate"), sep= "_") %>% select(-("name"))

UN_samples <- sort(unique(datapoints$time))
UN_samples

grouped_data <- datapoints_end %>% group_by(Time,condition,replicate)

grouped_data

for(i in 1:length(grouped_data$Time)){
if(grouped_data$Status[i] == "E")
 grouped_data$Concentration[i] <- 0
  else if(grouped_data$Status[i] == ">,")
    grouped_data$Concentration[i] <-  0
  else(grouped_data$Concentration <- grouped_data$Concentration)
    }

grouped_data
#mean calculations
mean_data <- grouped_data %>% summarise(cp=mean(Cp),concentration=mean(Concentration))
mean_data
#semi_join(mean_data,grouped_data,by=c("Time","replicate","condition"))

amplificated <- mean_data %>% filter(str_detect(Time,'7a|24a')) %>% group_by(Time)
non_amplificated <- anti_join(mean_data,amplificated)

UN_amplificated <- amplificated
UN_non_amplificated <- non_amplificated
ungroup(UN_amplificated)
ungroup(UN_non_amplificated)

UN_non_amplificated$Time <- UN_non_amplificated$Time %>% as.integer()
UN_non_amplificated <- arrange(UN_non_amplificated,-desc(Time))
UN_non_amplificated$Time <- UN_non_amplificated$Time %>% as.character()
#UN_non_amplificated$Time <- UN_non_amplificated$Time %>% as.character()

UN_non_amplificated <- UN_non_amplificated %>% ungroup(Time) %>% mutate(Time = parse_number(as.character(Time)))

UN_non_amplificated <- UN_non_amplificated %>% ungroup(condition) %>% mutate(condition= as.factor(condition))

UN_non_amplificated$condition <- relevel(UN_non_amplificated$condition, "ST")
UN_non_amplificated$condition <- relevel(UN_non_amplificated$condition, "AB")
UN_non_amplificated$condition <- relevel(UN_non_amplificated$condition, "inoc")

for (i in 1:length(UN_amplificated$Time)) {
  if(UN_amplificated$Time[i] == "24a")
    UN_amplificated$Time[i] <- 24
  
  if(UN_amplificated$Time[i] == "7a")
    UN_amplificated$Time[i] <- 7
  else( UN_amplificated$Time)
  
}

UN_amplificated$Time <- UN_amplificated$Time %>% as.integer()
UN_amplificated <- arrange(UN_amplificated,-desc(Time))
UN_amplificated$Time <- UN_amplificated$Time %>% as.character()

UN_amplificated <- UN_amplificated %>% ungroup(Time) %>% mutate(Time = parse_number(as.character(Time)))

UN_amplificated <- UN_amplificated %>% ungroup(condition) %>% mutate(condition= as.factor(condition))

UN_amplificated$condition <- relevel(UN_amplificated$condition, "ST")
UN_amplificated$condition <- relevel(UN_amplificated$condition, "AB")
UN_amplificated$condition <- relevel(UN_amplificated$condition, "inoc")

# ratios


ratio_amplificated <- left_join(bla_amplificated,UN_amplificated,by= c("Time","condition","replicate")) %>% 
  mutate(ratio = concentration.x/concentration.y) 
ratio_amplificated <- drop_na(ratio_amplificated)
ratio_non_amplificated <- full_join(bla_non_amplificated,UN_non_amplificated, by=c("Time","condition","replicate")) %>% 
  mutate(ratio = concentration.x/concentration.y)
ratio_non_amplificated <- drop_na(ratio_non_amplificated)

for(i in 1:length(ratio_non_amplificated$Time)){
  if(ratio_non_amplificated$ratio[i] == Inf)
    ratio_non_amplificated$ratio[i] <- 0
  else (ratio_non_amplificated$ratio)
}
#alternativ to for if loop:
#ratio_non_amplificated %>% mutate(ratio=ifelse(ratio==Inf,"0",ratio))

for(i in 1:length(ratio_amplificated$Time)){
  if(ratio_amplificated$ratio[i] == Inf)
    ratio_amplificated$ratio[i] <- 0
  else (ratio_amplificated$ratio)
}


# remove
rm(all_datapoints)
rm(datapoints)
rm(datapoints_end)
rm(datapoints1,datapoints2,datapoints3)
rm(grouped_data)
rm(mean_data)
rm(amplificated,non_amplificated)

