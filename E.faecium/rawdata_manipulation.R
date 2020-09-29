library(readr)
library(tidyverse)
datapoints1 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm5_e.faecium/van (10^2 not working)/plate1_datapoints.txt", na = "0", skip=1)
datapoints2 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm5_e.faecium/van (10^2 not working)/plate2_datapoints.txt", na = "0", skip=1)
datapoints3 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm5_e.faecium/van (10^2 not working)/plate3_datapoints.txt", na = "0", skip=1)


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
datapoints_end <- drop_na(datapoints) #drop 5 values incl. 24a_AB_2
View(datapoints_end)


datapoints_end <- datapoints_end%>%
  separate(time, into=c("Time","condition","replicate"), sep= "_") %>% select(-("name"))

samples_van <- sort(unique(datapoints$time))
samples_van

grouped_data <- datapoints_end %>% group_by(Time,condition,replicate)

grouped_data

#for(i in 1:length(grouped_data$Time)){
  #grouped_data$Concentration[i] <- grouped_data$Concentration[i]/100
#}

grouped_data
for(i in 1:length(grouped_data$Time)){
  if(grouped_data$Status[i] == ">")
    grouped_data$Concentration[i] <- 0
  else if(grouped_data$Status[i] == ">,")
   grouped_data$Concentration[i] <-  0
  else(grouped_data$Concentration <- grouped_data$Concentration)
}
grouped_data

#mean calculations
mean_data <- grouped_data %>% summarise(cp=mean(Cp),concentration=mean(Concentration))

amplificated <- mean_data %>% filter(str_detect(Time,'7a|24a')) %>% group_by(Time)
van_amplificated <- amplificated
van_non_amplificated <- anti_join(mean_data,amplificated)
ungroup(van_amplificated)
ungroup(van_non_amplificated)

van_non_amplificated$Time <- van_non_amplificated$Time %>% as.integer()
van_non_amplificated <- arrange(van_non_amplificated,-desc(Time))
van_non_amplificated$Time <- van_non_amplificated$Time %>% as.character()

van_non_amplificated <- van_non_amplificated %>% ungroup(Time) %>% mutate(Time = parse_number(as.character(Time)))

van_non_amplificated <- van_non_amplificated %>% ungroup(condition) %>% mutate(condition= as.factor(condition))

van_non_amplificated$condition <- relevel(van_non_amplificated$condition, "ST")
van_non_amplificated$condition <- relevel(van_non_amplificated$condition, "AB")
van_non_amplificated$condition <- relevel(van_non_amplificated$condition, "inoc")

for (i in 1:length(van_amplificated$Time)) {
  if(van_amplificated$Time[i] == "24a")
    van_amplificated$Time[i] <- "24"
  
  if(van_amplificated$Time[i] == "7a")
    van_amplificated$Time[i] <- "7"
  else( van_amplificated$Time)
  
}



van_amplificated$Time <- van_amplificated$Time %>% as.integer()
van_amplificated <- arrange(van_amplificated,-desc(Time))
van_amplificated$Time <- van_amplificated$Time %>% as.character()

van_amplificated <- van_amplificated %>% ungroup(Time) %>% mutate(Time = parse_number(as.character(Time)))

van_amplificated <- van_amplificated %>% ungroup(condition) %>% mutate(condition= as.factor(condition))

van_amplificated$condition <- relevel(van_amplificated$condition, "ST")
van_amplificated$condition <- relevel(van_amplificated$condition, "AB")
van_amplificated$condition <- relevel(van_amplificated$condition, "inoc")
#################################
##-------------für UN ----------------------------------------------------------
##################################
datapoints1 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm5_e.faecium/UN (10^8 not working)/plate1_datapoints.txt", na = "0", skip=1)
datapoints2 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm5_e.faecium/UN (10^8 not working)/plate2_datapoints.txt", na = "0", skip=1)
datapoints3 <- read_table2("/Users/moritzherrmann/Bioinformatic_MA/Data/qPCR/BatchFerm5_e.faecium/UN (10^8 not working)/plate3_datapoints.txt", na = "0", skip=1)


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
datapoints_end <- drop_na(datapoints) #drop 3 values all of 0_inoc3
View(datapoints_end)


datapoints_end <- datapoints_end%>%
  separate(time, into=c("Time","condition","replicate"), sep= "_") %>% select(-("name"))

samples <- sort(unique(datapoints$time))
samples

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

ratio_amplificated <- left_join(van_amplificated,UN_amplificated,by= c("Time","condition","replicate")) %>% 
  mutate(ratio = concentration.x/concentration.y) 
ratio_amplificated <- drop_na(ratio_amplificated)
ratio_amplificated <- arrange(ratio_amplificated,-desc(Time))
ratio_non_amplificated <- full_join(van_non_amplificated,UN_non_amplificated, by=c("Time","condition","replicate")) %>% 
  mutate(ratio = concentration.x/concentration.y)
ratio_non_amplificated <- drop_na(ratio_non_amplificated)
ratio_non_amplificated <- arrange(ratio_non_amplificated,-desc(Time))

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
