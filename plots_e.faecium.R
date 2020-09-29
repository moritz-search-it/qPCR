library(tidyverse)

library(modelr)

library(cowplot)
options(na.action = na.warn)

#---values------------

filtering <- function(x){
a <-  van_non_amplificated %>% filter(condition == x)
a$Time <- as_factor(a$Time)
b <- UN_non_amplificated %>% filter(condition == x)
b$Time <- as_factor(b$Time)
c <- van_amplificated%>% filter(condition== x)
c$Time <- as_factor(c$Time)

d <- UN_amplificated %>% filter(condition == x)
d$Time <- as_factor(d$Time)

b <- b %>% mutate(gene= "UN",enrichment=FALSE)
a<- a %>% mutate(gene="van",enrichment=FALSE)
c <- c %>% mutate(gene="van",enrichment=TRUE)
d <- d %>% mutate(gene="UN",enrichment=TRUE)

mod1 <- lm(concentration~Time,data = a)
mod2 <- lm(concentration~Time,data = b)
mod3 <- lm(concentration~Time,data = c)
mod4 <- lm(concentration~Time,data = d)
mod1
grid1 <- a %>%
  data_grid(Time) %>%
  add_predictions(mod1)
grid1
grid2 <- b %>%
  data_grid(Time) %>%
  add_predictions(mod2)
grid3 <- c %>%
  data_grid(Time) %>%
  add_predictions(mod3)
grid4 <- d %>%
  data_grid(Time) %>%
  add_predictions(mod4)

a <- full_join(a,grid1)
a <- a %>% mutate(Time = parse_number(as.character(Time)))
b <- full_join(b,grid2)
b <- b %>% mutate(Time = parse_number(as.character(Time)))
c <- full_join(c,grid3)
c <- c %>% mutate(Time = parse_number(as.character(Time)))
d <- full_join(d,grid4)
d <- d %>% mutate(Time = parse_number(as.character(Time)))

full_join(a,b) %>% full_join(c) %>% full_join(d) %>% 
ggplot(aes(Time,concentration, group=Time,color=enrichment)) + 
  geom_path(aes(alpha=0.5,color=enrichment,group=pred)) + facet_wrap(~gene,ncol=2,scales = "fixed") + 
  geom_point(aes(y=pred,color=enrichment))+ scale_y_log10() +guides(alpha="none") +theme(legend.position = "bottom") +
  labs(x= "Time (hrs)",y= " log 10 concentration" ,title="change in gene concentrations over time reveals wave pattern",
       subtitle = paste("condition inoculum with",x))
g1 <- ggsave2(paste(x,"_change.png"))
return(g1)
}

filtering("ST")
# filtering("inoc") no values in UN
#filtering ("AB") no values in UN
filtering("++")


#----------old graphs not used anymore, use now function filtering

#data_filters <- function(x= "AB"){
 a <-  van_non_amplificated %>% filter(condition == x)
  b <- UN_non_amplificated %>% filter(condition == x)
  
  c <- van_amplificated%>% filter(condition == x)
  d <- UN_amplificated %>% filter(condition == x)
  
  a1 <- ggplot(a,aes(Time,concentration,group=Time)) +  geom_violin(aes(alpha=0.5)) + scale_y_log10() +
    labs(x= "Time (hrs)", title="van gene",
         subtitle = x,
         caption = "no amplification")
  
  b1 <- ggplot(b,aes(Time,concentration,group=Time)) +  geom_violin(aes(alpha=0.5)) + scale_y_log10() +
    labs(x= "Time (hrs)", title="Unique sequence",
         subtitle = x,
         caption = "no amplification")
  
  plot_grid(a1,b1)
 y <-  ggsave2(paste(x,"_non_amplificated.png"))
  
c1 <- ggplot(c,aes(Time,concentration,group=Time)) +  geom_violin(aes(alpha=0.5)) + scale_y_log10() +
    labs(x= "Time (hrs)",title="van gene",
         subtitle =x,
         caption = "amplification")

d1 <- ggplot(d,aes(Time,concentration,group=Time)) + geom_violin(aes(alpha=0.5)) + scale_y_log10() +
    labs(x= "Time (hrs)",title="Unique sequence",
         subtitle = x,
         caption = "amplification")
  plot_grid(c1,d1)
  x <- ggsave2(paste(x,"_amplificated.png"))
  
  return(x)
  return(y)

}

#data_filters()
#data_filters("ST")
#data_filters("inoc")
#data_filters("++")


#-------------statistics
#ggpubr
#compare_means()

#------comparing between strain and strain + Antibiotic condition-------------------------
selecta <- function(x){
  df1 <- as.data.frame(x)
  df1 <- filter(df1, df1$condition != "AB")
  df1 <- filter(df1, df1$condition != "inoc")
  
  return(df1)
}


ggplot(selecta(van_non_amplificated),aes(Time,concentration,group=Time,fill=condition)) +
  geom_point(aes(alpha=0.5, size=5))+ facet_wrap(~condition,ncol=2,scales = "fixed") + scale_y_log10() +
  labs(x="Time in hrs", title="van gene concentration variates over time under both conditions",
       subtitle = "no amplification") + theme(legend.position= "none") 
ggsave2("van_non_amplification.png")

selecta(van_amplificated)
ggplot(selecta(van_amplificated),aes(Time,concentration,group=Time,fill=condition)) +
  geom_point(aes(alpha=0.5, size=5))+ facet_wrap(~condition,ncol=2,scales = "fixed") + scale_y_log10() +
  labs(x="Time in hrs", title="van gene concentration variates over time under both conditions",
       subtitle = "amplification") + theme(legend.position= "none") 
ggsave2("van_amplificated.png")

selecta(UN_non_amplificated)
ggplot(selecta(UN_non_amplificated),aes(Time,concentration,group=Time,fill=condition)) +
  geom_point(aes(alpha=0.5, size=5))+ facet_wrap(~condition,ncol=2,scales = "fixed") + scale_y_log10() +
  labs(x="Time in hrs", title = "Unique sequence concentration variates over time under both conditions",
       subtitle = "no amplification")+ theme(legend.position= "none") 
ggsave2("UN_non_amplificated.png")

selecta(UN_amplificated)
ggplot(selecta(UN_amplificated),aes(Time,concentration,group=Time,fill=condition)) +
  geom_point(aes(alpha=0.5, size=5))+ facet_wrap(~condition,ncol=2,scales = "fixed") + scale_y_log10() +
  labs(x="Time in hrs", title = "Unique sequence concentration variates over time under both conditions",
       subtitle = "amplification")+ theme(legend.position= "none") 
ggsave2("UN_amplificated.png")

#--------ratios
ratio1 <- selecta(ratio_non_amplificated) %>%  mutate(enrichment=FALSE)
ratio2 <- selecta(ratio_amplificated) %>% mutate(enrichment=TRUE)


full_join(ratio1,ratio2) %>% 
ggplot(aes(Time,ratio,group=Time,color=enrichment)) +
  geom_point(aes(alpha=0.5,color=enrichment,size=3)) +
  facet_wrap(~condition,ncol=2,scales = "fixed") + scale_y_log10()  +
  theme(legend.position= "bottom")+ guides(size="none",alpha="none")
ggsave2("ratio.png")


ggplot(selecta(ratio_non_amplificated),aes(Time,ratio,group=Time,fill=condition)) +
  geom_point(aes(alpha=0.5, size=5))+geom_abline(intercept=0,slope=0) +facet_wrap(~condition,ncol=2,scales = "fixed")+ scale_y_log10()+
  theme(legend.position= "none")  
ggsave2("ratio_non_amplificated.png")

ggplot(aes(Time,concentration, group=Time,color=enrichment)) + 
  geom_path(aes(alpha=0.5,color=enrichment,group=pred)) + facet_wrap(~gene,ncol=2,scales = "fixed") + 
  geom_point(aes(y=pred,color=enrichment))+ scale_y_log10() +guides(alpha="none") +theme(legend.position = "bottom") +
  labs(x= "Time (hrs)",y= " log 10 concentration" ,title="change in gene concentrations over time reveals wave pattern",
       subtitle = paste("condition inoculum with",x))
       