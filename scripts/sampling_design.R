############ Prepare complete dataset 

##################################################### Australia
rm(list=ls())

julie <- Sys.info()["nodename"] == "GCI-CY9ZK72"

if (julie){
  .libPaths("C:\\Julie\\Rpkgs")
}

if (julie){
  .libPaths("C:\\Julie\\Rpkgs")
}

source("R/packages.R")


db <- dbConnect(dbDriver("MySQL"),host = '203.101.226.180' , user = "julie", password = "spatial", dbname = "grr")

######### Step 1. Query from the database

query="SELECT s.transectid, YEAR(s.start_datetime) AS year, t.reef_name,s.reef_type,s.sub_region, s.transects_reef_section,
s.transects_length,s.dive_visibility, s.dive_temperature,i.depth,i.temperature as image_temperature, i.surveyid, 
c.* FROM australiacovers_quadrats c, quadrats q, images i, surveys s,transects t WHERE c.id=q.id  AND q.imageid=i.id 
AND s.id=i.surveyid AND t.id=s.transectid";

sql = sprintf(query)

######### Step 2. Send the query
t = dbSendQuery(db,sql)
dbColumnInfo(t)
t = fetch(t, n = -1)

######## Step 3. Stop the query
dbClearResult(dbListResults(db)[[1]])

######## # ## DISCONNECT FROM DATABASE
dbDisconnect(db)

################################### SELECT RADNOMLY ONE QUADRAT PER IMAGE
############ Add number of quadrats per images

t_sorted <- t %>%
  mutate(Transectid = substr(id, 0, 5)) %>%
  mutate(Image_ID = substr(id, 1, 9)) %>%
  mutate(Quadrat_id = substr(id, 10, 11))

year_keep <- c("2012")
t_sorted <- t_sorted %>%
  filter(year %in% year_keep) %>%
  group_by(Image_ID) %>%
  sample_n(1)


### Test

t_check <- t_sorted %>%
  group_by(Image_ID) %>%
  tally() %>%
  filter(n > 1)

t_check

################################### GATHER DATABLE 

Rec_t<-t_sorted%>%
  gather(key = label, Proportion,
         16:ncol(t))

############################################################################################ Merge with universal labelset
if (julie){
  labelset <- read.csv("C:\\Users\\uqjverce\\Dropbox\\GCI\\Catlin\\Research projects\\Global\\data\\labelset_full_edited.csv")
}else{
  labelset <- read.csv("/media/shallowreefs/papers/Global/data/labelset_full_edited.csv") 
}

############ Merge with labelset
lab <- labelset %>%
  dplyr:::select(label, global_label, fg_global) %>%
  mutate_if(is.factor, as.character) %>%
  group_by(label) %>%
  slice(1)

Rec <- inner_join(Rec_t, lab) %>%
  data.frame()

################### Mean per transect and global_label
to_remove <- c("NAN", "BL")

############################################################## Rescale proportion

Rec_rescale <- Rec %>%
  filter(!global_label %in% to_remove) %>%
  group_by(Image_ID) %>% mutate(Prop_check = sum(Proportion)) %>%
  ungroup() %>%
  mutate(New_prop = Proportion/Prop_check) %>%
  data.frame()

Rec_rescale$New_prop[is.na(Rec_rescale$New_prop)] <- 0
#Rec_rescale_check<-Rec_rescale%>%group_by(splitting.var)%>%mutate(Sum_new_prop=sum(New_prop))%>%dplyr:::select(splitting.var,Sum_new_prop)%>%data.frame()

Rec_morph <- Rec_rescale %>%
  group_by(Image_ID, global_label) %>% 
  mutate(Prop_morph = sum(New_prop)) %>%
  slice(1) %>%
  dplyr:::select(., -Proportion, -label, -New_prop) %>% 
  data.frame()

Rec_dom <- Rec_morph

############################# Get number of splitting.var per surveyid

Rec.numb <- Rec_dom %>%
  group_by(surveyid, Image_ID) %>%
  tally() %>%
  dplyr::select(., -n) %>%
  ungroup() %>%
  group_by(surveyid) %>%
  tally() %>%
  filter(n > 150) %>%
  dplyr::select(surveyid) %>%
  data.frame()


Rec_dom <- Rec_dom %>%
  filter(surveyid %in% Rec.numb$surveyid)

#### How many images per survey?
Rec.numb2 <- Rec_dom %>%
  group_by(surveyid, Image_ID) %>%
  tally() %>%
  dplyr::select(., -n) %>%
  ungroup() %>%
  group_by(surveyid) %>%
  tally()

ggplot(Rec.numb2, aes(x = n)) + geom_histogram()

########################################## Nest the longitude and latitude variables 
Rec.nest <- Rec_dom %>%
  dplyr::select(surveyid, Image_ID, Prop_morph) %>%
  nest(Prop_morph, .key = coords)
Rec.nest

##################################################### Pick-up transects randomly

survey_random <- unique(Rec.nest$surveyid)

##################################################### Spearman function
fun_run <- function(coords.x, coords.y) {
  cor(coords.x, coords.y, method = "spearman")
}

##################################################### Extract roots function
extractRoots <- function(parms){
  a <- parms[1]
  b <- parms[2]
  c <- parms[3]
  
  x1 <- (-b + sqrt(b^2 - 4*a*c))/(2*a)
  x2 <- (-b - sqrt(b^2 - 4*a*c))/(2*a)
  
  return(x1)
}

################### Save dataframe with minumum of image number for each transect

Rec_final <- data.frame(matrix(NA, length(survey_random), 4))
colnames(Rec_final) <- c("SurveyID", "Min_Number", "Nbre_combination", "Region")

Rec_final$Region <- rep("Australia", length(survey_random))

###################### Start loop 

for (j in 1:length(survey_random)){
  
  test <- Rec.nest %>%
    filter(surveyid == survey_random[j]) %>%
    droplevels()
  
  Rec.spear <- test %>%
    mutate(Image_ID1 = Image_ID) %>%
    dplyr:::select(.,-coords) %>%
    do(complete(., Image_ID, Image_ID1, fill = list(surveyid = .$surveyid))) %>%
    filter(Image_ID != Image_ID1) %>%
    left_join(Rec.nest, by = "Image_ID") %>%
    left_join(Rec.nest, by = c("Image_ID1" = "Image_ID")) %>%
    ungroup() %>%
    mutate(Spear = map2_dbl(.x = coords.x, .y = coords.y, .f = fun_run)) 
  
  
  ###################### See how much sub-transect is needed to reach a plateau
  Nbre_combination <- c(10, 100, 200, 500, 1000, 5000, 10000, 20000)
  
  Rec.spear.list <- list()
  
  for (i in 1:length(Nbre_combination)){
    Rec.spear.list[[i]] <- Rec.spear %>%
      sample_n(Nbre_combination[i]) %>%
      mutate(Nbre_combination = rep(Nbre_combination[i],
                                    Nbre_combination[i])) %>%
      mutate(Image_number = length(setdiff(unique(Image_ID), 
                                         unique(Image_ID1))))
    
  }
  Rec.spear.collapsed <- do.call(what = rbind, args = Rec.spear.list)
  
  # ggplot(Rec.spear.collapsed,aes(x=Nbre_combination,y=Spear,group=Nbre_combination))+geom_boxplot()+
  #   stat_summary(fun.y=mean, colour="darkred", geom="point", 
  #                shape=18, size=3)+
  #   scale_x_log10()+
  #   ylab("Spearman Correlation")+xlab("Number of image pairs")
  
  Rec.spear.list.summary <- list()
  
  for (i in 1:length(Nbre_combination)){
    Rec.spear.list.summary[[i]] <- Rec.spear %>%
      sample_n(Nbre_combination[i]) %>%
      filter(!is.na(Spear)) %>%
      summarise(Mean = mean(Spear),
                SD = sd(Spear), 
                N = Nbre_combination[i]) %>%
      mutate(SE = SD/sqrt(N)) %>%
      data.frame()
  }
  
  Rec.spear.collapsed.summary <- do.call(what = rbind, args = Rec.spear.list.summary)
  
  Min_sampling <- Rec.spear.collapsed.summary %>%
    filter(SE < 0.01)
  parms <- c(1, -1, -Min_sampling$N[1]*2)
  
  Nbre_image<-ceiling(extractRoots(parms))
  
  p1 <- ggplot(Rec.spear.collapsed.summary,
               aes(x = Nbre_combination, y = Mean)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin = Mean-1.96*SE, ymax = Mean+1.96*SE), alpha = 0.2) +
    ylim(-1,1) +
    ylab("Spearman Correlation") + 
    xlab("Number of image pairs") +
    ggtitle(survey_random[j]) +
    geom_label(aes(x = Nbre_combination[4],
                   y = -0.8,
                   label = paste("The minimum image number is", 
                                 Nbre_image, sep = " ")), size = 5.3) +
    scale_x_log10() +
    theme(axis.text.x = element_text(size = 10),
          legend.title = element_text(colour = "black", size = 110, face = "bold"),
          legend.text = element_text(colour = "black", size = 15),
          plot.title = element_text(hjust = 0.5, vjust = 2, size = 15, face = "bold"),
          axis.text.y = element_text(size = 13), 
          axis.title.y = element_text(size = 15),
          axis.title.x = element_text(size = 15),
          strip.text = element_text(size = 14),
          strip.background = element_rect(fill = "gray10"))
  
  
  if (julie){
    dir <- "C:\\Users\\uqjverce\\Dropbox\\GCI\\Catlin\\Research projects\\Global\\Reassembly_study\\Sampling_effect\\Australia"
  }else{
    dir <- "/media/shallowreefs/papers/Global/Reassembly_study/Sampling_effect/Australia2"
  }
  name_plot <- paste("Sampling_effect_images", survey_random[j], ".pdf", sep = "")
  ggsave(plot = p1, file = paste(dir, name_plot, sep = "/"))
  
  Rec_final[j, 1] <- survey_random[j]
  Rec_final[j, 2] <- Nbre_image
  Rec_final[j, 3] <- choose(Nbre_image, 2)
}

if (julie){
  write.csv(Rec_final,
            file = "C:\\Users\\uqjverce\\Dropbox\\GCI\\Catlin\\Research projects\\Global\\Reassembly_study\\Sampling_effect\\Sampling_effect_Australia.csv",row.names = F)
}else{
  save(Rec_final,
       file = "/media/shallowreefs/papers/Global/Reassembly_study/Sampling_effect/Sampling_effect_Australia2")
}