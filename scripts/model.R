#######  Boral model 

rm(list=ls())

julie <- Sys.info()["nodename"] == "SEF-PA00130783"

if (julie){
  .libPaths("C:\\Julie\\Rpkgs")
}


library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(mvabund)
library(boral)
library(tidyquant)
library(magrittr)
library(ggraph)
library(igraph)
library(tidygraph)
library(codatools)
library(forcats)


# Read the data
df2<-read.csv('data/reef_data.csv') %>%
  mutate(Count=ceiling(Mean*50))%>%mutate(Count_log=log(Count+1))%>%mutate(region=rep("Australia"))%>%mutate(Points=50)%>%
  filter(year=="2012")%>%filter(!global_label=="OTHER")

df.full<-df2

df.full.coord<-df.full%>%group_by(reef_name)%>%filter(row_number()==1)%>%
  dplyr::select(region,reef_name,year,lng,lat)

df.spread<-df.full%>%
  dplyr::select(reef_name,region,sub_region,Bleaching.cat,Density.cat,Cyclone.cat,global_label,Count)%>%
  gather(variable, value, -(reef_name:global_label)) %>%
  unite(temp, global_label, variable) %>%
  spread(temp, value)

df.spread<-df.spread%>%replace(., is.na(.), 0)
df.spread<-df.spread%>%arrange(sub_region)

############### responses
y<-ceiling(df.spread[,7:ncol(df.spread)])%>%as.matrix()
colnames(y)<-str_remove(colnames(y),'_Count')

############### covariates
df.spread$sub_region_quant<-as.numeric(df.spread$sub_region)
df.spread$Bleaching.cat<-as.numeric(df.spread$Bleaching.cat)
df.spread$Cyclone.cat<-as.numeric(df.spread$Cyclone.cat)
df.spread$Density.cat<-as.numeric(df.spread$Density.cat)


#X<-df.spread[,c(ncol(df.spread),4, 5, 6)]%>% mutate_all(.,as.integer)
X0<-df.spread[,c(4, 5, 6)]%>% mutate_all(.,as.integer)

n <- nrow(y)
p <- ncol(y)


################################# Formulation model: 2LVs, all, random effects on the reef
id.row<-1:nrow(y)%>%as.matrix

# Use this params for testing 
mcmc_control <- list(n.burnin = 300, n.iteration = 1000,
                     n.thin = 2)

# Use this control params for running the final model 
#mcmc_control <- list(n.burnin = 3000, n.iteration = 90000,
#                     n.thin = 50)

model.list<- boral(y,X=X0, family = "negative.binomial", row.eff = "random",row.ids=id.row,
                   mcmc.control = mcmc_control,save.model=TRUE,calc.ics = TRUE,
                   lv.control = list(num.lv = 2, type = "independent", distmat = NULL),
                   prior.control = list(type = c("normal","normal","normal","uniform"),
                                        hypparams = c(10, 10, 10, 10)))

############################### Variance partitioning

varpart<-calc.varpart(model.list, groupX=c(1,2,3,4))

varpart.lv<-varpart$varpart.lv%>%data.frame%>%cbind(names(varpart$varpart.lv))%>%mutate(Origin=rep("LVs"))%>%
  dplyr::select(2,1,3)

colnames(varpart.lv)<-c("Functional_group","Variance","Origin")

varpart.x<-varpart$varpart.X%>%data.frame%>%gather(key="Functional_group",value = "Variance")%>%
  mutate(Origin=rep(c("Intercept",colnames(model.list$X)),length(names(varpart$varpart.lv))))

colnames(varpart.x)<-c("Functional_group","Variance","Origin")


varpart.rand<-varpart$varpart.row%>%data.frame%>%cbind(names(varpart$varpart.row))%>%mutate(Origin=rep("Random effect"))%>%
  dplyr::select(2,1,3)

colnames(varpart.rand)<-c("Functional_group","Variance","Origin")

varpart.tot<-rbind(varpart.x,varpart.lv,varpart.rand)%>%data.frame()
varpart.tot$Functional_group<-str_remove(varpart.tot$Functional_group,'_Count')
varpart.tot$Functional_group<-str_replace(varpart.tot$Functional_group,"[.]","-")
varpart.tot<-varpart.tot%>%filter(!Origin=="Intercept")



## Add general group 
labelset<-read.csv("data/labelset_full_edited.csv")%>%
  group_by(global_label)%>%filter(row_number()==1)%>%
  dplyr::select(global_label,fg_global)

colnames(labelset)<-c("Functional_group","Group")

varpart.tot<-inner_join(varpart.tot,labelset)%>%arrange(Origin,Variance)

varpart.tot<-varpart.tot%>%mutate(Origin=ifelse(Origin=="Bleaching.cat","Bleaching occurence",
                                                ifelse(Origin=="Cyclone.cat","Cyclone exposure", 
                                                       ifelse(Origin=="Density.cat","Human Density",Origin))))

varpart.tot$Origin<-as.factor(varpart.tot$Origin)
varpart.tot$Origin<-factor(varpart.tot$Origin,levels=c("Bleaching occurence","Cyclone exposure",
                                                       "Human Density","LVs","Random effect"))


varpart.tot$Group<-factor(varpart.tot$Group,levels=c("Hard corals","Soft Coral","Other Invertebrates","Algae","Other"))


############################### Variance partitioning 2

varpart2<-calc.varpart(model.list)

varpart2.lv<-varpart2$varpart.lv%>%data.frame%>%cbind(names(varpart2$varpart.lv))%>%mutate(Origin=rep("LVs"))%>%
  dplyr::select(2,1,3)

colnames(varpart2.lv)<-c("Functional_group","Variance","Origin")

varpart2.x<-varpart2$varpart.X%>%data.frame%>%cbind(names(varpart2$varpart.lv))%>%mutate(Origin=rep("Environmnental predictors"))%>%
  dplyr::select(2,1,3)
colnames(varpart2.x)<-c("Functional_group","Variance","Origin")


varpart2.rand<-varpart2$varpart.row%>%data.frame%>%cbind(names(varpart2$varpart.row))%>%mutate(Origin=rep("Random effect"))%>%
  dplyr::select(2,1,3)

colnames(varpart2.rand)<-c("Functional_group","Variance","Origin")

varpart2.tot<-rbind(varpart2.x,varpart2.lv,varpart2.rand)%>%data.frame()
varpart2.tot$Functional_group<-str_remove(varpart2.tot$Functional_group,'_Count')
varpart2.tot$Functional_group<-str_replace(varpart2.tot$Functional_group,"[.]","-")
varpart2.tot<-varpart2.tot%>%filter(!Origin=="Intercept")



## Add general group 
labelset<-read.csv("data/labelset_full_edited.csv")%>%
  group_by(global_label)%>%filter(row_number()==1)%>%
  dplyr::select(global_label,fg_global)

colnames(labelset)<-c("Functional_group","Group")

varpart2.tot<-inner_join(varpart2.tot,labelset)%>%arrange(Origin,Variance)

varpart2.tot$Origin<-as.factor(varpart2.tot$Origin)

varpart2.tot$Group<-factor(varpart2.tot$Group,levels=c("Hard corals","Soft Coral","Other Invertebrates","Algae","Other"))


#p.varpart2<-ggplot(data=varpart2.tot, aes(x=fct_reorder(Functional_group,as.numeric(Group)),y=Variance, fill = Origin)) +
#  geom_bar(position = "stack",stat = "identity",alpha=.4)+ylab("Partition of variance (%)")+xlab("Functional group")+
#  theme_tq()+
#  theme(axis.text.x = element_text(size=10,angle=30,hjust=1),
#        legend.title = element_text(colour = "black", size = 13, face="bold"),
#        legend.text = element_text(colour = "black", size = 15),
#        plot.title = element_text(hjust=0.5,vjust=2,size=15, face="bold"),
#        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
#       strip.text=element_text(size=14),strip.background = element_rect(fill="white"))+
#  scale_fill_manual(values=c("darkcyan",
#                             palette_light()[[2]],palette_light()[[7]]))

#p.varpart2    

#ggsave(plot = p.varpart2,width=10, height=6, file = "C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Australia2//Variance_partioning//Var_2cat.pdf")

##################### Covariate effect

Functional_group<-rownames(model.list$X.coefs.mean)
Functional_group<-str_remove(Functional_group,'_Count')


X.coefs.<-cbind(model.list$X.coefs.mean%>%data.frame()%>%gather(key=Covariate, value = Mean)%>%mutate(Functional_group=rep(Functional_group,ncol(model.list$X.coefs.mean))),
                model.list$hpdintervals$X.coefs[,,type = "lower"]%>%data.frame()%>%gather(key=Covariate, value = Lower)%>%dplyr::select(Lower),
                model.list$hpdintervals$X.coefs[,,type = "upper"]%>%data.frame()%>%gather(key=Covariate, value = Upper)%>%dplyr::select(Upper))

X.coefs.<-X.coefs.%>%mutate(Cov=ifelse(Covariate=="Bleaching.cat","Bleaching occurence",
                                       ifelse(Covariate=="Cyclone.cat","Cyclone exposure", 
                                              ifelse(Covariate=="Density.cat","Human Density",
                                                     ifelse(Covariate=="sub_region_quant","Subregion",NA)))))

X.coefs.<-X.coefs.%>%mutate(Sig=ifelse(Lower<0 & Upper<=0 | Lower>=0 & Upper>0,1,0))  


## Latent variables effect
mcmc.raw<-get.mcmcsamples(model.list)

#save(mcmc.raw,file="C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Indian_ocean2//MCMC_outputs//MCMC.Rdata")

mcmc.nb<-coda_df(mcmc.raw)%>%
  gather(key=Parameter,value = value) 


### Grap LV1
lv1 <- grep("^lv.c.+1]$", unique(mcmc.nb$Parameter), value=TRUE)

mcmc.lv1<-mcmc.nb%>%filter(Parameter%in%lv1)%>%
  mutate(Functional_group=rep(Functional_group,each=dim(mcmc.raw)[[1]]))

## Sumarize

quantt<-mcmc.lv1%>%group_by(Functional_group)%>% 
  do(data.frame(t(quantile(.$value, probs = c(0.025, 0.975)))))

mean<-mcmc.lv1%>%group_by(Functional_group)%>% 
  do(data.frame(t(mean(.$value))))

sum_lv1<-inner_join(mean,quantt)
colnames(sum_lv1)<-c("Functional_group","Mean","Lower","Upper")

sum_lv1<-sum_lv1%>%mutate(Sig=ifelse(Lower<0 & Upper<=0 | Lower>=0 & Upper>0,1,0))

### Grap lv2
lv2 <- grep("^lv.c.+2]$", unique(mcmc.nb$Parameter), value=TRUE)

mcmc.lv2<-mcmc.nb%>%filter(Parameter%in%lv2)%>%
  mutate(Functional_group=rep(Functional_group,each=dim(mcmc.raw)[[1]]))

## Sumarize

quantt<-mcmc.lv2%>%group_by(Functional_group)%>% 
  do(data.frame(t(quantile(.$value, probs = c(0.025, 0.975)))))

mean<-mcmc.lv2%>%group_by(Functional_group)%>% 
  do(data.frame(t(mean(.$value))))

sum_lv2<-inner_join(mean,quantt)
colnames(sum_lv2)<-c("Functional_group","Mean","Lower","Upper")

sum_lv2<-sum_lv2%>%mutate(Sig=ifelse(Lower<0 & Upper<=0 | Lower>=0 & Upper>0,1,0))

xmin<-min(X.coefs.$Lower,sum_lv1$Lower,sum_lv2$Lower)
xmax<-max(X.coefs.$Upper,sum_lv1$Upper,sum_lv2$Upper)

X.coefs.<-inner_join(X.coefs.,labelset)
X.coefs.$Group<-factor(X.coefs.$Group,levels=c("Hard corals","Soft Coral","Other Invertebrates","Algae","Other"))

sum_lv1<-inner_join(sum_lv1,labelset)
sum_lv1$Group<-factor(sum_lv1$Group,levels=c("Hard corals","Soft Coral","Other Invertebrates","Algae","Other"))

sum_lv2<-inner_join(sum_lv2,labelset)
sum_lv2$Group<-factor(sum_lv2$Group,levels=c("Hard corals","Soft Coral","Other Invertebrates","Algae","Other"))


p_covariate<-ggplot()+geom_point(data=X.coefs.,aes(y=Mean,x=fct_reorder(Functional_group,as.numeric(Group)),col=as.factor(Sig)))+
  geom_errorbar(data=X.coefs.,aes(x=Functional_group,ymin=Lower,ymax=Upper,col=as.factor(Sig)),width=0.1)+
  facet_wrap(~Cov)+coord_flip()+
  theme(axis.text.x = element_text(size=13),legend.position="none",
        panel.grid = element_blank(),plot.title = element_text(hjust=0,vjust=2,size=18, face="bold"),
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
        strip.text=element_text(size=13),strip.background = element_rect(fill="white"))+ 
  scale_color_brewer(type = "qual", palette = "Paired")+ylim(xmin,xmax)+
  geom_hline(yintercept = 0,linetype="dashed")+xlab("")+ylab("Estimates")+xlab("Functional groups")



p_lv1<-ggplot()+geom_point(data=sum_lv1,aes(y=Mean,x=fct_reorder(Functional_group,as.numeric(Group)),col=as.factor(Sig)))+
  geom_errorbar(data=sum_lv1,aes(x=Functional_group,ymin=Lower,ymax=Upper,col=as.factor(Sig)),width=0.1)+
  coord_flip()+
  theme(axis.text.x = element_text(size=13),legend.position="none",
        panel.grid = element_blank(),plot.title = element_text(hjust=0,vjust=2,size=18, face="bold"),
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
        strip.text=element_text(size=14),strip.background = element_rect(fill="gray98"))+ 
  scale_color_brewer(type = "qual", palette = "Paired")+ylim(xmin,xmax)+
  geom_hline(yintercept = 0,linetype="dashed")+xlab("")+ylab("LV1")


p_lv2<-ggplot()+geom_point(data=sum_lv2,aes(y=Mean,x=fct_reorder(Functional_group,as.numeric(Group)),col=as.factor(Sig)))+
  geom_errorbar(data=sum_lv2,aes(x=Functional_group,ymin=Lower,ymax=Upper,col=as.factor(Sig)),width=0.1)+
  coord_flip()+
  theme(axis.text.x = element_text(size=13),legend.position="none",
        panel.grid = element_blank(),plot.title = element_text(hjust=0,vjust=2,size=18, face="bold"),
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
        strip.text=element_text(size=14),strip.background = element_rect(fill="gray98"))+ 
  scale_color_brewer(type = "qual", palette = "Paired")+ylim(xmin,xmax)+
  geom_hline(yintercept = 0,linetype="dashed")+xlab("")+ylab("LV2")

p_cat1<-plot_grid(p_lv1,p_lv2,nrow=2,align = 'vh',
                  hjust = -1,labels = c("B","C"))


p_cat<-plot_grid(p_covariate, p_cat1,
                 align = 'vh',labels=c("A",""),
                 hjust = -1,
                 nrow = 1,rel_widths = c(1/2,1/3))

ggsave(plot = p_cat,width=13, height=6, file = "C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Australia2//Model_estimates//Cov_Aus_new.png")

######### Dispersion 
disp_names <- grep("^lv.c.+4]$", unique(mcmc.nb$Parameter), value=TRUE)

mcmc.disp<-mcmc.nb%>%filter(Parameter%in%disp_names)%>%
  mutate(Functional_group=rep(Functional_group,each=dim(mcmc.raw)[[1]]))

## Sumarize

quantt<-mcmc.disp%>%group_by(Functional_group)%>% 
  do(data.frame(t(quantile(.$value, probs = c(0.025, 0.975)))))

mean<-mcmc.disp%>%group_by(Functional_group)%>% 
  do(data.frame(t(mean(.$value))))

sum_disp<-inner_join(mean,quantt)
colnames(sum_disp)<-c("Functional_group","Mean","Lower","Upper")

sum_disp<-sum_disp%>%mutate(Sig=ifelse(Lower<0 & Upper<=0 | Lower>=0 & Upper>0,1,0))


sum_disp<-inner_join(sum_disp,labelset)
sum_disp$Group<-factor(sum_disp$Group,levels=c("Hard corals","Soft Coral","Other Invertebrates","Algae","Other"))

p_sum_disp<-ggplot()+geom_point(data=sum_disp,aes(y=Mean,x=fct_reorder(Functional_group,as.numeric(Group))),col="#1F78B4")+
  geom_errorbar(data=sum_disp,aes(x=Functional_group,ymin=Lower,ymax=Upper),col="#1F78B4",width=0.1)+
  coord_flip()+
  theme(axis.text.x = element_text(size=13),legend.position="none",
        panel.grid = element_blank(),plot.title = element_text(hjust=0,vjust=2,size=18, face="bold"),
        axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
        strip.text=element_text(size=14),strip.background = element_rect(fill="gray98"))+ 
  geom_hline(yintercept = 0,linetype="dashed")+xlab("Functional groups")+ylab("Estimates")

ggsave(plot = p_sum_disp,width=6, height=6, file = "C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Australia2//Model_estimates//Dispersion_Aus_new.png")


################## extract covariances and correlations due to co-occurence

res.cor<-get.residual.cor(model.list, est="mean")
label<-str_remove(colnames(res.cor$sig.cor),'_Count')
label<-str_replace(label,"[.]","-")

### Order for matrix 

ord.mat<-sum_disp%>%dplyr::select(Functional_group,Group)%>%mutate_if(is.factor,as.numeric)%>%arrange(Group)

p.mat.prep<-res.cor$sig.cor
diag(p.mat.prep)=0
p.mat.prep<-p.mat.prep/100
p.mat.prep<- abs(p.mat.prep)
p.mat.prep.f<-p.mat.prep%>%replace(.==0,1)

p.mat.prep.f<-1-p.mat.prep.f

p.mat.prep.f<-p.mat.prep.f[ord.mat$Functional_group,ord.mat$Functional_group]

cor.mat<-res.cor$correlation
diag(cor.mat)=0
cor.mat<-cor.mat[ord.mat$Functional_group,ord.mat$Functional_group]

library(corrplot)
library(RColorBrewer)
# 
# # 
pdf(width=6, height=6,file="C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Indian_ocean2//Network//Correlation_matrix_Aus_new.pdf")
corrplot(cor.mat, type="upper", tl.col="black",p.mat=p.mat.prep.f, sig.level = 0.01,pch=16,pch.cex=.8,pch.col="red",
         col=brewer.pal(n=11, name="PRGn"))
dev.off()


################## extract correlations due to enviro


env.cor<-get.enviro.cor(model.list, est="mean")

p2.mat.prep<-env.cor$sig.cor
diag(p2.mat.prep)=0
p2.mat.prep<-p2.mat.prep/100
p2.mat.prep<- abs(p2.mat.prep)
p2.mat.prep.f<-p2.mat.prep%>%replace(.==0,1)

p2.mat.prep.f<-1-p2.mat.prep.f
p2.mat.prep.f<-p2.mat.prep.f[ord.mat$Functional_group,ord.mat$Functional_group]



cor2.mat<-env.cor$cor
diag(cor2.mat)=0
cor2.mat<-cor2.mat[ord.mat$Functional_group,ord.mat$Functional_group]


pdf(width=6, height=6,file="C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Australia2//Network//Correlation_matrix_enviro_Aus_new.pdf")
corrplot(cor2.mat, type="upper", tl.col="black",p.mat=p2.mat.prep.f, sig.level = 0.01,pch=16,pch.cex=.8,pch.col="red",
         col=brewer.pal(n=11, name="PRGn"))
dev.off()

########## Model predictions

if(model.list$num.X>0 && model.list$lv.control$type=="independent"){
  pred<-predict(model.list,est="mean")%>%data.frame()
}
if(model.list$num.X>0 && model.list$lv.control$type=="exponential"){
  pred<-predict(model.list,est="mean",distmat=distmat.reef)%>%data.frame() 
}



pred.mean<-exp(pred[str_detect(names(pred), "linpred")])
pred.lower<-exp(pred[str_detect(names(pred), "lower")])
pred.upper<-exp(pred[str_detect(names(pred), "upper")])
colnames(pred.mean)<-label;colnames(pred.lower)<-label;colnames(pred.upper)<-label
pred.gathr<-pred.mean%>%gather(global_label,pred.mean)
pred.gathr1<-pred.lower%>%gather(global_label,pred.lower)
pred.gathr2<-pred.upper%>%gather(global_label,pred.upper)

df.full.pred<-cbind(df.full%>%arrange(global_label),pred.gathr$pred.mean,pred.gathr1$pred.lower,pred.gathr2$pred.upper)%>%
  dplyr::select(reef_name,global_label,sub_region,Count,25:27)
colnames(df.full.pred)<-c("reef_name","global_label","sub_region","Observed count","Predicted count","Lower","Upper")

df.full.pred<-df.full.pred%>%mutate(Sign = ifelse(Lower<`Observed count` & Upper>=`Observed count`,1,0))
df.full.pred3<-df.full.pred%>%gather(key="Mod",value="Count",c(`Observed count`,`Predicted count`))

p2<-ggplot(df.full.pred3,aes(x=global_label,y=log(Count+1),fill=Mod))+geom_boxplot()+theme_bw()+xlab("Functional groups")+
  theme(axis.text.y = element_text(size=11),legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(hjust=0.5,vjust=2,size=17, face="bold"),
        axis.text.x = element_text(size=11,angle=40,hjust=1),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15))+ 
  scale_fill_manual(values = c("white", "grey"))


#ggsave(plot = p2,width=8, height=6, file = "C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Australia2//Predictions//Pred.pdf")


p3<-ggplot(df.full.pred,aes(x=`Observed count`,y=`Predicted count`))+geom_point()+theme_bw()+xlab("Observed count")+
  ylab("Predicted count")+geom_abline(col="red",linetype="dashed")+theme_bw()+
  theme(axis.text.y = element_text(size=13),legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", size = 12),
        plot.title = element_text(hjust=0.5,vjust=2,size=17, face="bold"),
        axis.text.x = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15))

ggsave(plot = p3,width=6, height=6, file = "C://Users//vercello//Dropbox//GCI//Catlin//Research projects//Global//Reassembly_study//Results//Boral_outputs//Australia2//Predictions//Model_fit_Aus_new.png")

