library(tidyr)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(survival)
library(rms)
library(survminer)
library(caret)
library(tableone)

### load SEER data from different registries (see SEER-Dictionary for details)
c9 <- fread("crc9_cust_2016.csv") %>% 
  data.frame(.,stringsAsFactors = F,check.names = F)
c5 <- fread("crc5_cust_2016.csv")  %>% 
  data.frame(.,stringsAsFactors = F,check.names = F)
c4 <- fread("crc4_cust_2016.csv")  %>% 
  data.frame(.,stringsAsFactors = F,check.names = F)
cl <- fread("crclo_cust_2016.csv")  %>% 
  data.frame(.,stringsAsFactors = F,check.names = F)

aa <-c("MAR_STAT","PUBCSNUM","REG","ST_CNTY","AGE_DX","SEX",
       "RAC_RECA","RAC_RECY","ORIGRECB","YEAR_DX","MDXRECMP","REPT_SRC","DX_CONF","SEQ_NUM","FIRSTPRM",
       "PRIMSITE","HISTO3V","GRADE","EOD10_PN","EOD10_NE","EOD10_ND","EOD10_SZ",
       "EOD10_EX","CSTUMSIZ","DAJCCT","DAJCCN","DAJCCM","DAJCCSTG",
       #"CSMETSDXB_PUB","CSMETSDXBR_PUB","CSMETSDXLIV_PUB","CSMETSDXLUNG_PUB",
       "SURGPRIF","NO_SURG","SURGSITF","SURGSITE","SS_SURG",
       "SRV_TIME_MON","SRV_TIME_MON_FLAG","VSRTSADX","ODTHCLASS",
       "RADIATNR", "RAD_BRNR", "RAD_SURG","CHEMO_RX_REC",
       "TUMSIZS","DASRCT","DASRCN","DASRCM"
)#"RADIATNR", "RAD_BRNR", "RAD_SURG","CHEMO_RX_REC"
table(aa %in%names(c9))
table(aa %in%names(c5))
aa[!aa%in%names(c9)]

crca <- rbind(c9[,aa],c5[,aa],c4[,aa],cl[,aa])
#crca <- rbind(c9[,aa],c5[,aa],c4[,aa])

### general patient selection (see SEER-Dictionary for details)
bc <- crca %>% 
  filter(YEAR_DX>=1988,
         YEAR_DX<=2015,
         AGE_DX>=18 & AGE_DX!=999,
         SEQ_NUM<=1,
         DX_CONF%in%c(1,2,4),
         HISTO3V%in%c(8140,8144,8210,8211,8221,8220,8261,8262,8263,8260,8480,8481,8490),
         SRV_TIME_MON<9999,
         PRIMSITE!="C260",
         (DAJCCM==10|(EOD10_EX==85)|(EOD10_ND==7)), #crca <- rbind(c9[,aa],c5[,aa],c4[,aa],cl[,aa])
         !SURGPRIF%in%c(99),
         (SURGSITF%in%c(0,3:5)|(SURGSITE%in%c(0,3:8)))
         
  )

### tumor location regrouping
bc$loc3 <- ifelse(bc$PRIMSITE%in%c("C180","C181","C182","C183","C184"),1,
                  ifelse(bc$PRIMSITE%in%c("C185","C186","C187"),2,
                         ifelse(bc$PRIMSITE%in%c("C199"),3,
                                ifelse(bc$PRIMSITE%in%c("C188","C189"),9,4))))

table(bc$loc3 )


### HSA data
HSA <- fread("HSA.csv")%>% 
  data.frame(.,stringsAsFactors = F,check.names = F)

bc_prepared <- inner_join(bc,HSA,by="ST_CNTY") ### merge in one
NA_num <- apply(bc_prepared,2,function(x){length(which(is.na(x)))})
sort(NA_num,decreasing = T)

## define covariates
CRC_prepared <- bc_prepared %>% 
  #dplyr::select(-c(EOD10_ND,EOD10_SZ,EOD10_EX,SURGSITE,SS_SURG)) %>% 
  mutate(region1=ifelse(REG<=1527,9,ifelse(REG>=1541,5,4)),
         region2=ifelse(REG%in%c(1501,1523,1525,1526,1529,1531,1535,1541),"west",
                        ifelse(REG%in%c(1520,1522,1542),"midwest",
                               ifelse(REG%in%c(1502,1544),"northeast","south"))),
         race=ifelse(RAC_RECA==1&ORIGRECB==0,1,
                     ifelse(RAC_RECA==2&ORIGRECB==0,2,
                            ifelse(ORIGRECB==1,3,4))),
         marry=ifelse(MAR_STAT==2,1,
                      ifelse(MAR_STAT==5,0,9)),
         grade2=ifelse(GRADE<=2,1,
                       ifelse(GRADE==9,9,2)),
         age6=I(AGE_DX>=50)+I(AGE_DX>=60)+I(AGE_DX>=70)+I(AGE_DX>=80),
         exitus = ifelse(VSRTSADX==1|ODTHCLASS==1,1,0),
         
  ) 


CRC_prepared$SURGSITE[is.na(CRC_prepared$SURGSITE)] <- 9999
CRC_prepared$SURGSITF[is.na(CRC_prepared$SURGSITF)] <- 9999
CRC_prepared$CSTUMSIZ[is.na(CRC_prepared$CSTUMSIZ)] <- 9999
CRC_prepared$EOD10_SZ[is.na(CRC_prepared$EOD10_SZ)] <- 9999
#CRC_prepared$TUMSIZS[is.na(CRC_prepared$TUMSIZS)] <- 9999

## PMTR ###
CC_df <- CRC_prepared %>% 
  #filter(race==1) %>% 
  #filter(loc3%in%c(2:4)) %>%
  filter(loc3!=9) %>% 
  filter(YEAR_DX%in%c(2005:2015)) %>% 
  filter(SURGPRIF%in%c(0,30:80)) %>% 
  mutate(age=AGE_DX,
         yr=YEAR_DX,
         time=SRV_TIME_MON,
         size4=ifelse(CSTUMSIZ<20|EOD10_SZ<20,1, # |TUMSIZS<20
                      ifelse(CSTUMSIZ<40|EOD10_SZ<40,2,
                             ifelse(CSTUMSIZ<60|EOD10_SZ<60,3,
                                    ifelse(CSTUMSIZ<990|EOD10_SZ<990,4,9)))),
         surgMet=ifelse(SURGSITF==0|SURGSITE==0,0,1),
         surgPrim=ifelse(SURGPRIF>0,1,0)
  ) %>% 
  mutate(surgPM=ifelse(surgPrim==1&surgMet==1,1,ifelse(is.na(surgPrim)|is.na(surgMet),0,0))) %>% 
  mutate(surgPM3=ifelse(surgPrim==1&surgMet==0,2,ifelse(surgPrim==0&surgMet==1,2,surgPM))) %>% 
  mutate(surgPM4=ifelse(surgPrim==1&surgMet==0,2,ifelse(surgPrim==0&surgMet==1,3,surgPM))) %>% 
  #.[,c(2:6,25,26,42,43,46:ncol(.))]
  .[,c(2:6,25,26,46,47,50:ncol(.))]

table(CC_df$surgPM);table(CC_df$surgPM4);table(CC_df$size4)

## baseline -- HSA PM
HSAnum <- data.frame(table(CC_df$HSA)) 
HSA_list <- sort(unique(CC_df$HSA))
for(ss in 1:length(HSA_list)){
  if(HSAnum$Freq[HSAnum$Var1==HSA_list[ss]]<50){
    CC_df$HSA[CC_df$HSA==HSA_list[ss]] <- HSA_list[ss+1]
  }
  if(ss==length(HSA_list) & HSAnum$Freq[HSAnum$Var1==HSA_list[ss]]<50){
    CC_df$HSA[CC_df$HSA==HSA_list[ss]] <- HSA_list[ss-1]
  }
}
sort(table(CC_df$HSA))

hsa_ls <- lapply(split(CC_df,CC_df$HSA),function(x){
  return(data.frame(HSA=unique(x$HSA),pc=length(which(x$surgPM3==1))/nrow(x)))
})
hsa_df <- data.frame(do.call(rbind,hsa_ls))
tt <- quantile(hsa_df$pc,probs = c(1/3,2/3))
qq <- quantile(hsa_df$pc,probs = seq(.25,.75,.25))
# group by median, tertiles, quartiles
CC_df <- inner_join(CC_df,hsa_df,by="HSA") %>% 
  mutate(pc2=ifelse(pc<median(pc),1,2)) %>% 
  mutate(pc3=ifelse(pc<tt[1],1,ifelse(pc<tt[2],2,3))) %>% 
  mutate(pc4=ifelse(pc<qq[1],1,ifelse(pc<qq[2],2,ifelse(pc<qq[3],3,4)))) 

## baseline -- HSA Prim
hsa_ls2 <- lapply(split(CC_df,CC_df$HSA),function(x){
  return(data.frame(HSA=unique(x$HSA),ppc=length(which(x$surgPM3==2))/nrow(x)))
})
hsa_df2 <- data.frame(do.call(rbind,hsa_ls2))
tt2 <- quantile(hsa_df2$ppc,probs = c(1/3,2/3))
qq2 <- quantile(hsa_df2$ppc,probs = seq(.25,.75,.25))

CC_df <- inner_join(CC_df,hsa_df2,by="HSA") %>% 
  mutate(ppc2=ifelse(ppc<median(ppc),1,2)) %>% 
  mutate(ppc3=ifelse(ppc<tt2[1],1,ifelse(ppc<tt2[2],2,3))) %>% 
  mutate(ppc4=ifelse(ppc<qq2[1],1,ifelse(ppc<qq2[2],2,ifelse(ppc<qq2[3],3,4)))) 

## factorial design
CC_df$pc22 <- factor(CC_df$pc2 + CC_df$ppc2*10,
                     labels = c("LowPM_LowPrim","HighPM_LowPrim","LowPM_HighPrim","HighPM_HighPrim"))
table(CC_df$pc22)

myVars <- names(CC_df)[c(5:8,10:13,16,18,23,ncol(CC_df))]
dtb_factor <- data.frame(lapply(CC_df[myVars],factor));str(dtb_factor)

tb_group<-CreateTableOne(vars = myVars, strata = c("pc22"), data = dtb_factor, # here choose tertile for grouping 
                         factorVars = myVars, includeNA = TRUE) 
print(tb_group)
#write.csv(print(tb_group),"results_demo/Table_by_HSA.csv")

### IV analysis
library(purrr)
CC_use <- CC_df %>%
  filter(loc3%in%c(2:4)) %>%  ## left, 2:4; right, 1
  mutate(race2=ifelse(race==1,1,2),
         race3=ifelse(race==1,1,ifelse(race==2,2,3)),
         loc2=ifelse(loc3==1,1,2),
         year4=I(yr>2007)+I(yr>2010))# %>% filter(surgPM3<=1)

myVars <- names(CC_use)[c(32,23,5,13,10,11,36,8,18,12,6,7)] # 5:8,10:13,16,18,23,32
dtb_factor <- data.frame(lapply(CC_use[myVars],factor));str(dtb_factor)
tb_group<-CreateTableOne(vars = myVars[-1], strata = c("pc22"), data = dtb_factor, # here choose tertile for grouping 
                         factorVars = myVars, includeNA = TRUE) 
print(tb_group)
#write.csv(print(tb_group),"results_demo3/Table_by_HSA_left.csv")

myVars <- names(CC_use)[c(21,5,13,10,11,36,8,18,12,6,7)] # 5:8,10:13,16,18,23,32
dtb_factor <- data.frame(lapply(CC_use[myVars],factor));str(dtb_factor)
tb_group<-CreateTableOne(vars = myVars[-1], strata = c("surgPM"), data = dtb_factor, # here choose tertile for grouping 
                         factorVars = myVars, includeNA = TRUE) 
print(tb_group)

### 2SRI
xls <- map(split(CC_use,CC_use$ppc2+CC_use$race3*10),function(x){
  ft_IV <- glm(surgPM ~ pc2,family = "binomial",data = x) # here use HSA as a continuous variable as IV 
  x$res <- x$surgPM-predict(ft_IV,x,type = "response")
  return(x)
})
CC_use <- data.frame(do.call(rbind,xls))

# F statistic
f1 <- lm(surgPM ~ pc2,data = CC_use)
f2 <- lm(surgPM ~ 1,data = CC_use)
anova(f1,f2) 

#KM plot
hr=exp(coef(f3)[1])
ll=exp(confint(f3)[1,1])
ul=exp(confint(f3)[1,2])
pv <- ifelse(sf3$coefficients[1,5]<.001,"P<0.001",paste0("P=",sprintf("%.3f",sf3$coefficients[1,5])))

f4 <- cph(Surv(time,exitus)~strat(surgPM)+res, data=CC_use,x=TRUE,y=TRUE,surv=TRUE)
tCut <- 120
p1 <- survminer::ggsurvplot(fkm, data = CC_use, risk.table = TRUE, censor=FALSE,
                            palette=ggsci::pal_nejm("default")(7)[2:1],size=1,
                            legend=c(.8,0.9), legend.title='',legend.labs=c('Other', 'Prim/Met resection'),
                            title="Stage IV CRC, left",
                            subtitle=paste0("(n=",nrow(CC_use),")"),
                            font.title=22, font.subtitle=20,
                            ylab='Survival Probability', xlab='Time (months)',
                            tables.y.text=TRUE,
                            tables.y.text.col=TRUE, risk.table.title='No. at Risk', 
                            break.time.by=12,xlim=c(0,tCut),
                            font.x=22, font.y=22, font.tickslab=16, font.legend=16, 
                            font.caption=20, risk.table.fontsize=5,
                            tables.theme = survminer::theme_survminer(font.main = 22, font.y=22,
                                                                      font.x=22, font.tickslab=16))#; p1

p2 <- p1
p2$plot <- p2$plot+
  theme(plot.margin = margin(1,1,0,1,"cm"))
p2$table <- p2$table+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0,1,1,1,"cm"))

## PTR ####

CC_df <- CRC_prepared %>% 
  #filter(race==1) %>% 
  #filter(loc3%in%c(1)) %>%
  filter(loc3!=9) %>% 
  filter(YEAR_DX%in%c(2005:2015)) %>% 
  filter(SURGPRIF%in%c(0,30:80)) %>% 
  mutate(age=AGE_DX,
         yr=YEAR_DX,
         time=SRV_TIME_MON,
         size4=ifelse(CSTUMSIZ<20|EOD10_SZ<20,1, # |TUMSIZS<20
                      ifelse(CSTUMSIZ<40|EOD10_SZ<40,2,
                             ifelse(CSTUMSIZ<60|EOD10_SZ<60,3,
                                    ifelse(CSTUMSIZ<990|EOD10_SZ<990,4,9)))),
         surgMet=ifelse(SURGSITF==0|SURGSITE==0,0,1),
         surgPrim=ifelse(SURGPRIF>0,1,0)
  ) %>% 
  mutate(surgPM=ifelse(surgPrim==1&surgMet==1,1,ifelse(is.na(surgPrim)|is.na(surgMet),0,0))) %>% 
  mutate(surgPM3=ifelse(surgPrim==1&surgMet==0,2,ifelse(surgPrim==0&surgMet==1,2,surgPM))) %>% 
  mutate(surgPM4=ifelse(surgPrim==1&surgMet==0,2,ifelse(surgPrim==0&surgMet==1,3,surgPM))) %>% 
  #.[,c(2:6,25,26,42,43,46:ncol(.))]
  .[,c(2:6,25,26,46,47,50:ncol(.))]

table(CC_df$surgPM);table(CC_df$surgPM4);table(CC_df$size4)

## baseline -- HSA PM
hsa_ls <- lapply(split(CC_df,CC_df$HSA),function(x){
  return(data.frame(HSA=unique(x$HSA),pc=length(which(x$surgPM3==1))/nrow(x)))
})
hsa_df <- data.frame(do.call(rbind,hsa_ls))
tt <- quantile(hsa_df$pc,probs = c(1/3,2/3))
qq <- quantile(hsa_df$pc,probs = seq(.25,.75,.25))
# group by median, tertiles, quartiles
CC_df <- inner_join(CC_df,hsa_df,by="HSA") %>% 
  mutate(pc2=ifelse(pc<median(pc),1,2)) %>% 
  mutate(pc3=ifelse(pc<tt[1],1,ifelse(pc<tt[2],2,3))) %>% 
  mutate(pc4=ifelse(pc<qq[1],1,ifelse(pc<qq[2],2,ifelse(pc<qq[3],3,4)))) 

## baseline -- HSA Prim
hsa_ls2 <- lapply(split(CC_df,CC_df$HSA),function(x){
  return(data.frame(HSA=unique(x$HSA),ppc=length(which(x$surgPM3==2))/nrow(x)))
})
hsa_df2 <- data.frame(do.call(rbind,hsa_ls2))
tt2 <- quantile(hsa_df2$ppc,probs = c(1/3,2/3))
qq2 <- quantile(hsa_df2$ppc,probs = seq(.25,.75,.25))

CC_df <- inner_join(CC_df,hsa_df2,by="HSA") %>% 
  mutate(ppc2=ifelse(ppc<median(ppc),1,2)) %>% 
  mutate(ppc3=ifelse(ppc<tt2[1],1,ifelse(ppc<tt2[2],2,3))) %>% 
  mutate(ppc4=ifelse(ppc<qq2[1],1,ifelse(ppc<qq2[2],2,ifelse(ppc<qq2[3],3,4)))) 

## factorial design
CC_df$pc22 <- factor(CC_df$pc2 + CC_df$ppc2*10,
                     labels = c("LowPM_LowPrim","HighPM_LowPrim","LowPM_HighPrim","HighPM_HighPrim"))
table(CC_df$pc22)

myVars <- names(CC_df)[c(5:8,10:13,16,18,23,ncol(CC_df))]
dtb_factor <- data.frame(lapply(CC_df[myVars],factor));str(dtb_factor)

tb_group<-CreateTableOne(vars = myVars, strata = c("pc22"), data = dtb_factor, # here choose tertile for grouping 
                         factorVars = myVars, includeNA = TRUE) 
print(tb_group)


### IV analysis
library(purrr)
CC_use <- CC_df %>% 
  filter(surgPM4%in%c(0,2)) %>% 
  filter(loc3%in%c(1)) %>% 
  #filter(race!=1) %>% 
  mutate(race2=ifelse(race==1,1,2),
         race3=ifelse(race==1,1,ifelse(race==2,2,3)),
         loc2=ifelse(loc3==1,1,ifelse(loc3%in%c(2:4),2,NA)),
         year4=I(yr>2007)+I(yr>2010)) #%>% 

myVars <- names(CC_use)[c(32,5,13,10,11,36,8,18,12,6,7)] # 5:8,10:13,16,18,23,32
dtb_factor <- data.frame(lapply(CC_use[myVars],factor));str(dtb_factor)
tb_group<-CreateTableOne(vars = myVars[-1], strata = c("pc22"), data = dtb_factor, # here choose tertile for grouping 
                         factorVars = myVars, includeNA = TRUE) 
print(tb_group)
#write.csv(print(tb_group),"results_demo3/Table_by_PTR_HSA_left.csv")
### 2SRI
xls <- map(split(CC_use,CC_use$pc2+CC_use$race3),function(x){
  ft_IV <- glm(surgPrim ~ ppc2,family = "binomial",data = x) # here use HSA as a continuous variable as IV 
  x$res <- x$surgPrim-predict(ft_IV,x, type = "response")
  return(x)
})
CC_use <- data.frame(do.call(rbind,xls))
f1 <- lm(surgPrim ~ ppc2,data = CC_use)
f2 <- lm(surgPrim ~ 1,data = CC_use)
anova(f1,f2) # F statistic

# KM plot
hr=exp(coef(f3)[1])
ll=exp(confint(f3)[1,1])
ul=exp(confint(f3)[1,2])
pv <- ifelse(sf3$coefficients[1,5]<.001,"P<0.001",paste0("P=",sprintf("%.3f",sf3$coefficients[1,5])))

f4 <- cph(Surv(time,exitus)~strat(surgPrim)+res, data=CC_use,x=TRUE,y=TRUE,surv=TRUE)
tCut <- 120
p1 <- survminer::ggsurvplot(fkm, data = CC_use, risk.table = TRUE, censor=FALSE,
                            palette=ggsci::pal_nejm("default")(7)[2:1],size=1,
                            legend=c(.8,0.9), legend.title='',legend.labs=c('No surgery', 'Prim resection'),
                            title="Stage IV CRC, right",
                            subtitle=paste0("(n=",nrow(CC_use),")"),
                            font.title=22, font.subtitle=20,
                            ylab='Survival Probability', xlab='Time (months)',
                            tables.y.text=TRUE,
                            tables.y.text.col=TRUE, risk.table.title='No. at Risk', 
                            break.time.by=12,xlim=c(0,tCut),
                            font.x=22, font.y=22, font.tickslab=16, font.legend=16, 
                            font.caption=20, risk.table.fontsize=5,
                            tables.theme = survminer::theme_survminer(font.main = 22, font.y=22,
                                                                      font.x=22, font.tickslab=16))#; p1

p2 <- p1
p2$plot <- p2$plot+
  theme(plot.margin = margin(1,1,0,1,"cm"))
p2$table <- p2$table+
  theme(axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0,1,1,1,"cm"))

