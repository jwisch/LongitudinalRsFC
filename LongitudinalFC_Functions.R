

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))}

MatchbyNearestDate<-function(dataframe1, dataframe2, IDcolumnName, datedf1, datedf2){
  df.combined<-merge(dataframe1, dataframe2, by = IDcolumnName)
  unique.id<-unique(df.combined[,IDcolumnName])
  k<-1
  df.result<-df.combined[FALSE,]
  for(i in 1:length(unique(df.combined[,IDcolumnName]))){
    df.working<-subset(df.combined, df.combined[,IDcolumnName] == unique.id[i])
    scandates<-unique(df.working[,datedf1])
    for(j in 1:length(scandates)){
      df.working2<-subset(df.working, df.working[,datedf1] == scandates[j])
      ind.val<-which.min(abs(df.working2[,datedf1]-df.working2[,datedf2]))
      df.result[k,]<-df.working2[ind.val,]
      k<-k+1
    }
  }
  return(df.result)
}

var_coord_func <- function(loadings, comp.sdev){
  loadings*comp.sdev}

contrib <- function(var.cos2, comp.cos2){var.cos2/comp.cos2}

MakeBarPlots<-function(Component, ComponentTitle){
  theme_set(theme_bw())  
  pcaPlot<-as.data.frame(Component)
  pcaPlot<-cbind(rownames(pcaPlot), pcaPlot[,1])
  colnames(pcaPlot)<-c("region", "contribution")
  pcaPlot<-as.data.frame(pcaPlot)
  pcaPlot$contribution<-as.numeric(as.character(pcaPlot$contribution))
  pcaPlot <- pcaPlot[order(-pcaPlot$contribution), ]  # sort
  # Diverging Barcharts
  p<-ggplot(pcaPlot, aes(x=reorder(region, contribution), y = contribution, label=contribution)) + 
    geom_bar(stat='identity', width=.5)  +
    #coord_cartesian(ylim = c(-0.4, 0.4)) +
    scale_y_continuous(limits = c(0, 0.15))+
    labs(subtitle="Explains 25% of the data variance", 
         title= ComponentTitle) + 
    coord_flip() + ylab("Proportion of Contribution") + xlab("")
  theme(legend.position = "none")
  return(p)}


LikelihoodRatioTests_RandInt<-function(MEASURE, df){
  significancetests<-list()
  Model.Total.null<-lmer(MEASURE ~ cdr05+years+apoe4+cdr05:years+GENDER+EDUC+apoe4:years+(1|ID), df, REML = FALSE)
  summary(Model.Total.null)
  Model.Total<-lmer(MEASURE ~ cdr05+years+cdr05:years+GENDER+EDUC+(1|ID), df, REML = FALSE)
  significancetests[[1]]<-anova(Model.Total.null, Model.Total) #no APOE effect
  Model.Total<-lmer(MEASURE ~ cdr05+years+apoe4+cdr05:years+GENDER+apoe4:years+(1|ID), df, REML = FALSE)
  significancetests[[2]]<-anova(Model.Total.null, Model.Total) #no ed effect
  Model.Total<-lmer(MEASURE ~ cdr05+years+apoe4+cdr05:years+EDUC+apoe4:years+(1|ID), df, REML = FALSE)
  significancetests[[3]]<-anova(Model.Total.null, Model.Total) #significant gender effect p = 7.8e-03
  Model.Total<-lmer(MEASURE ~ years+apoe4+GENDER+EDUC+apoe4:years+(1|ID), df, REML = FALSE)
  significancetests[[4]]<-anova(Model.Total.null, Model.Total) #converter or not is not significant
  Model.Total<-lmer(MEASURE ~ cdr05+apoe4+GENDER+EDUC+(1|ID), df, REML = FALSE)
  significancetests[[5]]<-anova(Model.Total.null, Model.Total) #significant time effect p = 1.48e-14
  Model.Total<-lmer(MEASURE ~ cdr05+years+apoe4+GENDER+EDUC+apoe4:years+(1|ID), df, REML = FALSE)
  significancetests[[6]]<-anova(Model.Total.null, Model.Total) #no significant interaction cdr*years
  Model.Total<-lmer(MEASURE ~ cdr05+years+apoe4+cdr05:years+GENDER+EDUC+(1|ID), df, REML = FALSE)
  significancetests[[7]]<-anova(Model.Total.null, Model.Total) #no significant interaction apoe*years
  return(significancetests)}

#getting z scores on a region by region basis
normalize<-function(PARAM){
  (PARAM - mean(PARAM))/sd(PARAM)
}

set.seed(23)
MultipleCorrelations <- function(DF, ResponseVariable, formula) {
  model.correct<-lmer(formula, DF)
  DF[,"Corrected_Parameter"]<-fitted(model.correct)
  #DF[,"ResponseVariable"]<-ResponseVariable
  result<-rmcorr(DF[,"ID"], DF[,"FCSummary"], ResponseVariable, DF, CIs = "bootstrap",  nreps = 1000)
  return(result)
}


#partial correlation without multiple measures
pcorr_eqn <- function(x,y, variablestocontrolfor = NULL, digits = 4) {
  corr_coef <- round(pcor.test(x, y, variablestocontrolfor, method = "spearman")$estimate, digits = digits)
  #p_val<-round(spearman.test(x, y)$p.value, digits = digits)
  p_val<-formatC(pcor.test(x, y, variablestocontrolfor, method = "spearman")$p.value, digits = digits)
  paste("R =  ", corr_coef, "\n p =  ", p_val)
}

pcorr_eqn_All <- function(x,y, variablestocontrolfor = NULL, digits = 4) {
  corr_coef <- round(pcor.test(x, y, variablestocontrolfor, method = "spearman")$estimate, digits = digits)
  #p_val<-round(spearman.test(x, y)$p.value, digits = digits)
  p_val<-formatC(pcor.test(x, y, variablestocontrolfor, method = "spearman")$p.value, digits = digits)
  paste("R =  ", corr_coef, ", p =  ", p_val)
}

ATNplot<-function(DATAFRAME, X_param, Y_param, rmcorrOUTPUT_r, rmcorrOUTPUT_p, YLABEL, YPOS = 0.9){
  pformat<-ifelse(rmcorrOUTPUT_p < 0.001, formatC(rmcorrOUTPUT_p, format = "e", digits = 3), round(rmcorrOUTPUT_p, 4))
  ggplot(data = DATAFRAME, aes(x = X_param, y = Y_param, group = ID, colour = cdr05, shape = cdr05))+ 
    geom_line(alpha = 0.6)+ 
    geom_point(alpha = 0.6)+scale_color_manual("", values = c( "#FFB000", "#1437AD"))+
    scale_shape_manual("", values = c(1, 17))+
    xlab("FC Summary Value") + ylab(YLABEL)+
    geom_smooth(aes(x = Measure1, y = Measure2, group = 123), method = "lm", colour = "gray38")+
    scale_x_continuous(limits = c(-5, 5))+
    geom_text(aes(x = -5, y = 0.99*max(Y_param), 
                  label = paste("R = ", round(rmcorrOUTPUT_r, 4), sep = ""), vjust = "inward", hjust = "inward"), colour = "black")+
    geom_text(aes(x = -5, y = YPOS*max(Y_param), 
                  label = paste("p = ", pformat, sep = ""), vjust = "inward", hjust = "inward"), colour = "black")
}

plotSetUp<-function(DATAFRAME, PARTICIPANT, MEASURE1, MEASURE2){
  data.results<-as.data.frame(cbind(PARTICIPANT, MEASURE1, MEASURE2))
  colnames(data.results)<-c("ParticipantNo", "Measure1", "Measure2")
  df.plot<-as.data.frame(cbind(DATAFRAME, data.results))
  return(df.plot)}


CleanUpDates<-function(DF){
  if("first05" %in% names(DF)){DF[,"first05"]<-as.Date(DF[,"first05"], origin="1960-01-01")}
  if("cdr05" %in% names(DF)){DF[,"cdr05"]<-as.factor(DF[,"cdr05"])}
  if("ID" %in% names(DF)){DF[,"ID"]<-as.factor(DF[,"ID"])}
  if("BIRTH" %in% names(DF)){as.Date(format(as.Date(DF[,"BIRTH"], format = "%Y-%m-%d"), "19%y-%m-%d"), format = "%Y-%m-%d")}
  if("PET_Date" %in% names(DF)){DF[,"PET_Date"]<-as.Date(DF[,"PET_Date"], format = "%m/%d/%Y")
  if("years" %in% names(DF)){DF[,"years"]<-difftime(DF[,"PET_Date"], df[,"first05"], units = "weeks")
  DF[,"years"]<-as.numeric(DF[,"years"]/52)}}
  if("Session_Date" %in% names(DF)){DF[,"Session_Date"]<-as.Date(DF[,"Session_Date"], format="%Y-%m-%d")}
  return(DF)}

normalize<-function(PARAM){(PARAM - mean(PARAM))/sd(PARAM)}