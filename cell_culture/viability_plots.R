#################
# Load libraries
#################
library(dplyr)
library(ggbeeswarm)
library(scales)
library(ggplot2)
library(ggpubr)
library(rstatix)

col_pal <- c('#960162','#ED5315','#249FCA','#1D3367')
setwd("/Users/jg9zk/Desktop/Grad/Papers/")

##################
# Functions
##################
plot_gelwise_metrics <- function(data,x='',y='',color='',facet_col='',
                                 summarise1='sample',summarise2='Condition',
                                 x_label='',y_label='',cmap='Blues',ylims=c(0,NA),
                                 plot_violin=TRUE,violin_bounds = c(-Inf,Inf)){
  if(facet_col == ''){
    agg_df <- data %>% 
      select(one_of(c(x,y,color,summarise1,summarise2)))
  }else{
    agg_df <- data %>% 
      select(one_of(c(x,y,color,facet_col,summarise1,summarise2)))
  }
  
  
  if(facet_col == ''){
    agg_df <- agg_df %>% 
      group_by((!!sym(summarise1))) %>% 
      summarise(y=mean((!!sym(y))),x=unique((!!sym(x))),color=unique((!!sym(color))),
                summarise2=unique((!!sym(summarise2))))
    
    agg2_df <- agg_df %>% 
      group_by(summarise2) %>% 
      summarise(mean=mean(y),sd=sd(y),x=unique(x),
                color=unique(color))
  }else{
    agg_df <- agg_df %>% 
      group_by((!!sym(summarise1))) %>% 
      summarise(y=mean((!!sym(y))),x=unique((!!sym(x))),color=unique((!!sym(color))),
                facet_col=unique((!!sym(facet_col))),
                summarise2=unique((!!sym(summarise2))))
    
    agg2_df <- agg_df %>% 
      group_by(summarise2) %>% 
      summarise(mean=mean(y),sd=sd(y),x=unique(x),
                color=unique(color),facet_col=unique(facet_col))
    data$facet_col <- data[,facet_col]
  }
  
  g <- ggplot(data,aes(x=.data[[x]],y=.data[[y]],color=.data[[color]]))
  
  if(plot_violin == TRUE){
    g <- g + geom_violin(aes(fill=.data[[color]]),color=NA,alpha=.2, bounds = violin_bounds)
  }
  
  g <- g +
    geom_jitter(data=agg_df,aes(x=x,y=y,color=color),width=.3,size=1.25, height = 0)+
    geom_point(data=agg2_df,aes(x=x,y=mean,color=color),color='black',size=2.2)+
    geom_errorbar(data=agg2_df,inherit.aes = FALSE, 
                  aes(x=x, ymin = mean-sd, ymax = mean+sd),
                  color='black', width = 0.2, show.legend = FALSE)+
    ylab(y_label)+
    xlab(x_label)+
    ylim(ylims)+
    theme_classic()+
    theme(strip.background = element_blank(), 
          #text=element_text(size=24,  family="sans"), 
          #strip.text.x = element_blank(),
          axis.text = element_text(color="black"),
          axis.ticks = element_line(color = "black"))+
    scale_color_manual(values = colorRampPalette(cmap)(length(unique(data[,color]))))+
    scale_fill_manual(values = colorRampPalette(cmap)(length(unique(data[,color]))))
  
  if(facet_col != ''){
    g <- g + facet_grid(cols=vars(facet_col))
  }
  
  return(list(gel_wise_stats = agg_df, g = g))
}


umap_preset <- function(color_var,cmap = col_pal){
  list(theme_classic()+
         theme(strip.background = element_blank(),
               axis.title = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line=element_blank(),
               legend.title=element_blank()),
       
       scale_color_manual(values = colorRampPalette(cmap)(length(unique(color_var)))))
}

############################
# Viability same formulation
############################
path <- '../Papers/norbornene homopolymerization/cell culture/same formulations/'

#read in live cell shape metrics
f <- 'live_cells.csv'
df <- read.csv2(paste0(path,f),sep = ',')

df$group <- '4-arm PEG'
df$group[grep('B',df$well)] = 'NorHA8'
df$group[grep('C',df$well)] = 'NorHA18'
df$group[grep('D',df$well)] = 'NorHA40'
df$group <- factor(df$group,levels=c('4-arm PEG','NorHA8','NorHA18','NorHA40'))

live <- df
live$ld <- 'live'
live$area <- as.numeric(live$area)

#remove erroneous masks that are too small 
hist(log10(live$area),breaks='FD')
live <- live[live$area > 1000,]
hist(log10(live$area),breaks='FD')

#read in dead cell shape metrics
f <- 'dead_cells.csv'
df <- read.csv2(paste0(path,f),sep = ',')

df$group <- '4-arm PEG'
df$group[grep('B',df$well)] = 'NorHA8'
df$group[grep('C',df$well)] = 'NorHA18'
df$group[grep('D',df$well)] = 'NorHA40'
df$group <- factor(df$group,levels=c('4-arm PEG','NorHA8','NorHA18','NorHA40'))

dead <- df
dead$ld <- 'dead'
dead$area <- as.numeric(dead$area)

#remove erroneous masks that are too small 
hist(log10(dead$area),breaks='FD')
dead <- dead[dead$area > 10^2.1,]
dead <- dead[dead$area < 10^2.8,]
hist(log10(dead$area),breaks='FD')

#combine live and dead dataframes into one
viability <- rbind(live[,match(colnames(dead),colnames(live))],dead)

viability %>% 
  group_by(well,ld,group) %>% 
  count() %>% 
  group_by(well,group) %>% 
  summarize(v = n[ld=='live']/(n[ld=='live']+n[ld=='dead'])*100) %>% 
  group_by(group) %>% 
  summarize(mean(v),sd(v))

viability$viable <- ifelse(viability$ld=='live',100,0) #add binary column so plot_gelwise_metrics can plot mean viability
viability$sample <- viability$well
viability$Condition <- viability$group

#Plot mean viability
pdf(paste0(path,'/FigS10_d1_viability.pdf'), width=2, height=3)
out <- plot_gelwise_metrics(data=viability,
                            x='group',y='viable',color='group',
                            facet_col = '', violin_bounds = c(-10,10),plot_violin = FALSE,
                            x_label='Polymer',y_label='Viability (%)',cmap=col_pal)
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>% 
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni")%>%
  add_xy_position(x = "Polymer") %>% 
  mutate(xmin = group1, xmax = group2)
out$g +theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1),
             legend.position = "none")+
  scale_x_discrete(labels = c('4-arm PEG',expression('NorHA '[8]),
                              expression('NorHA '[18]),expression('NorHA '[40])))+
  stat_pvalue_manual(pwc[c(1,4,6),], bracket.nudge.y = .1, step.increase = .03)+
  ylim(0,120)

dev.off()


#Plot form factor of live cells
live$form_factor <- as.numeric(live$form_factor)
live$aspect <- as.numeric(live$aspect)
live$sample <- live$well
live$Condition <- live$group

pdf(paste0(path,'/FigS10_d1_viability_formfactor.pdf'), width=2, height=3)
out <- plot_gelwise_metrics(data=live[live$boundary == 'False',],
                            x='group',y='form_factor',color='group',
                            facet_col = '', violin_bounds = c(-10,10),
                            x_label='Polymer',y_label='Form Factor',cmap=col_pal)
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>% 
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni")%>%
  add_xy_position(x = "Polymer") %>% 
  mutate(xmin = group1, xmax = group2)
out$g +theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1),
             legend.position = "none")+
  scale_x_discrete(labels = c('4-arm PEG',expression('NorHA '[8]),
                              expression('NorHA '[18]),expression('NorHA '[40])))+
  stat_pvalue_manual(pwc[c(1,4,6),], bracket.nudge.y = .1, step.increase = .03)+
  ylim(0,1.2)
dev.off()

summary(aov(y~x,formfactor_stats))

############################
# Viability same mechanics
############################
path <- '../Papers/norbornene homopolymerization/cell culture/same mechanics/'

#read in live cell shape metrics
f <- 'live_cells.csv'
df <- read.csv2(paste0(path,f),sep = ',')

df$group <- '4-arm PEG'
df$group[grep('B',df$well)] = 'NorHA8'
df$group[grep('C',df$well)] = 'NorHA18'
df$group[grep('D',df$well)] = 'NorHA40'
df$group <- factor(df$group,levels=c('4-arm PEG','NorHA8','NorHA18','NorHA40'))

live <- df
live$ld <- 'live'
live$area <- as.numeric(live$area)

#remove erroneous masks that are too small 
hist(log10(live$area),breaks='FD')
live <- live[live$area > 1000,]
hist(log10(live$area),breaks='FD')

#read in dead cell shape metrics
f <- 'dead_cells.csv'
df <- read.csv2(paste0(path,f),sep = ',')

df$group <- '4-arm PEG'
df$group[grep('B',df$well)] = 'NorHA8'
df$group[grep('C',df$well)] = 'NorHA18'
df$group[grep('D',df$well)] = 'NorHA40'
df$group <- factor(df$group,levels=c('4-arm PEG','NorHA8','NorHA18','NorHA40'))

dead <- df
dead$ld <- 'dead'
dead$area <- as.numeric(dead$area)

#remove erroneous masks that are too small 
hist(log10(dead$area),breaks='FD')
dead <- dead[dead$area > 10^2.1,]
dead <- dead[dead$area < 10^2.8,]
hist(log10(dead$area),breaks='FD')

#combine live and dead dataframes into one
viability <- rbind(live[,match(colnames(dead),colnames(live))],dead)

viability %>% 
  group_by(well,ld,group) %>% 
  count() %>% 
  group_by(well,group) %>% 
  summarize(v = n[ld=='live']/(n[ld=='live']+n[ld=='dead'])*100) %>% 
  group_by(group) %>% 
  summarize(mean(v),sd(v))

viability$viable <- ifelse(viability$ld=='live',100,0)#add binary column so plot_gelwise_metrics can plot mean viability
viability$sample <- viability$well
viability$Condition <- viability$group

#Plot mean viability
pdf(paste0(path,'/FigS11_d1_viability.pdf'), width=2, height=3)
out <- plot_gelwise_metrics(data=viability,
                            x='group',y='viable',color='group',
                            facet_col = '', violin_bounds = c(-10,10),plot_violin = FALSE,
                            x_label='Polymer',y_label='Viability (%)',cmap=col_pal)
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>% 
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni")%>%
  add_xy_position(x = "Polymer") %>% 
  mutate(xmin = group1, xmax = group2)
out$g +theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1),
             legend.position = "none")+
  scale_x_discrete(labels = c('4-arm PEG',expression('NorHA '[8]),
                              expression('NorHA '[18]),expression('NorHA '[40])))+
  stat_pvalue_manual(pwc[c(1,4,6),], bracket.nudge.y = .1, step.increase = .03)+
  ylim(0,120)

dev.off()

summary(aov(y~x,formfactor_stats))


#Plot form factor of live cells
live$form_factor <- as.numeric(live$form_factor)
live$aspect <- as.numeric(live$aspect)
live$sample <- live$well
live$Condition <- live$group

pdf(paste0(path,'/FigS11_d1_viability_formfactor.pdf'), width=2, height=3)
out <- plot_gelwise_metrics(data=live[live$boundary == 'False',],
                            x='group',y='form_factor',color='group',
                            facet_col = '', violin_bounds = c(-10,10),
                            x_label='Polymer',y_label='Form Factor',cmap=col_pal)
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>% 
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni")%>%
  add_xy_position(x = "Polymer") %>% 
  mutate(xmin = group1, xmax = group2)
pwc$y.position[6] = pwc$y.position[1]
out$g +theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1),
             legend.position = "none")+
  scale_x_discrete(labels = c('4-arm PEG',expression('NorHA '[8]),
                              expression('NorHA '[18]),expression('NorHA '[40])))+
  stat_pvalue_manual(pwc[c(1,4,6),], bracket.nudge.y = .1, step.increase = .03)+
  ylim(0,1.2)
dev.off()

summary(aov(y~x,formfactor_stats))