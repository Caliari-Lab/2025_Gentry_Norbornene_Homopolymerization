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


setwd('C:/Users/jg9zk/Desktop/Grad/')
path <- '../Grad/Papers/norbornene homopolymerization/cell culture/same formulations/'

###################################################
# Read in cell shape metrics
###################################################
f <- 'cyto_metrics_no_filtering.csv'
df <- read.csv2(paste0(path,f),sep = ',')

#assign gels to each well plate row
df$group <- '4-arm PEG'
df$group[grep('B',df$well)] = 'NorHA8'
df$group[grep('C',df$well)] = 'NorHA18'
df$group[grep('D',df$well)] = 'NorHA40'
df$group <- factor(df$group,levels=c('4-arm PEG','NorHA8','NorHA18','NorHA40'))

df$day <- 'D3'
df$day[grep('_d7_',df$fov)] = 'D7'
df$day[grep('_d14_',df$fov)] = 'D14'
df$day <- factor(df$day,levels=c('D3','D7','D14'))

#pixels to microns
area_multiplier <- (682/1960)^2 #1960 pixels = 682 microns

#text to numbers
df$aspect <- as.numeric(df$aspect)
df$area <- as.numeric(df$area)*area_multiplier
df$perimeter <- as.numeric(df$perimeter)
df$form_factor <- as.numeric(df$form_factor)

df$Condition <- paste0(df$group,'_',df$day)
df$sample <- unname(sapply(df$fov,function(x) strsplit(x,'ROI')[[1]][1]))

###################################################
# Read in CAJAL embedding
###################################################
f <- 'embedding_df_gw.csv'
cajal_df <- read.csv2(paste0(path,f),sep = ',')
cajal_df$file_cell <- paste0(cajal_df$file_id,'_',cajal_df$mask)

###################################################
# Merge cell shape metrics and CAJAL data frames
###################################################
file_map <- unique(cajal_df[,c('file_id','file_name')])

#add file and cell IDs to cell shape metrics data frame
df$file_id <- file_map$file_id[match(df$file,file_map$file_name)]
df$file_cell <- paste0(df$file_id,'_',df$labels)

#remove unnecessary columns from each data frame
df <- df[,-which(colnames(df) %in% c('X','file2','fov','labels','file','file2','file_id'))]
cajal_df <- cajal_df[,c('Dim1','Dim2','leiden','medoid_cell','file_cell','color_code')]

#merge data frames and change data types
df <- cbind(df,cajal_df[match(df$file_cell,cajal_df$file_cell),c(1:4,6)])
df$Dim1 <- as.numeric(df$Dim1)
df$Dim2 <- as.numeric(df$Dim2)
df$polymer <- df$group
df$leiden <- factor(df$leiden, levels = sample(unique(df$leiden)))
df$color_code <- factor(df$color_code)

#remove cells that were removed by CAJAL or are too small
hist(log10(df$area[(!is.na(df$leiden)) & (df$area>(10^2.35))]),breaks='fd')
abline(v=2.3)
df <- df[!is.na(df$leiden) & (df$area>(10^2.35)),]

#cell count per condition
df %>% group_by(Condition) %>% 
  count()

###################################################
# Plot form factor
###################################################
pdf(paste0(path,'/form_factor_same_formulations_diff_polyms.pdf'), width=6, height=4)
out <- plot_gelwise_metrics(data=df,x='group',y='form_factor',color='group',
                            facet_col = 'day', violin_bounds = c(-10,10),
                            x_label='Polymer',y_label='Form Factor',cmap=col_pal)
#out$g+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1))
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>% 
  group_by(facet_col) %>%
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Polymer") %>% 
  #select(-x) %>% 
  mutate(xmin = group1, xmax = group2)
out$g +theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1))+
  scale_x_discrete(labels = c('4-arm PEG',expression('NorHA '[8]),
                              expression('NorHA '[18]),expression('NorHA '[40])))+
  stat_pvalue_manual(pwc, bracket.nudge.y = .3) 

summary(aov(y~x+facet_col,formfactor_stats))

formfactor_stats %>% 
  group_by(facet_col) %>%
  group_map(~TukeyHSD(aov(y ~ x, .x))$x) %>% 
  do.call(rbind,.)

dev.off()

###################################################
# Plot area
###################################################
pdf(paste0(path,'/projected_area_same_formulations_diff_polyms.pdf'), width=6, height=4)
out <- plot_gelwise_metrics(data=df,x='group',y='area',color='group',
                            facet_col = 'day',
                            x_label='Polymer',y_label=expression('Projected Area ('*mu*'m'^2*')'),cmap=col_pal)
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>% 
  group_by(facet_col) %>%
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Polymer") %>% 
  #select(-x) %>% 
  mutate(xmin = group1, xmax = group2)

out$g +theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1))+
  scale_x_discrete(labels = c(expression('Nor'[4]*'PEG'),expression('Nor'[8]*'HA'),
                              expression('Nor'[18]*'HA'),expression('Nor'[40]*'HA')))+
  stat_pvalue_manual(pwc, tip.length = .0075, bracket.nudge.y = 1500, step.increase = .015) 
dev.off()

###################################################
# Plot aspect ratio
###################################################
pdf(paste0(path,'/aspect_ratio_same_formulations_diff_polyms.pdf'), width=6, height=5)
out <- plot_gelwise_metrics(data=df,x='group',y='aspect',color='group',
                            facet_col = 'day', ylims = c(1,18), violin_bounds = c(1,20),
                            x_label='Polymer',y_label='Aspect Ratio',cmap=col_pal)
#out$g <- out$g + ylim(c(0,10))
formfactor_stats <- out$gel_wise_stats
pwc <- formfactor_stats %>%
  group_by(facet_col) %>%
  pairwise_t_test(y ~ x, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "Polymer") %>%
  #select(-x) %>%
  mutate(xmin = group1, xmax = group2)
out$g + theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1)) +
  scale_x_discrete(labels = c(expression('Nor'[4]*'PEG'),expression('Nor'[8]*'HA'),
                              expression('Nor'[18]*'HA'),expression('Nor'[40]*'HA')))+
  stat_pvalue_manual(pwc, tip.length = .0075, bracket.nudge.y = 15)

formfactor_stats %>% 
  group_by(facet_col) %>%
  group_map(~TukeyHSD(aov(y ~ x, .x))$x) %>% 
  do.call(rbind,.)

dev.off()

###################################################
# Assign leiden clusters to morphological category
###################################################
df$morphology <- 'Spherical'
df$morphology[df$leiden %in% c(26,1,4,6,23,8,14)] <- 'Star'
df$morphology[df$leiden %in% c(16,19,20,13)] <- 'Spindle-1'
df$morphology[df$leiden %in% c(9,7,18,22)] <- 'Spindle-2'
df$morphology <- factor(df$morphology,levels=c('Spindle-2','Spindle-1','Star','Spherical'))

df$sample <- paste0(df$well,'_',df$day)

#Calculate fraction of cells in each morphological category in each hydrogel group at each time point
cluster_dat <- df %>% 
  group_by(sample,day,polymer, .drop=FALSE) %>% 
  dplyr::count(morphology) %>% 
  mutate(n_cells = sum(n)) %>% 
  mutate(frac_n = n/n_cells)
cluster_dat$condition <- paste0(cluster_dat$morphology,'_',cluster_dat$polymer)
cluster_dat$morphology <- factor(cluster_dat$morphology,
                                 levels=rev(c('Spindle-2','Spindle-1','Star','Spherical')))

###################################################
# Plot cluster dynamics
###################################################
pdf(paste0(path,'/Fig5_morph_cluster_dynamics.pdf'), width=8, height=2.5)
g <- cluster_dat %>% filter(!is.na(frac_n)) %>% 
  ggplot(aes(x=day, y=frac_n, group=condition, color=polymer))+
  stat_summary(aes(y = frac_n, group=condition, fill=polymer),
               alpha = .15,
               color=NA,
               fun.ymin = function(x) mean(x) - sd(x), 
               fun.ymax = function(x) mean(x) + sd(x), 
               geom = "ribbon")+
  stat_summary(aes(y = frac_n,group=condition,color=polymer), 
               fun.y=mean, geom="line", linewidth=.8)+
  stat_summary(aes(y = frac_n,group=condition,color=polymer,shape=polymer), 
               fun.y=mean, geom="point",size=3)+
  facet_grid(cols=vars(morphology))+
  xlab('')+
  ylab(paste0('Proportion of cells'))+
  theme_classic()+
  theme(strip.background = element_blank(),
        legend.title=element_blank(),
        axis.text = element_text(color="black"),
        axis.ticks = element_line(color = "black"))+
  scale_color_manual(values = colorRampPalette(c('white',col_pal[3:4]))(length(unique(cluster_dat$polymer))+1)[2:5])+
  scale_fill_manual(values = colorRampPalette(c('white',col_pal[3:4]))(length(unique(cluster_dat$polymer))+1)[2:5])+
  scale_shape_manual(values = c(15,17,8,19))
g
dev.off()

###################################################
# Cluster dynamics statistics
###################################################
#pair-wise comparisons between hydrogel groups for each morph category at each time point
p_mat <- cluster_dat %>% 
  filter(!is.na(frac_n)) %>% 
  group_by(morphology,day) %>% 
  group_map(~ TukeyHSD(aov(frac_n ~ polymer, .x))[[1]]) %>% 
  do.call(rbind,.)
p_mat

#2-way ANOVA for each morph category
p_mat <- cluster_dat %>% 
  filter(!is.na(frac_n)) %>% 
  group_by(morphology) %>% 
  group_map(~ summary(aov(frac_n ~ day * polymer, .x))[[1]]) %>% 
  do.call(rbind,.)
p_mat

#3-way ANOVA
summary(aov(frac_n ~ day * morphology * polymer, cluster_dat))[[1]]


###################################################
# Retrieving medoids
###################################################
library(magick)
library(ggimage)

in_dir <- '../Grad/Papers/norbornene homopolymerization/cell culture/same formulations/cluster_medoids/'
medoid_df <- df[df$medoid_cell=='True',]
medoid_df$image <- sapply(medoid_df$leiden, function(x){
  image_path <- paste0(in_dir,'cluster_',x,'.png')
  img <- image_read(image_path)
  img <- image_transparent(img, "white")
  image_write(img, gsub('.png','_transparent.png',image_path), format = "png")
  return(gsub('.png','_transparent.png',image_path))
})
medoid_df$rel_size <- sapply(medoid_df$leiden, function(x){
  image_path <- paste0(in_dir,'cluster_',x,'.png')
  return(image_info(image_read(image_path))$width)
})
medoid_df$rel_size <- medoid_df$rel_size/max(medoid_df$rel_size)#ensures that masks are scaled correctly

###################################################
# Manually assign colors to all clusters
###################################################
df$color_code2 <- 'black'#initialize column
#Spindle-2
#"#F7BAA1" "#F07543" "#BD4210"
df$color_code2[df$leiden == 18] = "#D599C0"
df$color_code2[df$leiden == 7] = "#AB3381"
df$color_code2[df$leiden == 9] = "#77004E"
df$color_code2[df$leiden == 22] = "#D599C0"

#Spindle-1
#"#D599C0" "#AB3381" "#77004E" "#3B0027"
df$color_code2[df$leiden == 13] = "#F07543"
df$color_code2[df$leiden == 16] = "#F7BAA1"
df$color_code2[df$leiden == 19] = "#F07543"
df$color_code2[df$leiden == 20] = "#F7BAA1"

#Star
#"#A7D8E9" "#4FB2D4" "#1C7FA1" "#0E3F50"
df$color_code2[df$leiden == 26] = "#A7D8E9"
df$color_code2[df$leiden == 1] = "#1C7FA1" 
df$color_code2[df$leiden == 4] = "#A7D8E9"
df$color_code2[df$leiden == 8] = "#A7D8E9"
df$color_code2[df$leiden == 23] = "#1C7FA1"
df$color_code2[df$leiden == 6] = "#4FB2D4"
df$color_code2[df$leiden == 14] = "#4FB2D4"

#Spherical
#"#A4ADC2" "#4A5B85" "#172852" "#0B1429"
df$color_code2[df$leiden == 17] = "#A4ADC2"
df$color_code2[df$leiden == 3] = "#4A5B85"
df$color_code2[df$leiden == 15] = "#172852"
df$color_code2[df$leiden == 12] = "#A4ADC2"
df$color_code2[df$leiden == 28] = "#0B1429"
df$color_code2[df$leiden == 2] = "#4A5B85"
df$color_code2[df$leiden == 24] = "#172852"
df$color_code2[df$leiden == 0] = "#A4ADC2"
df$color_code2[df$leiden == 5] = "#4A5B85"
df$color_code2[df$leiden == 10] = "#A4ADC2"
df$color_code2[df$leiden == 27] = "#4A5B85"
df$color_code2[df$leiden == 11] = "#0B1429"
df$color_code2[df$leiden == 21] = "#A4ADC2"

###################################################
# Plotting the CAJAL embedding with medoids
###################################################
pdf(paste0(path,'/Fig5_umap_leiden_4color_morph.pdf'), width=8, height=6)
ggplot(df, aes(Dim1, Dim2)) +
  geom_point(color=df$color_code2,size=.75,alpha=1)+
  geom_image(data = medoid_df, aes(image = image, size = I(rel_size/4)), 
             by='width', asp=1,color='black') +
  xlim(-5,16.5)+
  umap_preset(df$color_code2)+
  theme(legend.position = "none")+
  guides(colour = guide_legend(override.aes = list(size=2)))
dev.off()

###################################################
# Contour plots showing dynamics
###################################################
pdf(paste0(path,'/Fig5_umap_contour_polymer_v_day_black.pdf'), width=3.5, height=3.5)
ggplot(df, aes(x=Dim1, y=Dim2, fill=day))+
  #geom_point(size=.1,alpha=.25)+
  umap_preset(df$day,cmap = c("#AB3381","#172852","#4FB2D4"))+
  facet_grid(cols=vars(day),rows=vars(polymer),scales = 'free')+
  geom_density2d(color='black',alpha=1,linewidth=.1,adjust=1)+
  scale_x_reverse(limits=c(17,-7))+
  ylim(-.5,16)+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=Inf)) +
  geom_hline(aes(yintercept=Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  #xlim(-7,17)+
  theme(legend.position = "none",panel.spacing = unit(-.5, "lines"),axis.line=element_line(color='black'))
dev.off()

###################################################
# Others
###################################################
pdf(paste0(path,'/FigS10_cajal_colored_by_polymer.pdf'), width=5, height=4)
ggplot(df, aes(x=Dim1, y=Dim2, color=polymer))+
  geom_point(size=.5)+
  umap_preset(df$polymer)+
  guides(colour = guide_legend(override.aes = list(size=2)))
dev.off()

pdf(paste0(path,'/FigS10_cajal_contour_polymer_by_day.pdf'), width=6, height=4)
ggplot(df, aes(x=Dim1, y=Dim2, color=polymer))+
  geom_point(aes(color=polymer),size=.01)+
  geom_density_2d(color='black',linewidth=.2, adjust=.33, )+
  facet_grid(cols=vars(polymer),rows=vars(day))+
  xlim(-4,16)+
  ylim(-2,20)+
  umap_preset(df$polymer)
dev.off()

ggplot(df, aes(x=morphology, y=form_factor, fill=morphology))+
  geom_violin()+
  theme_classic()

ggplot(df, aes(x=morphology, y=area, fill=morphology))+
  geom_violin()+
  theme_classic()

ggplot(df, aes(x=morphology, y=aspect, fill=morphology))+
  geom_violin()+
  theme_classic()
