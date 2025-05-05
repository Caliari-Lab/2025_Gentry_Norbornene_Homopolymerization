library(APRheoPlotR)

col_pal <- c('#960162','#ED5315','#249FCA','#1D3367')
path <- '../../Desktop/Grad/Papers/norbornene homopolymerization/rheology/'

##################################################################
# Figure 1. Spontaneous vs Photocure
##################################################################
f <- 'cdha_adha_vs_uv_gh.csv'

#read in and format data
dat <- readData(paste0(path,f), meta_from_file_name = TRUE)
df <- dat@long #retrieves data in long table format
meta <- dat@metadata #retrieves meta data about samples
meta$condition <- meta$V1
meta$condition <- ifelse(meta$condition=='0wt%','Spontaneous','Photocurable') #samples without NorHA (0wt%) are Spontaneous

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample']) #append meta data to df
df$measurement[df$measurement == 'UV Cure '] = 'GH time' # change measurement label for plotting purposes
df$condition <- factor(df$condition, levels = c('Spontaneous','Photocurable'))

pdf(file = paste0(path,'Fig1_uv_vs_immediate_storage.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'GH time',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Storage_Modulus_Pa',
                    color_variable = 'condition',
                    cmap = col_pal[c(1,3)],
                    dot_size = 2,
                    y_scale = 'linear')
out$g+theme(axis.text.x = element_text(angle = 90))
dev.off()

#t test of G'
ggplot_build(out$g)$data[[1]][,c('y','group')] %>%
  group_map(~ t.test(y ~ group, .x, paired = FALSE))


pdf(file = paste0(path,'Fig1_uv_vs_immediate_tanDelta.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'GH time',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Loss_Factor_1',
                    color_variable = 'condition',
                    cmap = col_pal[c(1,3)],
                    dot_size = 2,
                    y_scale = 'linear')
out$g+theme(axis.text.x = element_text(angle = 90))
dev.off()

pdf(file = paste0(path,'Fig1_uv_vs_immediate_stress_relax.pdf'), width = 4, height = 3)
out$g <- plot_stress_relax(df, selected_measurement = 'Stress Relaxation',
                       color_variable = 'condition',
                       cmap = col_pal[c(1,3)],
                       normalize = TRUE,
                       plot_mean=TRUE, plot_sd = TRUE,
                       x_scale='linear',y_scale='linear')
out$g
dev.off()


df %>%
  filter(measurement=='GH time') %>%
  group_by(sample) %>%
  do(tail(.,n=5)) %>%
  summarize(mean = mean(Loss_Factor_1),cond=unique(condition)) %>%
  group_map(~ t.test(mean ~ cond, .x, paired = FALSE))


##################################################################
# Figure 2. UV Cure NorHA + LAP
##################################################################
f <- '2wt%_23mod_norha_1mm_LAP_r1.csv'
dat <- readData(paste0(path,f),meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta$condition <- meta$V4

meta$condition[meta$condition=='1mm'] = '1 mM'
meta$condition[meta$condition=='2mm'] = '2 mM'
meta$condition[meta$condition=='4mm'] = '4 mM'
meta$condition[meta$condition=='8mm'] = '8 mM'

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'Fig2a_norha_lap_concs.pdf'), width = 4.5, height = 3)#, units = 'in')#, res = 300)
out$g <- plot_time_sweep(df, selected_measurement = 'UV Cure ',
                     plot_mean=TRUE,
                     color_variable = 'condition',
                     cmap = col_pal,
                     y_scale='log10', ylim = c(.01,NA),
                     xlim=c(5,NA),
                     legend_title = 'LAP',
                     highlight_intervals = list(c(10,10+300)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

out <- plot_plateau(df, selected_measurement = 'UV Cure ',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Storage_Modulus_Pa',
                    color_variable = 'condition',
                    cmap = col_pal,
                    y_scale = 'log10', ylim = c(10,1600))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))


##################################################################
# Figure S10. UV Cure GH excess thiol
##################################################################
f <- '2wt%_norha_xwt%_aldha_1_0_ad-to-cd_3_0adpep_0_0peg_r1.csv'
dat <- readData(paste0(path,f),meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta <- meta[meta$V8=='0.0peg',]
df <- df[df$sample %in% meta$sample,]

meta$V7 <- gsub('adpep','x',meta$V7)
meta$condition <- meta$V7
meta$condition <- factor(meta$condition, levels = c('1.0x','2.0x','3.0x'))

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'FigS5_excess_ad_GH_cure.pdf'), width = 8, height = 2.5)
out$g <- plot_time_sweep(df, selected_measurement = 'UV Cure ',
                     facet_col  = 'condition',
                     color_variable = 'condition',
                     cmap = col_pal[c(1,3)],
                     plot_mean=TRUE,
                     ylim = c(1,1500), y_scale='log10',
                     highlight_intervals = list(c(10,10+120)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'FigS5_excess_ad_GH_storage.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'UV Cure ',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Storage_Modulus_Pa',
                    color_variable = 'condition',
                    cmap = col_pal[c(1,3)],
                    dot_size = 2,
                    y_scale = 'linear', ylim = c(0,NA),
                    x_label = '', xlim = c(1,NA))
out$g
dev.off()

pdf(file = paste0(path,'FigS5_excess_ad_GH_tanDelta.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'UV Cure ',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Loss_Factor_1',
                    color_variable = 'condition',
                    cmap = col_pal[c(1,3)],
                    dot_size = 2,
                    y_scale = 'linear', ylim = c(0,NA),
                    x_label = '', xlim = c(1,NA))
out$g
dev.off()

pdf(file = paste0(path,'FigS5_excess_ad_GH_relaxation.pdf'), width = 4, height = 3)
out$g <- plot_stress_relax(df, selected_measurement = 'Measurement 2',
                       normalize=TRUE,
                       plot_mean=TRUE, plot_sd = TRUE,
                       color_variable = 'condition',
                       cmap = col_pal[c(1,3)],
                       xlim = c(0,30), x_scale = 'linear',
                       ylim = c(0,1), y_scale='linear')
out$g
dev.off()


##################################################################
# Figure 3. UV Cure GH diff MWs
##################################################################
f <- 'polymer_f_gh.csv'
dat <- readData(paste0(path,f), meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta$condition <- meta$V2
meta$condition[meta$condition=='4armPEG'] = 'Nor4PEG'
meta$condition[meta$condition=='10kNorHA'] = 'Nor8HA'
meta$condition[meta$condition=='40kNorHA'] = 'Nor18HA'
meta$condition[meta$condition=='60kNorHA'] = 'Nor40HA'
meta$condition <- factor(meta$condition, levels = c('Nor4PEG','Nor8HA','Nor18HA','Nor40HA'))

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'Fig3_VE_curing_diff_mws.pdf'), width = 7.8, height = 2.25)
out$g <- plot_time_sweep(df, selected_measurement = 'UV Cure',
                     facet_col='condition',
                     color_variable = 'condition',
                     cmap = col_pal,
                     plot_mean=TRUE,
                     y_scale='log10',
                     highlight_intervals = list(c(10,10+300)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'Fig3_VE_stress_relaxation.pdf'), width = 3.5, height = 2)#, units = 'in')#, res = 300)
out <- plot_stress_relax(df, selected_measurement = 'Stress Relaxation',
                       normalize=TRUE,
                       plot_mean=TRUE, plot_sd = TRUE,
                       color_variable = 'condition',
                       cmap = col_pal,
                       xlim = c(-1,30), x_scale = 'linear',
                       ylim = c(0,NA), y_scale='linear')
out$g
dev.off()

##################################################################
# Figure S12. UV Cure NorHA cysteine
##################################################################
f <- 'polymer_f_cysteine.csv'
dat <- readData(paste0(path,f), meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta$V2[meta$V2=='4armPEG'] = 'Nor4PEG'
meta$V2[meta$V2=='10kNorHA'] = 'Nor8HA'
meta$V2[meta$V2=='40kNorHA'] = 'Nor18HA'
meta$V2[meta$V2=='60kNorHA'] = 'Nor40HA'
meta$V2 <- factor(meta$V2, levels = c('Nor4PEG','Nor8HA','Nor18HA','Nor40HA'))
meta$V5 <- factor(meta$V5, levels = c('0x',"0.25x","0.325x","0.35x","0.375x","0.5x",'1x'))
meta$condition <- paste0(meta$V2,'_',meta$V5)

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'all_cysteine_curing_diff_mws.pdf'), width = 8, height = 10)#, units = 'in')#, res = 300)
out <- plot_time_sweep(df, selected_measurement = 'UV Cure',
                      color_variable = 'V2',
                      cmap = col_pal,
                      plot_mean=TRUE,
                      facet_col='V2',
                      facet_row = 'V5',
                      y_scale='log10',
                      highlight_intervals = list(c(10,10+300)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'FigS15_cysteine_curing_diff_mws.pdf'), width = 5, height = 3)#, units = 'in')#, res = 300)
out <- plot_time_sweep(df[df$V5 %in% c('0x','0.25x','0.5x') & df$V2 !='Nor4PEG',],
                       selected_measurement = 'UV Cure',
                       color_variable = 'V2',
                       cmap = col_pal,
                       plot_mean=TRUE,
                       facet_col='V2',
                       facet_row = 'V5',
                       y_scale='log10', ylim = c(.1,1000),
                       xlim = c(5,NA),
                       highlight_intervals = list(c(10,10+300)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'FigS15_cysteine_curing_plateau_storage.pdf'), width = 3, height = 3)
out <- plot_plateau(df[df$V5 %in% c('0x','0.25x','0.5x') & df$V2 !='Nor4PEG',],
                    selected_measurement = 'UV Cure',
                    color_variable = 'V2',
                    cmap = col_pal[2:4],
                    xaxis_column_name = 'V2',
                    facet_row='V5',
                    y_scale='linear')
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_x_discrete(labels = c(expression('Nor'[8]*'HA'),
                              expression('Nor'[18]*'HA'),expression('Nor'[40]*'HA')))+
  ylim(0,650)
dev.off()



##################################################################
# Figure 4. UV Cure same formulation
##################################################################
f <- 'norb_homopolymerization_same_formulations.csv'
dat <- readData(paste0(path,f), meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta$condition <- meta$V2
meta$condition[meta$condition=='4armPEG'] = 'Nor4PEG'
meta$condition[meta$condition=='10kNorHA'] = 'Nor8HA'
meta$condition[meta$condition=='40kNorHA'] = 'Nor18HA'
meta$condition[meta$condition=='60kNorHA'] = 'Nor40HA'
meta$condition <- factor(meta$condition, levels = c('Nor4PEG','Nor8HA','Nor18HA','Nor40HA'))

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'Fig4_same_formulation_mean_curing.pdf'), width = 10, height = 2.5)
out <- plot_time_sweep(df, selected_measurement = 'UV Cure',
                     facet_col='condition',
                     color_variable = 'condition',
                     cmap = col_pal,
                     plot_mean=TRUE,
                     ylim = c(1,NA), y_scale='log10',
                     highlight_intervals = list(c(10,10+300)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'Fig4_same_formulation_plateau_storage.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'UV Cure',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Storage_Modulus_Pa',
                    color_variable = 'condition',
                    cmap = col_pal,
                    y_scale = 'linear')
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

p_mat <- ggplot_build(out$g)$data[[1]][,c('y','group')] %>%
  mutate(group=factor(group))%>%
  group_map(~ aov(y ~ group, .x))
summary(p_mat[[1]])
TukeyHSD(p_mat[[1]])$group

pdf(file = paste0(path,'Fig4_same_formulation_stress_relaxation.pdf'), width = 3, height = 3)#, units = 'in')#, res = 300)
out <- plot_stress_relax(df, selected_measurement = 'Stress Relaxation',
                       normalize=TRUE,
                       plot_mean=TRUE, plot_sd = TRUE,
                       color_variable = 'condition',
                       cmap = col_pal,
                       xlim = c(-1,30), x_scale = 'linear',
                       ylim = c(0,NA), y_scale='linear')
out$g
dev.off()

p_mat <- out$filtered_df %>%
  filter(measurement=='Stress Relaxation') %>%
  group_by(sample,condition) %>%
  mutate(Norm_Shear_Stress=as.numeric(Norm_Shear_Stress))%>%
  do(tail(.,n=1)) %>%
  ungroup() %>%
  group_map(~ aov(Norm_Shear_Stress ~ condition, .x))
summary(p_mat[[1]])
TukeyHSD(p_mat[[1]])$condition


##################################################################
# Figure 6 & S13. UV Cure same formulation
##################################################################
f <- 'norb_homopolymerization_mechanics_match.csv'
dat <- readData(paste0(path,f), meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta$condition <- meta$V2
meta$condition[meta$condition=='4armPEG'] = 'Nor4PEG'
meta$condition[meta$condition=='10kNorHA'] = 'Nor8HA'
meta$condition[meta$condition=='40kNorHA'] = 'Nor18HA'
meta$condition[meta$condition=='60kNorHA'] = 'Nor40HA'
meta$condition <- factor(meta$condition, levels = c('Nor4PEG','Nor8HA','Nor18HA','Nor40HA'))

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'Fig6_same_mechanics_mean_curing.pdf'), width = 10, height = 2.5)
out <- plot_time_sweep(df, selected_measurement = 'UV Cure',
                     facet_col='condition',
                     color_variable = 'condition',
                     cmap = col_pal,
                     plot_mean=TRUE,
                     ylim = c(1,NA), y_scale='log10',
                     highlight_intervals = list(c(10,10+300)))
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'Fig6_same_mechanics_plateau_storage.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'UV Cure',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Storage_Modulus_Pa',
                    color_variable = 'condition',
                    cmap = col_pal,
                    y_scale = 'linear')
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

p_mat <- ggplot_build(out$g)$data[[1]][,c('y','group')] %>%
  mutate(group=factor(group))%>%
  group_map(~ aov(y ~ group, .x))
summary(p_mat[[1]])
TukeyHSD(p_mat[[1]])$group

pdf(file = paste0(path,'Fig6_same_mechanics_loss_factor.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'UV Cure',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Loss_Factor_1',
                    color_variable = 'condition',
                    cmap = col_pal,
                    y_scale = 'linear')
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

pdf(file = paste0(path,'Fig6_same_mechanics_stress_relaxation.pdf'), width = 3, height = 3)
out <- plot_stress_relax(df, selected_measurement = 'Stress Relaxation',
                         normalize=TRUE,
                         plot_mean=TRUE, plot_sd = TRUE,
                         color_variable = 'condition',
                         cmap = col_pal,
                         xlim = c(-1,30), x_scale = 'linear',
                         ylim = c(0,NA), y_scale='linear')
out$g
dev.off()

p_mat <- out$filtered_df %>%
  filter(measurement=='Stress Relaxation') %>%
  group_by(sample,condition) %>%
  mutate(Norm_Shear_Stress=as.numeric(Norm_Shear_Stress))%>%
  do(tail(.,n=1)) %>%
  ungroup() %>%
  group_map(~ aov(Norm_Shear_Stress ~ condition, .x))
summary(p_mat[[1]])
TukeyHSD(p_mat[[1]])$condition


f <- 'mechanics_match_post_swelling.csv'
dat <- readData(paste0(path,f), meta_from_file_name = TRUE)

df <- dat@long
meta <- dat@metadata

meta$condition <- meta$V1
meta$condition[meta$condition=='peg'] = 'Nor4PEG'
meta$condition[meta$condition=='10k'] = 'Nor8HA'
meta$condition[meta$condition=='40k'] = 'Nor18HA'
meta$condition[meta$condition=='60k'] = 'Nor40HA'
meta$condition <- factor(meta$condition, levels = c('Nor4PEG','Nor8HA','Nor18HA','Nor40HA'))

df <- cbind(df, meta[match(df$sample,meta$sample),colnames(meta) != 'sample'])

pdf(file = paste0(path,'Fig6_same_mechanics_swollen_plateau_storage.pdf'), width = 3, height = 3)
out <- plot_plateau(df, selected_measurement = 'Time Sweep',
                    xaxis_column_name = 'condition',
                    yaxis_column_name = 'Storage_Modulus_Pa',
                    color_variable = 'condition',
                    cmap = col_pal,
                    y_scale = 'linear')
out$g+theme(axis.text.x = element_text(angle = 45,hjust=1))
dev.off()

p_mat <- ggplot_build(out$g)$data[[1]][,c('y','group')] %>%
  mutate(group=factor(group))%>%
  group_map(~ aov(y ~ group, .x))
summary(p_mat[[1]])
TukeyHSD(p_mat[[1]])$group

out <- plot_stress_relax(df, selected_measurement = 'Stress Relaxation',
                         normalize=TRUE,
                         plot_mean=TRUE, plot_sd = TRUE,
                         color_variable = 'condition',
                         cmap = col_pal,
                         xlim = c(-1,30), x_scale = 'linear',
                         ylim = c(0,NA), y_scale='linear')
out$g
