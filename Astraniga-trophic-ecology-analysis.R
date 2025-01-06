
# Astrangia trophic ecology analysis 
# Taylor Lindsay 
# DOI: 
# Date: Jan 2025

# SET UP ------------------------------------------------------------------

# set seed
set.seed(500)

#install packages 
library(knitr)       #for rmd functions 
library(tidyr)
library(tidyverse)
library("dplyr")
library(ggplot2)
library(ggpubr)
library(lubridate)    #environmental data 
library("survival")   #Kaplan-meier curves
library("survminer")  #survival curves 
library(ggfortify)    #pca plots 
library('gridExtra')  #arrange plots 
library('grid')       #arranging plots 
library(SIBER)        #SIBER plots  
library(rjags)        #bayesian mixing models 
#devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)     # sankey plots 
library(cowplot)      #arranging plots 
library(vegan)        #MANOVAs
library('ggnewscale')   # for SIBER plot scale 
library(scales)       # sum v 

# set working directory 
setwd('~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/')

### Set up comparisons and parameters for later graphs 

# define comparisons for boxplot t-tests 
treatment_comparisons <- list( c("apo_deep","sym_deep"), c("apo_control","sym_control"), c("apo_shade","sym_shade"), c("apo_control","apo_shade"), c("apo_control","apo_deep"), c("sym_control","sym_shade"), c("sym_control","sym_deep"))
treatment_comparisons2 <- list( c("apo_deep","sym_deep"), c("apo_control","sym_control"))
treatment_comparisons3 <- list(c("Host","Sym"))
treatment_comparisons4 <- list(c("APO","SYM"))
treatment_comparisons_combined <- list( c("control","shade"), c("control","deep"))

# define orders for different graphs 
measurement_order <- c('control','shade','deep') 
c_d_order <- c('apo_control','sym_control','apo_deep', 'sym_deep') 
x_order <- c('apo_control','sym_control','apo_shade', 'sym_shade', 'apo_deep', 'sym_deep') 
x_order_combined <- c('control','shade', 'deep') 
a_s_combined <- c('APO','SYM')
a_s_list <- list(c('APO','SYM'))

# define jitter width 
pj <- position_jitterdodge(jitter.width=0.2, seed=9,
                           jitter.height = 0)

# ENVIRONMENTAL: data prep--------------------------------------------------

# Master data set 
enviro <- read.csv('DATA/TLAP_results_enviro(C).csv') %>%
  dplyr::select(!X)

# make it datetime format 
enviro$datetime <- mdy_hm(enviro$datetime,tz=Sys.timezone())

colnames(enviro) <- c('datetime','temp.Control','light.Control','temp.Shade','light.Shade','temp.Deep','light.Deep') 

# save only earlier dates 
enviro <- enviro %>%
  filter(datetime <= "2023-08-08") %>%
  filter(datetime >= "2023-06-21")

# current means 
colMeans(enviro[,2:7], na.rm = TRUE)

# pivot longer 
long <- pivot_longer(enviro,cols = c(temp.Control,temp.Shade,temp.Deep,light.Control,light.Shade,light.Deep),
                     names_to = "treatment", values_to = "value") %>%
  separate(.,treatment, c("variable","treatment"))

#calculate mean and sd 
enviro_summary <- long %>%
  group_by(variable, treatment) %>%
  summarise(mean = mean(value), sd = sd(value))

# save summary file 
write.csv(enviro_summary, "STATS/TLAP_STATS_Enviro_means.csv", row.names = FALSE)

# reorganize data and get mean daily values 
long_temp <- long %>% filter(variable == "temp")
long_temp$datetime <- as.Date(long_temp$datetime)
daily_temp <-long_temp %>% group_by(datetime,treatment) %>%
  summarize(mean_temp = mean(value))
long_light <- long %>% filter(variable == "light")
long_light$datetime <- as.Date(long_light$datetime)
daily_light <-long_light %>% group_by(datetime,treatment) %>%
  summarize(mean_light = mean(value))

# ENVIRONMENTAL: Statistics (Table S1) --------------------------------

# write function 
perform_paired_t_test <- function(group1, group2, group_label) {
  t_test <- t.test(group1, group2, paired = TRUE)
  data.frame(
    Comparison = group_label,
    t = t_test$statistic,
    df = t_test$parameter,
    p = formatC(t_test$p.value, format = "e", digits = 3)
  )
}

# Perform the t-tests and collect results
test_results_list <- list(
  perform_paired_t_test(enviro$temp.Control, enviro$temp.Shade, "Temp Control vs. Shade"),
  perform_paired_t_test(enviro$temp.Control, enviro$temp.Deep, "Temp Control vs. Deep"),
  perform_paired_t_test(enviro$light.Control, enviro$light.Shade, "Light Control vs. Shade"),
  perform_paired_t_test(enviro$light.Control, enviro$light.Deep, "Light Control vs. Deep")
)

# Combine all results into one data frame
all_enviro <- do.call(rbind, test_results_list)

write.csv(all_enviro, "STATS/TLAP_STATS_Enviro_paired-t-test.csv", row.names = FALSE)

# ENVIRONMENTAL: Graph temp & light (Fig 2------------------------------------

# Temperature plot
daily_temp_plot <- ggplot(daily_temp, aes(x=datetime,y=mean_temp, color=treatment)) +
  geom_line(size=1.5) +
  labs(y= "Mean Daily Temperature (˚C)", x = "", color="Treatment") +
  theme_bw() + 
  theme(text = element_text( size=20),legend.position = c(.5, .25), 
        legend.box.background = element_rect(color="black", size=1.5), 
        axis.text.x = element_blank()) +
  scale_color_manual(values = c("#C0C0E1", "#210124", "#8A4594"))

# Light plot 
daily_light_plot <- ggplot(daily_light, aes(x=datetime,y=mean_light, color=treatment)) +
  geom_line(size=1.5) +
  labs(y= "Mean Daily Light (lum)", x = "Date", color="Treatment")  +
  theme_bw() + 
  theme(text = element_text( size=20), legend.position = "none") +
  scale_color_manual(values = c("#C0C0E1", "#210124", "#8A4594"))

#daily_temp_plot
#daily_light_plot

temp1 <- ggplotGrob(daily_temp_plot)
light1 <- ggplotGrob(daily_light_plot)

#combine ecotype and treatment 
enviro_arrange <- grid.arrange(temp1, light1, nrow=2,
                               bottom = grobTree(
                                 textGrob("A", x = 0.08, y = 97, just = "left", gp = gpar(fontsize = 18, fontface = "bold")),
                                 textGrob("B", x = 0.1, y = 47, just = "right", gp = gpar(fontsize = 18, fontface = "bold"))))

#save graphs 
ggsave("TLAP_FIG_2_Enviro.jpg", plot = enviro_arrange, path = 'FIGURES/', width = 10, height = 10)


# SURVIVAL: data prep & Graph (Fig 3) --------------------------------------------------------------

# load data and filter out only dead and alive individuals 
raw_s <- read.csv('DATA/TLAP_results_survival.csv') %>%
  filter(survival_10.10.23 == 0 | survival_10.10.23 == 1) %>%
  filter(end_status != "Missing")

# make the final column numeric 
raw_s$survival_10.10.23 <- as.numeric(raw_s$survival_10.10.23)

# create a surv object 
surv_object <- Surv(raw_s$days_alive, raw_s$survival_10.10.23)

# plot by treatment 
fit1 <- survfit(surv_object ~ treatment, data = raw_s)
plot_treatment <- ggsurvplot(fit1, data = raw_s, 
                             #pval = TRUE, 
                             #pval.coord = c(20, 0.905),
                             size=1.5,
                             legend = c(.25, 0.15), 
                             legend.title = "Treatment", 
                             legend.labs = c("Control", "Deep", "Shade"), 
                             palette = c("#C0C0E1", "#210124", "#8A4594"), 
                             ylim = c(0.91, 1),
                             ylab = c(""), 
                             ggtheme = theme_classic2(base_size=20)
)
plot_treatment <- plot_treatment +  xlab("Days")
plot_treatment$plot <- plot_treatment$plot + theme(axis.text.y = element_blank(), legend.box.background = element_rect(color="black", size=1.5) )
#plot_treatment

# plot by apo/sym
fit2 <- survfit(surv_object ~ apo_sym, data = raw_s)
plot_sym_apo <- ggsurvplot(fit2, data = raw_s, 
                           #pval = TRUE, 
                           #pval.coord = c(20, 0.905),  
                           size = 1.5, 
                           legend = c(.3, 0.15), 
                           legend.title = "Initial ecotype", 
                           legend.labs = c("Aposymbiotic", "Symbiotic"), 
                           palette =c("#e4d2ba", "#724a29"), 
                           ylim = c(0.91, 1), 
                           ggtheme = theme_classic2(base_size=20))
plot_sym_apo <- plot_sym_apo + xlab("Days")
plot_sym_apo$plot <- plot_sym_apo$plot + theme(legend.box.background = element_rect(color="black", size=1.5) )
#plot_sym_apo

surv1 <- ggplotGrob(plot_sym_apo$plot)
surv2 <- ggplotGrob(plot_treatment$plot)

#combine ecotype and treatment 
survival_arrange <- grid.arrange(surv1, surv2, nrow=1,
                                 bottom = grobTree(
                                   textGrob("A", x = 0.105, y = 78, just = "left", gp = gpar(fontsize = 18, fontface = "bold")),
                                   textGrob("B", x = 0.57, y = 78, just = "right", gp = gpar(fontsize = 18, fontface = "bold"))))

#save graphs 
ggsave("TLAP_FIG_3_surival.jpg", plot = survival_arrange , path = 'FIGURES/', width = 12, height = 8)

# SURVIVAL: Statistics ----------------------------------------------------

# Define a function to perform a log-rank test and extract the results
perform_log_rank_test <- function(surv_object, grouping_variable, data) {
  # Perform the log-rank test
  log_rank_test <- survdiff(surv_object ~ get(grouping_variable), data = data)
  
  # Extract test statistic, degrees of freedom, and p-value
  test_statistic <- log_rank_test$chisq
  degrees_of_freedom <- length(log_rank_test$n) - 1
  p_value <- 1 - pchisq(test_statistic, df = degrees_of_freedom)
  
  # Return a data frame with the results
  data.frame(
    Grouping_Variable = grouping_variable,
    chisq = test_statistic,
    df = degrees_of_freedom,
    p = formatC(p_value, digits = 3)
  )
}

# Example: List of survival objects and grouping variables
tests <- list(
  list(surv_object = surv_object, grouping_variable = "treatment", data = raw_s),
  list(surv_object = surv_object, grouping_variable = "apo_sym", data = raw_s)
)

# Apply the function to all tests
all_log_rank_results <- do.call(rbind, lapply(tests, function(test) {
  perform_log_rank_test(test$surv_object, test$grouping_variable, test$data)
}))

# Save 
write.csv(all_log_rank_results, "STATS/TLAP_STATS_survival_kaplan-meier.csv", row.names = FALSE)

# How many from each group survived

# Function to count alive and dead
count_alive_dead <- function(data, filter_condition = NULL, remove_na_apo_sym = FALSE, remove_na_treatment = FALSE) {
  # Remove rows with NA in the apo_sym column if specified
  if (remove_na_apo_sym) {
    data <- subset(data, !is.na(apo_sym))}
  
  # Remove rows with NA in the apo_sym column if specified
  if (remove_na_treatment) {
    data <- subset(data, !is.na(treatment))}
  
  # Apply filter condition if specified
  if (!is.null(filter_condition)) {
    data <- subset(data, eval(parse(text = filter_condition)))}
  
  # Count occurrences of "alive" and "dead"
  counts <- table(data$end_status)
  
  # Ensure counts default to 0 if "alive" or "dead" are missing
  alive_count <- ifelse("Alive" %in% names(counts), counts["Alive"], 0)
  dead_count <- ifelse("Dead" %in% names(counts), counts["Dead"], 0)
  
  # Create a data frame with the results
  result <- data.frame(
    Category = ifelse(is.null(filter_condition), "All Data", filter_condition),
    Alive = alive_count,
    Dead = dead_count
  )
  return(result)
}

# Count for each category
count_results <- rbind(
  count_alive_dead(raw_s),  # 1. All data
  count_alive_dead(raw_s, "treatment == 'control'", remove_na_treatment = TRUE),  # 2. Treatment = control
  count_alive_dead(raw_s, "treatment == 'shade'", remove_na_treatment = TRUE),  # 3. Treatment = shallow
  count_alive_dead(raw_s, "treatment == 'deep'", remove_na_treatment = TRUE),  # 4. Treatment = deep
  count_alive_dead(raw_s, "apo_sym == 'APO'", remove_na_apo_sym = TRUE),  # 5. Apo_sym = APO
  count_alive_dead(raw_s, "apo_sym == 'SYM'", remove_na_apo_sym = TRUE)  # 6. Apo_sym = SYM
)

# add total column 
count_results <- count_results %>%
  mutate(total = (Alive + Dead))

# Save 
write.csv(count_results, "STATS/TLAP_STATS_survival_counts.csv", row.names = FALSE)


# ECOTYPE SWITCHING: analysis & graph (Fig 4) -------------------------------------


# read phys data 
master <- read.csv('DATA/TLAP_results_meta_phys.csv') %>%
  filter(full_treatment != "NA")

# create final ecotype categories
master <- master %>% 
  mutate(new = if_else(sym.cm2 > 500000, "SYM", "APO"))

# sankey plot 
sankey_data <- master %>%
  filter(!is.na(sym.cm2)) %>%
  mutate(full = paste0(Apo_Sym,"_",new)) %>%
  mutate(full_start = paste0(treatment, "_", Apo_Sym)) %>%
  mutate(full_end = paste0(treatment, "_", new)) 

### one at a time 

# CONTROL 

# filter for control data 
df_cont <- sankey_data %>% filter(treatment == "control")

# Table for percents that switched
sum_cont <- sankey_data %>% group_by(treatment) %>% count(full)

# count the total 
Totalcontrol = nrow(df_cont)

# make a sankey object 
df_control <- df_cont %>% make_long(full_start, full_end) %>% 
  mutate(full_percent= paste0(x, node))

# count each 
dagg_control <- df_control%>%
  dplyr::group_by(full_percent, node)%>%
  tally()

# calculate percents and then generate label column 
dagg_control <- dagg_control%>%
  group_by(full_percent, node)%>%
  mutate(pct = signif((n / Totalcontrol) * 100, 2)) %>%
  separate(node, c("treatment", "ecotype"), remove = FALSE) %>%
  mutate(per = paste0(ecotype, " ", pct, "%" ) )

# merge new labels back to sankey object 
df2_control <- merge(df_control, dagg_control, all.x = TRUE)

sankey_control <- ggplot(df2_control, aes(x = x, 
                                          next_x = next_x, 
                                          node = node, 
                                          next_node = next_node,
                                          fill = factor(node),
                                          label = ecotype,
)) +
  geom_sankey(flow.alpha = 0.7, # this changes the pigmentation of the flow 
              node.color="black") + # this node color changes the outline 
  scale_fill_manual(values = c("#e4d2ba", "#8a6136")) + 
  theme_sankey(base_size = 16) + 
  theme(#text = element_text(size=25), 
    legend.position="none",
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank()) +
  geom_sankey_label(size = 6, color = "black", angle = 90) + #hjust = .1
  # coord_fixed(ratio = .008) + 
  # Add a colored box using annotation_custom
  annotate("rect", xmin = 0.8, xmax = 0.9, ymin = -Inf, ymax = Inf, fill = "#C0C0E1", color="black") +  # Box properties
  annotate("text", x = .85, y = 0, label = "Control", size = 6, angle = 90, color = "black", hjust = 0.5)

# SHADE 

# filter for control data 
df_shade <- sankey_data %>% filter(treatment == "shade")

# count the total 
Totalshade = nrow(df_shade)

# make a sankey object 
df_shade <- df_shade %>% make_long(full_start, full_end) %>%
  mutate(full_percent = paste0(x, node))

# count each 
dagg_shade <- df_shade%>%
  dplyr::group_by(full_percent, node) %>%
  tally()

# calculate percents and then generate label column 
dagg_shade <- dagg_shade%>%
  group_by(full_percent, node)%>%
  mutate(pct = signif((n / Totalshade) * 100, 2)) %>%
  separate(node, c("treatment", "ecotype"), remove = FALSE) %>%
  mutate(per = paste0(ecotype, " ", pct, "%" ) )

# merge new labels back to sankey object 
df2_shdae <- merge(df_shade, dagg_shade, all.x = TRUE)

sankey_shade <- ggplot(df2_shdae , aes(x = x, 
                                       next_x = next_x, 
                                       node = node, 
                                       next_node = next_node,
                                       fill = factor(node),
                                       label = ecotype,
)) +
  geom_sankey(flow.alpha = 0.7, # this changes the pigmentation of the flow 
              node.color="black") + # this node color changes the outline 
  scale_fill_manual(values = c("#e4d2ba", "#8a6136")) + 
  theme_sankey(base_size = 16) + 
  theme(#text = element_text(size=25), 
    legend.position="none",  
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank()) +
  geom_sankey_label(size = 6, color = "black", angle = 90) + #hjust = .1
  #coord_fixed(ratio = .008) + 
  # Add a colored box using annotation_custom
  annotate("rect", xmin = 0.8, xmax = 0.9, ymin = -Inf, ymax = Inf, fill = "#8A4594", color="black") +  # Box properties
  annotate("text", x = .85, y = 0, label = "Shade", size = 6, angle = 90, color = "white", hjust = 0.5)
sankey_shade

# DEEP 

# filter for deep data 
df_deep <- sankey_data %>% filter(treatment == "deep")

# count the total 
Totaldeep = nrow(df_deep)

# make a sankey object 
df_deep <- df_deep %>% make_long(full_start, full_end) %>%
  mutate(full_percent = paste0(x, node))

# count each 
dagg_deep<- df_deep%>%
  dplyr::group_by(full_percent, node) %>%
  tally()

# calculate percents and then generate label column 
dagg_deep <- dagg_deep%>%
  group_by(full_percent, node)%>%
  mutate(pct = signif((n / Totaldeep) * 100, 2)) %>%
  separate(node, c("treatment", "ecotype"), remove = FALSE) %>%
  mutate(per = paste0(ecotype, " ", pct, "%" ) )

# merge new labels back to sankey object 
df2_deep <- merge(df_deep, dagg_deep, all.x = TRUE)

sankey_deep <- ggplot(df2_deep, aes(x = x, 
                                    next_x = next_x, 
                                    node = node, 
                                    next_node = next_node,
                                    fill = factor(node),
                                    label = ecotype,
)) +
  geom_sankey(flow.alpha = 0.7, # this changes the pigmentation of the flow 
              node.color="black") + # this node color changes the outline 
  scale_fill_manual(values = c("#e4d2ba", "#8a6136")) + 
  theme_sankey(base_size = 16) + 
  theme(#text = element_text(size=25), 
    legend.position="none",  
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    axis.title.x=element_blank()) +
  geom_sankey_label(size = 6, color = "black", angle = 90) + #hjust = .1
  #coord_fixed(ratio = .008) + 
  # Add a colored box using annotation_custom
  annotate("rect", xmin = 0.8, xmax = 0.9, ymin = -Inf, ymax = Inf, fill = "#210124", color="white") +  # Box properties
  annotate("text", x = .85, y = 0, label = "Deep", size = 6, angle = 90, color = "white", hjust = 0.5)

sankey_control
sankey_shade
sankey_deep

#combine ecotype and treatment 
all <- grid.arrange(sankey_control, sankey_shade, sankey_deep, nrow=3,
                    bottom = grobTree(
                      textGrob("Initial", x = 0.255, y = 2, just = "left", gp = gpar(fontsize = 18)),
                      textGrob("Final", x = 0.745, y = 2, just = "right", gp = gpar(fontsize = 18)),
                      textGrob("", x = .59, y = 1.5, just = "right", gp = gpar(fontsize = 20)),
                      # add percents 
                      textGrob("72%", x = .66, y = 93.5, just = "right", gp = gpar(fontsize = 16)),
                      textGrob("31%", x = .66, y = 87.5, just = "right", gp = gpar(fontsize = 16)),
                      textGrob("28%", x = .66, y = 80.6, just = "right", gp = gpar(fontsize = 16)),
                      textGrob("69%", x = .66, y = 74.5, just = "right", gp = gpar(fontsize = 16)),
                      
                      textGrob("23%", x = .66, y = 63, just = "right", gp = gpar(fontsize = 16)),
                      textGrob("78%", x = .66, y = 53, just = "right", gp = gpar(fontsize = 16)),
                      textGrob("100%", x = .66, y = 43, just = "right", gp = gpar(fontsize = 16)),
                      
                      textGrob("100%", x = .66, y = 23, just = "right", gp = gpar(fontsize = 16)),
                      textGrob("100%", x = .66, y = 12, just = "right", gp = gpar(fontsize = 16))
                    ))

ggsave("TLAP_FIG_4_ecotype_switch.jpg", all, path = 'FIGURES/', height=10, width = 12)


# PHYSIOLOGY:  ------------------------------------------------------------

# create new columns 
master <- master %>%
  mutate(log_chl_sym = log(chl_ug.sym)) %>%
  mutate(log_chl_sym = ifelse(is.infinite(log_chl_sym), NA, log_chl_sym))

# Graph sym x treatment 
AP23_sym_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=sym.cm2, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,0.5),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) + 
  #specifics 
  #scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #ylim(0, 1900000) +
  scale_y_continuous(
    limits = c(0, 1900000),  # Set y-axis limits
    labels = function(x) {scales::label_scientific()(signif(x, 1))} ) + 
  #labs(x = "", y = expression(atop("Symbiont Density", paste("(cells/cm"^2, ")")))) + 
  labs(x = "", y = bquote("Symbiont Density (cells/" ~ cm^2 ~ ")")) 

# Graph Sym x ecotype 
AP23_sym_ecotype <- ggplot(master,aes(x=Apo_Sym, y=sym.cm2, fill=Apo_Sym)) + 
  #DATA 
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  #ylim(0, 1900000)+
  scale_y_continuous(limits = c(0, 1900000))+ 
  stat_compare_means(comparisons = a_s_list, method = "wilcox.test", size=8) + 
  #specifics 
  #scale_x_discrete(labels=c('Apo','Sym') ) +
  scale_fill_manual(values = c("#e4d2ba", "#724a29"), labels=c('Apo', 'Sym')) + 
  labs(x = "", y = "") 

# Graph CHL x treatment 
AP23_chl_treat <- ggplot(master,aes(x=factor(treatment, level=x_order_combined), y=chla.ug.cm2, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,0.5),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) +
  #specifics 
  #scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "", y = expression(atop("Chlorophyll a", paste("(ug/cm"^2, ")")))) + 
  labs(x = "", y = bquote("Chlorophyll a (ug/" ~ cm^2 ~ ")")) + 
  ylim(0, 3.7)

# Graph Chl x ecotype 
AP23_chl_ecotype <- ggplot(master,aes(x=sym.cm2, y=chla.ug.cm2)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "", y = "") +
  ylim(0, 3.7)

# Graph chl/sym by sym 
AP23_chl_sym_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=log_chl_sym, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() + 
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.ticks.x= element_blank(),
        axis.text.x = element_blank(),
        plot.margin=unit(c(0,0,0,0.5),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) +  
  #specifics 
  scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "", y = expression(atop("Chlorophyll per Symbiont", paste("(ug/cell)")))) +
  labs(x = "", y = bquote("Chlorophyll per Symbiont (log(ug/cell))")) +
  ylim(-15.5, -8)

# Chl/sym by sym
AP23_chl_sym_ecotype <- ggplot(master,aes(x=sym.cm2, y=log_chl_sym)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "", y = "") +
  ylim(-15.5, -8)


# Graph  calice density by treatment 
AP23_cal_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=calice_density, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_text(size=30),
        plot.margin=unit(c(0,0,0,0.5),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=10) + 
  #specifics 
  scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "Treatment", y = expression(atop("Calice Density", paste("(Polyps per cm2)")))) +
  labs(x = "Treatment", y = bquote("Calice Density (Polyps per" ~ cm^2 ~ ")")) + 
  ylim(0.5,5)

# calice density by treatment 
AP23_cal_ecotype <- ggplot(master,aes(x=sym.cm2, y=calice_density)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        #axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "Symbiont Density", y = "") +
  ylim(0.5,5)

# Graph antioxidants x treatment 
AP23_AO_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=cre.umol.mgprot, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) + 
  #specifics 
  #scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "", y = expression(atop("Antioxidants/protein", paste("(umol/mg)")))) + 
  labs(x = "", y = bquote("Antioxidants per protein (umol/mg)")) + 
  ylim(0,1.6)

# Graph antioxidants x sym
AP23_AO_ecotype <- ggplot(master,aes(x=sym.cm2, y=cre.umol.mgprot)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "", y = "")+ 
  ylim(0,1.6)

# Graph protein x treatment 
AP23_prot_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=prot_mg.cm2, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        #plot.margin=unit(c(0,0,0,1),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) + 
  #specifics 
  #scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "", y = expression(atop("Soluble Protein", paste("(mg/cm"^2, ")")))) + 
  labs(x = "", y = bquote("Soluble Protein (mg/" ~ cm^2~ ")")) + 
  ylim(0,3.8)

# Graph prot x syms 
AP23_prot_ecotype <- ggplot(master,aes(x=sym.cm2, y=prot_mg.cm2)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "", y = "")+ 
  ylim(0,3.8)

# Graph AFDW x treatment 
AP23_AFDW_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=AFDW.mg.cm2, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text(size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        #plot.margin=unit(c(0,0,0,1),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) + 
  #specifics 
  #scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "", y = expression(atop("Biomass", paste("(mg/cm"^2, ")")))) +
  labs(x = "", y = bquote("Biomass (mg/" ~ cm^2~ ")")) + 
  ylim(0.001,0.006)

AP23_AFDW_ecotype <- ggplot(master,aes(x=sym.cm2, y=AFDW.mg.cm2)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "", y = "")+ 
  ylim(0.001,0.006)

# Graph lipids x treatment 
AP23_lip_treat <- ggplot(master, aes(x=factor(treatment, level=x_order_combined), y=lipids.mg.cm2, fill=factor(treatment, level=measurement_order))) + 
  #DATA
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=1, show.legend = FALSE)+
  #AESTHETICS
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        axis.text.y = element_text(size=30), 
        axis.text.x = element_text(size=30),
        #plot.margin=unit(c(0,0,0.5,1),"cm")
  ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=8) + 
  #specifics 
  scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  scale_fill_manual(values = c("#C0C0E1", "#9e7bb5", "#3e2f84")) + 
  #labs(x = "Treatment", y = expression(atop("Total Lipids", paste("(mg/cm"^2, ")")))) + 
  labs(x = "Treatment", y = bquote("Total Lipids (mg/" ~ cm^2 ~ ")")) + 
  ylim(0,0.003)

AP23_lip_ecotype <- ggplot(master,aes(x=sym.cm2, y=lipids.mg.cm2)) + 
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_vline(xintercept = 500000, linetype="dashed") + 
  geom_point( aes(color=Apo_Sym), size=4, show.legend = FALSE) +
  # AESTHETICS 
  theme_bw() +
  theme(legend.position="none", 
        text = element_text( size=40), 
        #axis.text.x = element_blank(), 
        axis.ticks.x= element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y= element_blank(), 
        plot.margin=unit(c(0,1.5,0,0),"cm")
  ) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  #specifics 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  labs(x = "Symbiont Density", y = "")+ 
  ylim(0,0.003)

# save the treatment and ecotype plots into objects 
g_treat1 <- plot_grid(AP23_sym_treat, AP23_chl_treat,AP23_chl_sym_treat, AP23_cal_treat, 
                      ncol = 1, align = "v", 
                      labels = c("A", "E", "I", "M"),  label_size = 35, label_x = 0.36, label_y = 0.99)
g_eco1 <-plot_grid( AP23_sym_ecotype,AP23_chl_ecotype,AP23_chl_sym_ecotype, AP23_cal_ecotype, 
                    ncol = 1,
                    labels = c("B", "F", "J", "N"),  label_size = 35, label_x = 0.08, label_y = 0.99)
g_treat2 <-plot_grid(AP23_AO_treat, AP23_prot_treat,AP23_AFDW_treat, AP23_lip_treat, 
                     ncol = 1, align = "v", 
                     labels = c("C", "G", "K", "O"),  label_size = 35, label_x = 0.33, label_y = 0.99)
g_eco2 <-plot_grid(AP23_AO_ecotype, AP23_prot_ecotype,AP23_AFDW_ecotype,AP23_lip_ecotype,
                   ncol = 1, align = "v",
                   labels = c("D", "H", "L", "P"),  label_size = 35, label_x = 0.08, label_y = 0.99)

# combine 
all_metab <-plot_grid(g_treat1, g_eco1, g_treat2, g_eco2, 
                      ncol = 4, rel_widths = c(1, 1.3, 1, 1.3))

#save 2x6 graphs 
ggsave("TLAP_Fig_5_phys.jpg", plot = all_metab, path = 'FIGURES/', width = 30, height = 40)






# ISOTOPES: data & basic plots (Fig S1) ----------------------------------------------------------------

# load raw isotope data 
raw <- read.csv('~/Desktop/GITHUB/TL_Astrangia/TLAP_CSIA/TLAP_CSIA_results.csv') %>%
  filter(remove == "N")

metadata <- master %>%
  dplyr::select(c(sample_id, Apo_Sym, treatment,full_treatment,sym.cm2))

# unify with metadata 
raw2 <-  merge(metadata, raw, by="sample_id") %>%
  dplyr::select(full_id, sample_id, round, vial, element, fraction, Apo_Sym, treatment, full_treatment, sym.cm2, Ala, Gly, Thr, Ser, Val, Leu, Ile, Nle, Pro, Asx, Met, Glx, Phe, Lys, Arg) 

#write.csv(file = "STATS/TLAP_ALL_Isotope_Results.csv", raw2)

# reorder treatments 
raw2$treatment <- factor(raw2$treatment, levels = c("control", "shade", "deep"))

# pivot data to be longer 
c_long <- raw2 %>% 
  filter(element == "C")%>%
  pivot_longer(cols = c(Ala, Gly, Thr, Ser, Val, Leu, Ile, Pro, Asx, Met, Glx, Phe, Lys), # no Nle or Arg
               names_to = "amino_acid", values_to = "delt") 

# Carbon EAAs only
c_long_EAAs <- c_long %>%
  filter(amino_acid %in% c("Ile", "Lys", "Leu", "Phe", "Thr", "Val"))

n_long <- raw2 %>% 
  filter(element == "N") %>%
  pivot_longer(cols = c(Ala, Gly, Thr, Ser, Val, Leu, Ile, Pro, Asx, Met, Glx, Phe, Lys, Arg),# no Nle
               names_to = "amino_acid", values_to = "delt") 

# reorder the AAs for the x axis 
C_order <- c("Gly","Ser","Asx","Glx","Pro","Ala","Thr","Lys","Ile","Met","Val","Phe","Leu") 
N_order <- c("Glx","Asx","Pro","Ala","Leu","Ile","Val","Phe","Met","Lys","Arg","Gly","Ser","Thr") 

custom_c_colors <- c("Gly" = "dimgray", "Ser"= "dimgray","Asx"= "dimgray", "Glx"= "dimgray", "Pro"= "dimgray", "Ala"= "dimgray", "Thr" = "white", "Lys"= "white", "Ile"= "white", "Met"= "white", "Val"= "white", "Phe"= "white", "Leu"= "white")


#Carbon boxplot x fraction
c_boxplot_frac <- ggplot(c_long) +
  geom_boxplot(aes(x=factor(amino_acid, level=C_order),y=delt, fill=amino_acid))+
  labs(y= ~ δ^13 ~ "C (‰)", x = "Amino Acid", color="") +
  scale_fill_manual(values = custom_c_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x=element_text(size=25),
    axis.title.y=element_text(size=25),
    plot.margin=unit(c(0,0,0.5,0),"cm")
  )
#c_boxplot_frac

custom_n_colors <- c("Gly" = "white", "Ser"= "white","Asx"= "dimgray", "Glx"= "dimgray", "Pro"= "dimgray", "Ala"= "dimgray", "Thr" = "gray", "Lys"= "white", "Ile"= "dimgray", "Met"= "white", "Val"= "dimgray", "Phe"= "white", "Leu"= "dimgray", "Arg" = "white")

# N BASIC GGPLOT 
N_boxplot_frac <- ggplot(n_long) +
  geom_boxplot(aes(x=factor(amino_acid, level=N_order),y=delt, fill=amino_acid))+
  scale_fill_manual(values = custom_n_colors) +
  labs(y=  ~ δ^15 ~ "N (‰)", x = "Amino Acid", color="") +
  theme_bw() + 
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=20),
    axis.text.y=element_text(size=20),
    axis.title.x=element_text(size=25),
    axis.title.y=element_text(size=25),
    plot.margin=unit(c(0,0,0,0),"cm")
  )
#N_boxplot_frac

# ARRANGE 

supp_arrange <- grid.arrange(c_boxplot_frac, N_boxplot_frac, nrow=2, 
                             bottom = grobTree(
                               textGrob("A", x = 0.08, y = 146, just = "left", gp = gpar(fontsize = 30, fontface = "bold")),
                               textGrob("B", x = 0.1, y = 71, just = "right", gp = gpar(fontsize = 30, fontface = "bold"))
                             ))

ggsave("TLAP_FIG_S1_Isotopes.jpg", plot=supp_arrange, path = "FIGURES/",width = 15, height = 15)

# ISOTOPES: amino acid stats  ---------------------------------------------

### SOURCE AMINO ACID MANOVAs

# SOURCE AAs x Fraction 
n_source_frac <- raw2 %>% 
  filter(element == "N") %>%
  select("fraction", "Gly", "Lys", "Phe", "Ser")
matrix1<- as.matrix(n_source_frac[,2:5]) # response variables
manova1<-manova(matrix1~as.factor(fraction), data=n_source_frac)
x <- summary(manova1) # p = < 0.0001

# SOURCE AAs x Ecotype 
n_source_ecotype <- raw2 %>% 
  filter(element == "N") %>%
  select("Apo_Sym", "Gly", "Lys", "Phe", "Ser")
matrix2<- as.matrix(n_source_ecotype[,2:5]) # response variables
manova2<-manova(matrix2~as.factor(Apo_Sym), data=n_source_ecotype)
summary(manova2) # p = 0.16

# SOURCE AAs x treatment 
n_source_treat <- raw2 %>% 
  filter(element == "N") %>%
  select("treatment", "Gly", "Lys", "Phe", "Ser")
matrix3<- as.matrix(n_source_treat[,2:5]) # response variables
manova3<-manova(matrix3~as.factor(treatment), data=n_source_treat)
summary(manova3) # p = 0.11

# TROPHIC AAs x fraction
n_trophic_frac <- raw2 %>% 
  filter(element == "N") %>%
  select("fraction", "Glx", "Asx", "Pro", "Ala", "Leu", "Ile", "Val")
matrix4<- as.matrix(n_trophic_treat[,2:5]) # response variables
manova4<-manova(matrix4~as.factor(fraction), data=n_trophic_frac)
summary(manova4) # p = 2.127e-12

# TROPHIC AAs x treatment 
n_trophic_treat <- raw2 %>% 
  filter(element == "N") %>%
  select("treatment", "Glx", "Asx", "Pro", "Ala", "Leu", "Ile", "Val")
matrix5<- as.matrix(n_trophic_treat[,2:5]) # response variables
manova5<-manova(matrix5~as.factor(treatment), data=n_trophic_treat)
summary(manova5) # p = 0.4147

# TROPHIC AAs x ecotype
n_trophic_ecotype <- raw2 %>% 
  filter(element == "N") %>%
  select("Apo_Sym", "Glx", "Asx", "Pro", "Ala", "Leu", "Ile", "Val")
matrix6<- as.matrix(n_trophic_ecotype[,2:5]) # response variables
manova6<-manova(matrix6~as.factor(Apo_Sym), data=n_trophic_ecotype)
summary(manova6) # p = 0.445


# ISOTOPES: SIBER analysis ---------------------------------------------------------

# Set up & parameters for SIBER analysis & Baysean models 

EAAs <- raw2 %>%
  filter(element == "C") %>% 
  mutate(new = if_else(sym.cm2 > 500000, "SYM", "APO")) %>%
  filter(new != "NA") %>%
  dplyr::select(full_id, treatment, fraction, new, Ile, Lys, Leu, Phe, Thr, Val) 

# Create lists of plotting arguments to be passed onwards to each of the three plotting functions.
#community.hulls.args <- list(col = 1, lty = 1, lwd = 1) # Only need if using group & community 
group.ellipses.args  <- list(n = 100, p.interval = 0.4, lty = 1, lwd = 6) #p.interval = % data included, n = points made to smooth ellipse  
group.hulls.args     <- list(lty = 2, lwd = 2, col = "darkgray") # hulls 

# (Just Another Gibbs Sampler): program used to analyze Bayesian hierarchical models using Markov Chain Monte Carlo (MCMC) simulation.
# parms hold the parameters defining how the sampling algorithm is to run
parms <- list() # create empty list 
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10 # thin the posterior by this many
parms$n.chains <- 2 # run this many chains

## PRIORS 

# a prior is the initial belief or assumption about a model's parameters before any data is observed
# priors can typically be left alone, they z-score the data and means will inherently be near 0
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# separate data 
data_fraction <- EAAs %>% filter(treatment=="control")
data_ecotype <- EAAs %>% filter(fraction == "Host")
data_ecotype_s <- EAAs %>% filter(fraction == "Sym")
data_treatment <- EAAs %>% filter(fraction == "Host")
data_treatment_s <- EAAs %>% filter(fraction == "Sym")

# Create the PCA 
C_pca_fraction <- prcomp(data_fraction[,5:10], scale. = TRUE, center = TRUE)
C_pca_ecotype <- prcomp(data_ecotype[,5:10], scale. = TRUE, center = TRUE)
C_pca_ecotype_s <- prcomp(data_ecotype_s[,5:10], scale. = TRUE, center = TRUE)
C_pca_treatment <- prcomp(data_treatment[,5:10], scale. = TRUE, center = TRUE)
C_pca_treatment_s <- prcomp(data_treatment_s[,5:10], scale. = TRUE, center = TRUE)

# save PCA data 
pca_data_fraction <-data.frame(C_pca_fraction$x, fraction=data_fraction$fraction)
pca_data_ecotype <-data.frame(C_pca_ecotype$x, new=data_ecotype$new)
pca_data_ecotype_s <-data.frame(C_pca_ecotype_s$x, new=data_ecotype_s$new)
pca_data_treatment <-data.frame(C_pca_treatment$x, treatment=data_treatment$treatment)
pca_data_treatment_s <-data.frame(C_pca_treatment_s$x, treatment=data_treatment_s$treatment)

# Save the PCA data in the correct format (SIBER needs specific headings to work)
pca_data_fraction$comm <- 'z' #assign community group 
pca_data_ecotype$comm <- 'z' 
pca_data_ecotype_s$comm <- 'z' 
pca_data_treatment$comm <- 'z' 
pca_data_treatment_s$comm <- 'z' 

### MAKE PCAs 

# Graph PCA of each comparison 
graph_pca_fraction <- autoplot(C_pca_fraction, data=data_fraction ,color = 'fraction', frame = TRUE, frame.type = 't', loadings=TRUE,loadings.label=TRUE) +
  # AESTHETICS
  theme_bw() + 
  scale_color_manual(values = c("Black","#1b9e29"))+
  scale_fill_manual(values = c("Black","#1b9e29"))

# Graph PCA of each comparison for host
graph_pca_ecotype <- autoplot(C_pca_ecotype, data=data_ecotype ,color = 'new', frame = TRUE, frame.type = 't', loadings=TRUE,loadings.label=TRUE) +
  # AESTHETICS
  theme_bw() + 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+
  scale_fill_manual(values = c("#e4d2ba", "#724a29"))

# Graph PCA of each comparison for sym
graph_pca_ecotype_s <- autoplot(C_pca_ecotype_s, data=data_ecotype_s ,color = 'new', frame = TRUE, frame.type = 't', loadings=TRUE,loadings.label=TRUE) +
  # AESTHETICS
  theme_bw() + 
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+
  scale_fill_manual(values = c("#e4d2ba", "#724a29"))

# Graph PCA of each comparison 
graph_pca_treatment <- autoplot(C_pca_treatment, data=data_treatment,color = 'treatment', frame = TRUE, frame.type = 't', loadings=TRUE,loadings.label=TRUE) +
  # AESTHETICS
  theme_bw() + 
  scale_color_manual(values = c("#C0C0E1", "#210124", "#8A4594"))+
  scale_fill_manual(values = c("#C0C0E1", "#210124", "#8A4594"))

# Graph PCA of each comparison 
graph_pca_treatment_s <- autoplot(C_pca_treatment_s, data=data_treatment_s,color = 'treatment', frame = TRUE, frame.type = 't', loadings=TRUE,loadings.label=TRUE) +
  # AESTHETICS
  theme_bw() + 
  scale_color_manual(values = c("#C0C0E1", "#210124", "#8A4594"))+
  scale_fill_manual(values = c("#C0C0E1", "#210124", "#8A4594"))

### Generate SIBER datasets 

sib_fraction <- pca_data_fraction %>% dplyr::select(PC1, PC2, fraction, comm) 
colnames(sib_fraction) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_ecotype <- pca_data_ecotype %>% dplyr::select(PC1, PC2, new, comm) 
colnames(sib_ecotype) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_ecotype_s <- pca_data_ecotype_s %>% dplyr::select(PC1, PC2, new, comm) 
colnames(sib_ecotype_s) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_treatment <- pca_data_treatment %>% dplyr::select(PC1, PC2, treatment, comm) 
colnames(sib_treatment) <- c("iso1", "iso2", "group", "community") # rename column labels

sib_treatment_s <- pca_data_treatment_s %>% dplyr::select(PC1, PC2, treatment, comm) 
colnames(sib_treatment_s) <- c("iso1", "iso2", "group", "community") # rename column labels

# Plot Ellipses and SEAs

# create siber object 
siber_fraction <- createSiberObject(sib_fraction)
siber_ecotype <- createSiberObject(sib_ecotype)
siber_ecotype_s <- createSiberObject(sib_ecotype_s)
siber_treatment <- createSiberObject(sib_treatment)
siber_treatment_s <- createSiberObject(sib_treatment_s)

# Siber SUMMARY STATS 

# Calculate summary statistics for each group: TA, SEA and SEAc (single metrics)
summary_fraction <- groupMetricsML(siber_fraction)
summary_ecotype <- groupMetricsML(siber_ecotype)
summary_ecotype_s <- groupMetricsML(siber_ecotype_s)
summary_treatment <- groupMetricsML(siber_treatment)
summary_treatment_s <- groupMetricsML(siber_treatment_s)

# Ellipse overlap

# the overlap betweeen the corresponding 40% prediction ellipses is given by:
# returns 3 numbers: size of ellipse 2, size of overlap 
overlap_fraction <- maxLikOverlap("z.Host", "z.Sym", siber_fraction, p.interval = 0.4, n = 100)
overlap_ecotype <- maxLikOverlap("z.APO", "z.SYM", siber_ecotype, p.interval = 0.4, n = 100)
overlap_ecotype_s <- maxLikOverlap("z.APO", "z.SYM", siber_ecotype_s, p.interval = 0.4, n = 100)

overlap_treatment_deep <- maxLikOverlap("z.control", "z.deep", siber_treatment, p.interval = 0.4, n = 100)
overlap_treatment_na <- maxLikOverlap("z.deep", "z.shade", siber_treatment, p.interval = 0.4, n = 100)
overlap_treatment_shade <- maxLikOverlap("z.control", "z.shade", siber_treatment, p.interval = 0.4, n = 100)

overlap_treatment_deep_s <- maxLikOverlap("z.control", "z.deep", siber_treatment_s, p.interval = 0.4, n = 100)
overlap_treatment_na_s <- maxLikOverlap("z.deep", "z.shade", siber_treatment_s, p.interval = 0.4, n = 100)
overlap_treatment_shade_s <- maxLikOverlap("z.control", "z.shade", siber_treatment_s, p.interval = 0.4, n = 100)

# calculate overlap as a proportion of the non-overlapping area of both ellipses 
percent_overlap_fraction <- overlap_fraction[3] / 
  (overlap_fraction[2] + overlap_fraction[1] - overlap_fraction[3])
percent_overlap_ecotype <- overlap_ecotype[3] / 
  (overlap_ecotype[2] + overlap_ecotype[1] - overlap_ecotype[3])
percent_overlap_ecotype_s <- overlap_ecotype_s[3] / 
  (overlap_ecotype_s[2] + overlap_ecotype_s[1] - overlap_ecotype_s[3])
percent_overlap_treatment_deep <- overlap_treatment_deep[3] / 
  (overlap_treatment_deep[2] + overlap_treatment_deep[1] - overlap_treatment_deep[3])
percent_overlap_treatment_shade <- overlap_treatment_shade[3] / 
  (overlap_treatment_shade[2] + overlap_treatment_shade[1] - overlap_treatment_shade[3])
percent_overlap_treatment_na <- overlap_treatment_na[3] / 
  (overlap_treatment_na[2] + overlap_treatment_na[1] - overlap_treatment_na[3])

percent_overlap_treatment_deep_s <- overlap_treatment_deep_s[3] / 
  (overlap_treatment_deep_s[2] + overlap_treatment_deep_s[1] - overlap_treatment_deep_s[3])
percent_overlap_treatment_shade_s <- overlap_treatment_shade_s[3] / 
  (overlap_treatment_shade_s[2] + overlap_treatment_shade_s[1] - overlap_treatment_shade_s[3])
percent_overlap_treatment_na_s <- overlap_treatment_na_s[3] / 
  (overlap_treatment_na_s[2] + overlap_treatment_na_s[1] - overlap_treatment_na_s[3])

# BAYESIAN MODEL 

model_fraction <- siberMVN(siber_fraction, parms, priors)
model_ecotype <- siberMVN(siber_ecotype, parms, priors)
model_ecotype_s <- siberMVN(siber_ecotype_s, parms, priors)
model_treatment <- siberMVN(siber_treatment, parms, priors)
model_treatment_s <- siberMVN(siber_treatment_s, parms, priors)

# SEAb
SEAb_fraction <- siberEllipses(model_fraction)
SEAb_ecotype <- siberEllipses(model_ecotype)
SEAb_ecotype_s <- siberEllipses(model_ecotype_s)
SEAb_treatment <- siberEllipses(model_treatment)
SEAb_treatment_s <- siberEllipses(model_treatment_s)

# Display the modes from the SEAb analysis 
SEAb_modes_fraction <- lapply(as.data.frame(SEAb_fraction), 
                              function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_ecotype  <- lapply(as.data.frame(SEAb_ecotype), 
                              function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_ecotype_s  <- lapply(as.data.frame(SEAb_ecotype_s), 
                                function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_treatment <- lapply(as.data.frame(SEAb_treatment), 
                               function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)
SEAb_modes_treatment_s <- lapply(as.data.frame(SEAb_treatment_s), 
                                 function(x,...){tmp<-hdrcde::hdr(x)$mode}, prob = cr.p, all.modes=T)

# display the credibility intervals (similar to confidence intervals)
SEAb_cred_fraction <- lapply(as.data.frame(SEAb_fraction), 
                             function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_ecotype <- lapply(as.data.frame(SEAb_ecotype), 
                            function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_ecotype_s <- lapply(as.data.frame(SEAb_ecotype_s), 
                              function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_treatment <- lapply(as.data.frame(SEAb_treatment), 
                              function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_treatment_s <- lapply(as.data.frame(SEAb_treatment_s), 
                                function(x,...){tmp<-hdrcde::hdr(x)$hdr},prob = cr.p)
SEAb_cred_fraction

# Centroid distances 

# Graph prep 

f_clrs <- matrix(c("#1b9e29","#1b9e29","#1b9e29", 
                   "#353535", "#353535", "#353535"), nrow = 3, ncol = 2)

e_clrs <- matrix(c("#e4d2ba", "#e4d2ba", "#e4d2ba", 
                   "#724a29", "#724a29", "#724a29"), nrow = 3, ncol = 2)


t_clrs <- matrix(c("#4d2f51", "#4d2f51", "#4d2f51",
                   "#C0C0E1", "#C0C0E1", "#C0C0E1",
                   "#8A4594", "#8A4594", "#8A4594"), nrow = 3, ncol = 3)


# ISOTOPES: SIBER plots (FIG 6) ---------------------------------------------------

# SIBER BY FRACTION

# stats
summary_fraction # TA, SEA, SEAc 
overlap_fraction # ellipse overlap fraction
percent_overlap_fraction # percent ellipse overlap 
SEAb_modes_fraction # modes
SEAb_cred_fraction # credibilities 

#plot siber 
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_fraction.png", width = 700, height = 700)
# Print siber plot 
palette(c("#353535","#1b9e29"))
plotSiberObject(siber_fraction,
                  ax.pad = 2,  #determines the padding applied around the extremes of the data.
                  ellipses = T, #Ellipses are drawn for each group independently with 
                  group.ellipses.args = group.ellipses.args, # ellipse data 
                  group.hulls = T, #makes hulls for each group (polygon around all ) 
                  group.hulls.args = group.hulls.args, # group hull data 
                  bty = "L", # axes lines
                  iso.order = c(1,2),
                  xlab = "",
                  ylab = "", 
                  cex.lab = 1.5,  # Axis label font size
                  cex.axis = 2,  # Axis tick font size
                  cex.main = 1.8, 
                  #cex = 10,  # change size of points 
                  points.order = (19), # change point type
                  x.limits = c(-3.5,3.5),
                  y.limits = c(-2,3.5)
  )
dev.off()

# SEAb 
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_fraction_seab.png", width = 400, height = 700)
# SEAb density plot 
siberDensityPlot(
  SEAb_fraction,
  clr = f_clrs,
  xticklabels = c("", ""),              
  xlab = c(""),
  ylab = c(""),
  bty = "L",
  las = 1,
  cex.axis = 2,  # Axis tick font size
  ylim = c(2, 11), 
  main = ""
)

# Step 2: Add points to the existing plot
points(1:ncol(SEAb_fraction), summary_fraction[3,], col = "red", pch = "x", lwd = 2)   
dev.off()

# SIBER BY ECOTYPE - HOST 

# stats
summary_ecotype # TA, SEA, SEAc 
overlap_ecotype# ellipse overlap fraction
percent_overlap_ecotype # percent ellipse overlap 
SEAb_modes_ecotype # modes
SEAb_cred_ecotype # credibilities 

# plot siber 
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_ecotype.png", width = 700, height = 700)
# Print siber plot 
palette(c("#e4d2ba", "#724a29"))
plotSiberObject(siber_ecotype,
                ax.pad = 2,  #determines the padding applied around the extremes of the data.
                ellipses = T, #Ellipses are drawn for each group independently with 
                group.ellipses.args = group.ellipses.args, # ellipse data 
                group.hulls = T, #makes hulls for each group (polygon around all ) 
                group.hulls.args = group.hulls.args, # group hull data 
                bty = "L", # axes lines
                iso.order = c(1,2),
                xlab = "",
                ylab = "", 
                cex = 0.11, # change size of points 
                cex.lab = 1.5,  # Axis label font size
                cex.axis = 2,  # Axis tick font size
                cex.main = 1.8,
                points.order = (19), # change point type
                x.limits = c(-4,5),
                y.limits = c(-2,4)
)
dev.off()

# plot SEAb
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_ecotype_seab.png", width = 400, height = 700)
# SEAb density plot 
siberDensityPlot(
  SEAb_ecotype,                                          
  clr = e_clrs, 
  xticklabels = c("", ""),            
  xlab = c(""),
  ylab = c(""),
  bty = "L",
  las = 1,
  cex.axis = 2,  # Axis tick font size
  ylim = c(2, 14), 
  main = ""
)

# Step 2: Add points to the existing plot
points(1:ncol(SEAb_ecotype), summary_ecotype[3,], col = "red", pch = "x", lwd = 2)     
dev.off()

# SIBER BY ECOTYPE - HOST 

# stats
summary_ecotype_s # TA, SEA, SEAc 
overlap_ecotype_s # ellipse overlap fraction
percent_overlap_ecotype_s # percent ellipse overlap 
SEAb_modes_ecotype_s # modes
SEAb_cred_ecotype_s # credibilities 

# plot SIBER 
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_ecotype_s.png", width = 700, height = 700)
# Print siber plot 
palette(c("#e4d2ba", "#724a29"))
plotSiberObject(siber_ecotype_s,
                ax.pad = 2,  #determines the padding applied around the extremes of the data.
                ellipses = T, #Ellipses are drawn for each group independently with 
                group.ellipses.args = group.ellipses.args, # ellipse data 
                group.hulls = T, #makes hulls for each group (polygon around all ) 
                group.hulls.args = group.hulls.args, # group hull data 
                bty = "L", # axes lines
                iso.order = c(1,2),
                xlab = "",
                ylab = "", 
                cex = 0.11, # change size of points 
                cex.lab = 1.5,  # Axis label font size
                cex.axis = 2,  # Axis tick font size
                cex.main = 1.8,
                points.order = (19), # change point type
                x.limits = c(-3.5,3.5),
                y.limits = c(-3,3)
)
dev.off()

# plot SEAb
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_ecotype_s_seab.png", width = 400, height = 700)
# SEAb density plot 
siberDensityPlot(
  SEAb_ecotype_s,  
  clr = e_clrs, 
  xticklabels = c("", ""),              
  xlab = c(""),
  ylab = c(""),
  bty = "L",
  las = 1,
  cex.axis = 2,  # Axis tick font size
  ylim = c(0, 11), 
  main = ""
)
# Step 2: Add points to the existing plot
points(1:ncol(SEAb_ecotype_s), summary_ecotype_s[3,], col = "red", pch = "x", lwd = 2)     
dev.off()

# SIBER BY TREATMENT - HOST 

# stats
summary_treatment # TA, SEA, SEAc 
  # control vs. deep 
overlap_treatment_deep # ellipse overlap fraction
percent_overlap_treatment_deep # percent ellipse overlap 
  # control vs. shade
overlap_treatment_shade # ellipse overlap fraction
percent_overlap_treatment_shade # percent ellipse overlap 
SEAb_modes_treatment # modes
SEAb_cred_treatment # credibilities 

# plot SIBER 
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_treat_host.png", width = 700, height = 700)
# Print siber plot 
palette(c("#C0C0E1", "#8A4594", "#210124"))
plotSiberObject(siber_treatment,
                ax.pad = 2,  #determines the padding applied around the extremes of the data.
                ellipses = T, #Ellipses are drawn for each group independently with 
                group.ellipses.args = group.ellipses.args, # ellipse data 
                group.hulls = T, #makes hulls for each group (polygon around all ) 
                group.hulls.args = group.hulls.args, # group hull data 
                bty = "L", # axes lines
                iso.order = c(1,2),
                xlab = "",
                ylab = "", 
                cex = 0.11, # change size of points 
                cex.lab = 1.5,  # Axis label font size
                cex.axis = 2,  # Axis tick font size
                cex.main = 1.8,
                points.order = (19), # change point type
                x.limits = c(-4,5),
                y.limits = c(-2,4)
)
dev.off()

# Plot SEAb
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_treat_host_seab.png", width = 400, height = 700)
# SEAb density plot 
siberDensityPlot(
  SEAb_treatment,                                          
  clr = t_clrs, 
  xticklabels = c("", "", ""),              
  xlab = c(""),
  ylab = c("" ),
  bty = "L",
  las = 1,
  cex.axis = 2,  # Axis tick font size
  ylim = c(2, 13), 
  main = ""
)
# Step 2: Add points to the existing plot
points(1:ncol(SEAb_treatment), summary_treatment[3,], col = "red", pch = "x", lwd = 2)     
dev.off()


# SIBER BY TREATMENT - SYMBIONTS

# stats
summary_treatment_s # TA, SEA, SEAc 
# control vs. deep 
overlap_treatment_deep_s # ellipse overlap fraction
percent_overlap_treatment_deep_s # percent ellipse overlap 
# control vs. shade
overlap_treatment_shade_s # ellipse overlap fraction
percent_overlap_treatment_shade_s # percent ellipse overlap 
SEAb_modes_treatment_s # modes
SEAb_cred_treatment_s # credibilities 

# plot SIBER 
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_treat_sym.png", width = 700, height = 700)
# Print siber plot 
palette(c( "#C0C0E1", "#8A4594", "#210124"))
plotSiberObject(siber_treatment_s,
                ax.pad = 2,  #determines the padding applied around the extremes of the data.
                ellipses = T, #Ellipses are drawn for each group independently with 
                group.ellipses.args = group.ellipses.args, # ellipse data 
                group.hulls = T, #makes hulls for each group (polygon around all ) 
                group.hulls.args = group.hulls.args, # group hull data 
                bty = "L", # axes lines
                iso.order = c(1,2),
                xlab = "",
                ylab = "", 
                cex = 0.11, # change size of points 
                cex.lab = 1.5,  # Axis label font size
                cex.axis = 2,  # Axis tick font size
                cex.main = 1.8, 
                points.order = (19), # change point type
                x.limits = c(-3,3.5),
                y.limits = c(-3,3)
)

dev.off()

# Plot SEAb
png("~/Desktop/GITHUB/Pub-Astrangia-trophic-ecology/FIGURES/SIBER/siber_plot_treat_sym_seab.png", width = 400, height = 700)
# SEAb density plot 
siberDensityPlot(
  SEAb_treatment_s,                                          
  clr = t_clrs, 
  xticklabels = c("", "", ""),            
  xlab = c(""),
  ylab = expression("" ),
  bty = "L",
  las = 1,
  ylim = c(0, 15), 
  cex.axis = 2,  # Axis tick font size
  main = ""
)

# Step 2: Add points to the existing plot
points(1:ncol(SEAb_treatment_s), summary_treatment_s[3,], col = "red", pch = "x", lwd = 2)  
dev.off()

# PLOT LEGEND 

legend_data <- data.frame(
  category = factor(
    c("Host", "Symbiont", "Control", "Shade", "Deep", "Apo", "Sym"),
    levels = c("Host", "Symbiont", "Control", "Shade", "Deep", "Apo", "Sym")  # Ensure correct order
  ),
  group = factor(
    c("Fraction", "Fraction", "Treatment", "Treatment", "Treatment", "Final Ecotype", "Final Ecotype"),
    levels = c("Fraction", "Treatment", "Final Ecotype")  # Order of the groups
  ),
  value = c(1, 1, 1, 1, 1, 1, 1)  # Dummy values for plotting
)

# Create the ggplot
legend_plot <- ggplot() +
  # First group with its own scale
  geom_bar(
    data = subset(legend_data, group == "Fraction"),
    aes(x = group, y = value, fill = category),
    stat = "identity",
    position = "dodge"
  ) +
  scale_fill_manual(
    name = "Fraction",
    values = c("Host" = "#353535", "Symbiont" = "#1b9e29"),
    breaks = c("Host", "Symbiont")  # Explicit order of categories in the legend
  ) +
  new_scale_fill() +  # Allow new scale for Group 2
  
  # Second group with its own scale
  geom_bar(
    data = subset(legend_data, group == "Treatment"),
    aes(x = group, y = value, fill = category),
    stat = "identity",
    position = "dodge"
  ) +
  scale_fill_manual(
    name = "Treatment",
    values = c("Control" = "#C0C0E1", "Shade" = "#8A4594", "Deep" = "#210124"),
    breaks = c("Control", "Shade", "Deep")  # Explicit order of categories in the legend
  ) +
  new_scale_fill() +  # Allow new scale for Group 3
  
  # Third group with its own scale
  geom_bar(
    data = subset(legend_data, group == "Final Ecotype"),
    aes(x = group, y = value, fill = category),
    stat = "identity",
    position = "dodge"
  ) +
  scale_fill_manual(
    name = "Final Ecotype",
    values = c("Apo" = "#e4d2ba", "Sym" = "#724a29"),
    breaks = c("Apo", "Sym")  # Explicit order of categories in the legend
  ) +
  theme_void() +  # Remove the axes and background
  theme(
    legend.position = "right",  # Position all legends
    legend.title = element_text(face = "bold", size = 15),  # Bold legend titles
    legend.text = element_text(size = 15),  # Adjust legend text size
    legend.box = "horizontal"  # Ensure that the legends are in separate boxes
  ) +
  guides(
    fill = guide_legend(ncol = 1)  # Arrange legends in a single column
  )

# Display the legend-only plot
legend_plot

ggsave("TLAP_CSIA_C_legend.jpg", plot=legend_plot, path = "FIGURES/SIBER/",width = 14, height = 5)


# ISOTOPES: N source mean   -----------------------------------------------------------------


##### Figure out mean values of source AAs 

NB_vals <- raw2 %>% 
  filter(element == "N")%>% 
  dplyr::select(sample_id, Phe, Lys, Arg, Gly, Ser) %>%
  pivot_longer(cols = c(Phe, Lys, Arg, Gly, Ser), values_to = "deltN", names_to = "AA") %>%
  filter(deltN != "NA") %>%
  group_by(AA) %>%
  summarise(mean = mean(deltN), sd = sd(deltN))
# TROPHIC POSITION --------------------------------------------------------

# Trophic position (tp) can be calculated using the following equation:  tp = 1+((Glx-Phe-beta)/TDF)
# where 
  # Glx = the delt15N of Glutamine + Glutamic Acid (Trophic AA)
  # Phe = the delt15N of Phenylalanine (Source AA)
  # TDF = Trophic discrimination factor = 7.6 (because coral is an amonia producer)
  # beta =  Glx-phe of baseline primary producers. = 3.3 (Phytoplankton baseline)

# isoloate data we need to calculate TP 
tp <- raw2 %>% 
  filter(element == "N") %>% 
  dplyr::select(sample_id, Apo_Sym, treatment, fraction, sym.cm2, Glx, Phe) %>%
  mutate(glxphe = Glx-Phe) %>%
  mutate(trophic_position = 1+((glxphe-3.3)/7.6)) %>%
  mutate(new = if_else(sym.cm2 > 500000, "SYM", "APO"))

tp2 <- tp %>%
  filter(sample_id != "G5-C1-49")

tp_host <- tp2 %>% filter(fraction=="Host")
tp_sym <- tp2 %>% filter(fraction=="Sym")
tp_control <- tp2 %>% filter(treatment == "control")

# mean TP
tp_means <-tp2 %>%
  group_by(fraction) %>%
  summarise(mean = mean(trophic_position), SD = sd(trophic_position)) 
tp_means

tp2_control <- tp2 %>% filter(treatment == "control")

# TP by fraction 
tp_frac <- ggplot(tp2_control, aes(x=fraction, y=trophic_position, fill=fraction)) +
  #DATA 
  geom_boxplot(alpha=0.8, outlier.size=0) + 
  geom_point(position = pj, size=4, show.legend = FALSE)+
  #AESTHETICS 
  theme_bw() + 
  labs(y= "", x = "Fraction", fill="") +
  scale_fill_manual(values = c("#333333","#1b9e29")) +
  theme(
    legend.position = "none",
    axis.text.y=element_text(size=37), 
    axis.title=element_text(size=40), 
    axis.text.x=element_text(size=35),
    plot.margin=unit(c(0,1,0.5,0),"cm")
  )+
  stat_compare_means(comparisons = treatment_comparisons3, method = "wilcox.test", size=12) +
  scale_x_discrete(labels=c("Host" = "Host", "Sym" = "Symbiont")) 
#ylim(1.95,2.9)

### HOST 

# TP by treatment 
tp_treat <- ggplot(tp_host, aes(x=treatment, y=trophic_position, fill=treatment)) +
  #DATA 
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=4, show.legend = FALSE)+
  #AESTHETICS 
  theme_bw() + 
  labs(y= "", x = "", fill="") +
  scale_fill_manual(values = c("#C0C0E1", "#8A4594", "#210124")) +
  theme(
    legend.position = "none", 
    axis.text.y=element_text(size=37), 
    axis.title=element_text(size=40), 
    axis.text.x=element_blank(),
    plot.margin=unit(c(0,0,0,0),"cm")
  )+
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=12) +
  scale_x_discrete(labels=c("APO" = "Aposymbiotic", "SYM" = "Symbiotic")) +
  ylim(1.95, 3)

tp_symbionts <- ggplot(tp_host, aes(x=sym.cm2, y=trophic_position)) +
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_point(aes(color=Apo_Sym), size=5) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  geom_vline(xintercept = 500000, linetype="dashed") + 
  #AESTHETICS 
  theme_bw()+
  labs(y= "", x = "", color="Ecotype") +
  theme( 
    legend.position = "none", 
    axis.text.y=element_blank(),
    axis.title=element_text(size=40), 
    axis.text.x=element_blank(),
    plot.margin=unit(c(0,1,0,0),"cm")
  )+
  ylim(1.95, 3) +
  xlim(-1000, 1530000)

#### SYM 

# TP by fraction 
tp_treat_s <- ggplot(tp_sym, aes(x=treatment, y=trophic_position, fill=treatment)) +
  #DATA 
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=4, show.legend = FALSE)+
  #AESTHETICS 
  theme_bw() + 
  labs(y= "", x = "Treatment", fill="") +
  scale_fill_manual(values = c("#C0C0E1", "#8A4594", "#210124")) +
  theme( 
    legend.position = "none", 
    axis.text.y=element_text(size=37), 
    axis.title=element_text(size=40), 
    axis.text.x=element_text(size=35),
    plot.margin=unit(c(0,0,0.5,0),"cm")
  ) +
  scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=12) + 
  ylim(1.9, 2.9)

tp_sym_s <- ggplot(tp_sym, aes(x=sym.cm2, y=trophic_position)) +
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_point(aes(color=Apo_Sym), size=5) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  geom_vline(xintercept = 500000, linetype="dashed") + 
  #AESTHETICS 
  theme_bw()+
  labs(y= "", x = "Symbiont Density" ~(cells/cm^2), color="Ecotype") +
  theme(
    legend.position = "none", 
    axis.text.y= element_blank(),
    axis.title=element_text(size=36), 
    axis.text.x=element_text(size=34),
    plot.margin=unit(c(0,1,0,0),"cm")
  ) +
  ylim(1.9, 2.9) +
  xlim(-1000, 1530000)

## Arrange host graphs 

# designate the y axis label
yleft <- textGrob("Trophic Position (glx-phe)", rot=90, gp = gpar(fontsize = 55))
yleft_h <- textGrob("Host", rot=90, gp = gpar(fontsize = 45, fontface = "bold"))
yleft_s <- textGrob("Symbiont", rot=90, gp = gpar(fontsize = 45, fontface = "bold"))

## Arrange host graphs 
tp_arrange <- grid.arrange( tp_treat, tp_symbionts, nrow=1, left=yleft_h, widths=c(1,1.5))

## Arrange symbiont graphs 
tp_arrange_s <- grid.arrange(tp_treat_s, tp_sym_s, nrow=1, left=yleft_s, widths=c(1,1.5))

tp_arrange_both <- grid.arrange(tp_arrange, tp_arrange_s, nrow=2)

tp_arrange_all <- grid.arrange(tp_frac, tp_arrange_both, nrow=1, left=yleft, widths=c(0.5,1.5),
                               bottom = grobTree(
                                 textGrob("A", x = 0.06, y = 245, just = "left", gp = gpar(fontsize = 42, fontface = "bold")),
                                 textGrob("B", x = 0.35, y = 245, just = "right", gp = gpar(fontsize = 42, fontface = "bold")),
                                 textGrob("C", x = 0.605, y = 245, just = "right", gp = gpar(fontsize = 42, fontface = "bold")),
                                 textGrob("D", x = 0.35, y = 120, just = "right", gp = gpar(fontsize = 42, fontface = "bold")),
                                 textGrob("E", x = 0.605, y = 120, just = "right", gp = gpar(fontsize = 42, fontface = "bold"))
                               ))
tp_arrange_all

ggsave("TLAP_FIG_7_trophic_position.jpg", plot=tp_arrange_all, path = "FIGURES/", width = 30, height = 25)

# SUMV --------------------------------------------------------------------

# calculate sum v
sumV <- raw2 %>% 
  filter(element == "N") %>% 
  select(sample_id, Apo_Sym, treatment, fraction, sym.cm2, Glx,Asx,Pro,Ala,Leu,Ile,Val,Phe,Met,Lys,Arg,Gly,Ser,Thr) %>%
  mutate(mean_trophic = (Ala + Glx + Leu + Ile + Pro + Asx)/6) %>%
  mutate(sum_v = 1/6*(abs(Ala-mean_trophic)+abs(Glx-mean_trophic)+abs(Leu-mean_trophic)+abs(Ile-mean_trophic)+abs(Pro-mean_trophic)+abs(Asx-mean_trophic))) %>%
  mutate(new = if_else(sym.cm2 > 500000, "SYM", "APO")) %>%
  dplyr::select(sample_id, Apo_Sym, treatment, fraction, sym.cm2 ,new, sum_v)

# remove outliers 
sumV <- sumV %>% 
  filter(!is.na(sum_v)) %>%
  mutate(zscore = (sum_v - mean(sum_v))/sd(sum_v)) %>%
  filter(abs(zscore)<3) %>%
  dplyr::select(sample_id, Apo_Sym, treatment, fraction, sym.cm2 ,new, sum_v)

# sum v host only 
sumV_host <- sumV %>%
  filter(fraction == "Host")

# sum v host only 
sumV_sym <- sumV %>%
  filter(fraction == "Sym")

## GRAPH 
sumV_control <- sumV %>% filter(treatment == "control")

# TP by fraction 
sumV_frac <- ggplot(sumV_control, aes(x=fraction, y=sum_v, fill=fraction)) +
  #DATA 
  geom_boxplot(alpha=0.8, outlier.size=0) + 
  geom_point(position = pj, size=4, show.legend = FALSE)+
  #AESTHETICS 
  theme_bw() + 
  labs(y= "", x = "Fraction", fill="") +
  scale_fill_manual(values = c("#333333","#1b9e29")) +
  theme(
    legend.position = "none",
    axis.text.y=element_text(size=37), 
    axis.title=element_text(size=40), 
    axis.text.x=element_text(size=35),
    plot.margin=unit(c(0,1,0.5,0),"cm")
  )+
  stat_compare_means(comparisons = treatment_comparisons3, method = "wilcox.test", size=12) +
  scale_x_discrete(labels=c("Host" = "Host", "Sym" = "Symbiont")) 
#ylim(1.95,2.9)

### HOST 

# TP by treatment 
sumV_treat <- ggplot(sumV_host, aes(x=treatment, y=sum_v, fill=treatment)) +
  #DATA 
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=4, show.legend = FALSE)+
  #AESTHETICS 
  theme_bw() + 
  labs(y= "", x = "", fill="") +
  scale_fill_manual(values = c("#C0C0E1", "#8A4594", "#210124")) +
  theme(
    legend.position = "none", 
    axis.text.y=element_text(size=37), 
    axis.title=element_text(size=40), 
    axis.text.x=element_blank(),
    plot.margin=unit(c(0,0,0,0),"cm")
  )+
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=12) +
  scale_x_discrete(labels=c("APO" = "Aposymbiotic", "SYM" = "Symbiotic")) +
  ylim(1.1, 3.3)

sumV_symbionts <- ggplot(sumV_host, aes(x=sym.cm2, y=sum_v)) +
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_point(aes(color=Apo_Sym), size=5) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  geom_vline(xintercept = 500000, linetype="dashed") + 
  #AESTHETICS 
  theme_bw()+
  labs(y= "", x = "", color="Ecotype") +
  theme( 
    legend.position = "none", 
    axis.text.y=element_blank(),
    axis.title=element_text(size=40), 
    axis.text.x=element_blank(),
    plot.margin=unit(c(0,1,0,0),"cm")
  ) +
  ylim(1.1, 3.3)+
  xlim(-1000, 1100000)

#### SYM 

# TP by fraction 
sumV_treat_s <- ggplot(sumV_sym, aes(x=treatment, y=sum_v, fill=treatment)) +
  #DATA 
  geom_boxplot(alpha=0.9, outlier.size=0) + 
  geom_point(position = pj, size=4, show.legend = FALSE)+
  #AESTHETICS 
  theme_bw() + 
  labs(y= "", x = "Treatment", fill="") +
  scale_fill_manual(values = c("#C0C0E1", "#8A4594", "#210124")) +
  theme( 
    legend.position = "none", 
    axis.text.y=element_text(size=37), 
    axis.title=element_text(size=40), 
    axis.text.x=element_text(size=35),
    plot.margin=unit(c(0,0,0.5,0),"cm")
  ) +
  scale_x_discrete(labels=c('Control','Shade', 'Deep') ) +
  stat_compare_means(comparisons = treatment_comparisons_combined, method = "wilcox.test", size=12) +
  #scale_x_discrete(labels=c("APO" = "Aposymbiotic", "SYM" = "Symbiotic")) +
  scale_y_continuous(labels = number_format(accuracy = 0.1), limits = c(1.1, 4.0)) 

sumV_sym_s <- ggplot(sumV_sym, aes(x=sym.cm2, y=sum_v)) +
  #DATA 
  geom_smooth(method="lm", formula = y ~ x, color="black", fill="gray") +
  geom_point(aes(color=Apo_Sym), size=5) + 
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           size=8, label.x = 600000) +
  scale_color_manual(values = c("#e4d2ba", "#724a29"))+ 
  geom_vline(xintercept = 500000, linetype="dashed") + 
  #AESTHETICS 
  theme_bw()+
  labs(y= "", x = "Symbiont Density" ~(cells/cm^2), color="Ecotype") +
  theme(
    legend.position = "none", 
    axis.text.y= element_blank(),
    axis.title=element_text(size=36), 
    axis.text.x=element_text(size=34),
    plot.margin=unit(c(0,1,0,0),"cm")
  )+
  ylim(1.1, 4.0) +
  xlim(-1000, 1100000)

## Arrange host graphs 

# designate the y axis label
yleft <- textGrob("ΣV", rot=90, gp = gpar(fontsize = 55, fontface = "bold"))
yleft_h <- textGrob("Host", rot=90, gp = gpar(fontsize = 45, fontface = "bold"))
yleft_s <- textGrob("Symbiont", rot=90, gp = gpar(fontsize = 45, fontface = "bold"))

## Arrange host graphs 
sumV_arrange <- grid.arrange(sumV_treat, sumV_symbionts, nrow=1, left=yleft_h, widths=c(1,1.5))

## Arrange symbiont graphs 
sumV_arrange_s <- grid.arrange(sumV_treat_s, sumV_sym_s, nrow=1, left=yleft_s, widths=c(1,1.5))

sumV_arrange_both <- grid.arrange(sumV_arrange, sumV_arrange_s, nrow=2)

sumV_arrange_all <- grid.arrange(sumV_frac, sumV_arrange_both, nrow=1, left=yleft, widths=c(0.5,1.5),
                                 bottom = grobTree(
                                   textGrob("A", x = 0.04, y = 245, just = "left", gp = gpar(fontsize = 42, fontface = "bold")),
                                   textGrob("B", x = 0.34, y = 245, just = "right", gp = gpar(fontsize = 42, fontface = "bold")),
                                   textGrob("C", x = 0.6, y = 245, just = "right", gp = gpar(fontsize = 42, fontface = "bold")),
                                   textGrob("D", x = 0.34, y = 120, just = "right", gp = gpar(fontsize = 42, fontface = "bold")),
                                   textGrob("E", x = 0.6, y = 120, just = "right", gp = gpar(fontsize = 42, fontface = "bold"))
                                 ))
sumV_arrange_all

ggsave("TLAP_FIG_S2_sumV.jpg", plot=sumV_arrange_all, path = "FIGURES/", width = 30, height = 25)

# PHYSIOLOGY + ISOTOPES: statistics --------------------------------------------------

# Add sum v and tp to the big data set
tp_prep <- tp2 %>%
  filter(sample_id != "G5-C1-49") %>%
  mutate(
    tp_s = ifelse(fraction == "Sym", trophic_position, NA),
    tp_h = ifelse(fraction == "Host", trophic_position, NA)) %>%
  select("sample_id", "Apo_Sym", "fraction", "treatment", "new","sym.cm2", "trophic_position", "tp_s", "tp_h")

sumV_prep <- sumV %>%
  mutate(
    sumV_s = ifelse(fraction == "Sym", sum_v, NA),
    sumV_h = ifelse(fraction == "Host", sum_v, NA) )%>%
  select("sample_id", "Apo_Sym", "fraction", "treatment", "new", "sym.cm2", "sum_v", "sumV_s", "sumV_h")

tp_sumv <- merge(tp_prep, sumV_prep, all = TRUE)

master2 <- merge(master, tp_sumv,all = TRUE)

# WILCOXON'S TESTS

# prep sepecific ones 
master_shade <- master2 %>% filter(treatment != "deep")
master_deep <- master2 %>% filter(treatment != "shade")
columns_to_test <- c("sym.cm2", "chla.ug.cm2", "log_chl_sym", "calice_density", "cre.umol.mgprot", "prot_mg.cm2", "AFDW.mg.cm2", "lipids.mg.cm2", "tp_s", "tp_h", "sumV_h", "sumV_s")  

apply_wilcoxon <- function(data, columns, group_col) {
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each column in the specified columns
  for (col in columns) {
    # Perform Wilcoxon test for each column
    test_result <- wilcox.test(data[[col]] ~ data[[group_col]], data = data, p.adjust.method = "none")
    
    # Create a summary for the test result
    summary <- data.frame(
      metric = col,      
      w = test_result$statistic, 
      p = test_result$p.value)
    
    # Append the result to the results list
    results_list[[col]] <- summary}
  
  # Combine the results into a single data frame
  final_results <- do.call(rbind, results_list)
  
  return(final_results)
}

# apply to control vs. shade then control vs. deep 
shade_results <- apply_wilcoxon(master_shade, columns_to_test, "treatment")
deep_results <- apply_wilcoxon(master_deep, columns_to_test, "treatment")
ecotype_results <- apply_wilcoxon(master2, columns_to_test, "Apo_Sym")

# merge 
order_index <- match(shade_results$metric, shade_results$metric)
wilcox_2 <- merge(shade_results , deep_results, by = "metric", suffixes = c("_shade", "_deep"), all = TRUE)
wilcox_3 <- merge(wilcox_2, ecotype_results, by ="metric", all = TRUE)

# Linear regressions 
run_linear_regressions <- function(data, columns_to_test) {
  # Initialize an empty list to store results
  results_list <- list()
  
  # Loop through each column in columns_to_test
  for (col in columns_to_test) {
    # Run the linear regression with 'sym.cm2' as the response variable
    model <- lm(sym.cm2 ~ get(col), data = data)
    
    # Extract R-squared and p-value from the model summary
    summary_model <- summary(model)
    r_squared <- summary_model$r.squared
    p_value <- summary_model$coefficients[2,4]   # p-value for the independent variable
    
    # Store the results in a dataframe
    results <- data.frame(
      metric = col,
      R2 = r_squared,
      p = p_value
    )
    
    # Append the results to the results list
    results_list[[col]] <- results
  }
  
  # Combine the results into a single dataframe
  final_results <- do.call(rbind, results_list)
  
  return(final_results)
}

# Example usage:
#columns_to_test <- c("column1", "column2", "column3")  # Replace with your actual column names
lm_results <- run_linear_regressions(master2, columns_to_test)

phys_stats <- merge(wilcox_3, lm_results, by = "metric")

# save 
write.csv(phys_stats, "STATS/TLAP_STATS_phys.csv", row.names = FALSE)

# CALCULATE PHYSIOLOGY MEANS 

master_long <- master2 %>%
  pivot_longer(cols = c("sym.cm2", "chla.ug.cm2", "chl_ug.sym", "calice_density", "cre.umol.mgprot", "prot_mg.cm2", "AFDW.mg.cm2", "lipids.mg.cm2","tp_s", "tp_h", "sumV_h", "sumV_s"), names_to = "metric", values_to = "value") %>%
  filter(value != "NA") %>%
  select(sample_id, treatment, Apo_Sym, metric, value)

phys_means <- master_long %>%
  group_by(metric, treatment) %>%
  summarise(mean = signif(mean(value),3), SD = signif(sd(value),3)) %>%
  pivot_wider(names_from = treatment, values_from = c(mean, SD)) %>%
  select(metric, mean_control, SD_control, mean_shade, SD_shade, mean_deep, SD_deep)

ecotype_means <-master_long %>%
  group_by(metric, Apo_Sym) %>%
  summarise(mean = signif(mean(value),3), SD = signif(sd(value),3)) %>%
  pivot_wider(names_from = Apo_Sym, values_from = c(mean, SD)) %>%
  select(metric, mean_APO, SD_APO, mean_SYM, SD_SYM)

# merge 
physiology_means <- merge(phys_means, ecotype_means, all.x = TRUE)

write.csv(physiology_means, "STATS/TLAP_STATS_phys_means.csv", row.names = FALSE)

### STATS FOR SYM VS HOST 

columns_to_test2 <- c("trophic_position", "sum_v")  

master3 <- master2 %>% filter(treatment == "control")

# apply to control vs. shade then control vs. deep 
fraction_results <- apply_wilcoxon(master3, columns_to_test2, "fraction")

# find means 
master_long2 <- master3 %>%
  pivot_longer(cols = c("trophic_position", "sum_v"), names_to = "metric", values_to = "value") %>%
  filter(value != "NA") %>%
  select(sample_id, fraction, metric, value)

phys_means2 <- master_long2 %>%
  group_by(fraction, metric) %>%
  summarise(mean = signif(mean(value),3), SD = signif(sd(value),3)) %>%
  pivot_wider(names_from = fraction, values_from = c(mean, SD)) %>%
  select(metric, mean_Host, SD_Host, mean_Sym, SD_Sym)

sumv_tp_means <- merge(fraction_results, phys_means2, all = TRUE)

write.csv(sumv_tp_means, "STATS/TLAP_STATS_sumv_tp_fractions.csv", row.names = FALSE)

# FIDELITY ----------------------------------------------------------------

# Calculate physiological fidelity to original ecotype 

fid <- master %>%
  dplyr::select(sample_id, Apo_Sym, treatment, full_treatment, new, sym.cm2, chla.ug.cm2, chl_ug.sym, calice_density, cre.umol.mgprot,
                prot_mg.cm2,AFDW.mg.cm2, lipids.mg.cm2) %>%
  pivot_longer(cols=c(sym.cm2, chla.ug.cm2, chl_ug.sym, calice_density, cre.umol.mgprot, prot_mg.cm2,
                      AFDW.mg.cm2, lipids.mg.cm2), values_to = "vals", names_to = "metric") %>%
  filter(vals != "NA") %>%
  dplyr::select(sample_id, Apo_Sym, treatment, full_treatment, new, metric, vals) %>%
  pivot_wider(names_from = treatment, values_from = vals)
 
t_test_p <- function(x, y) {
  result <- tryCatch({
    test <- wilcox.test(x, y)  # Perform Wilcoxon test
    data.frame(p.value = test$p.value, W = as.numeric(test$statistic))  # Return as data frame
  }, error = function(e) {
    # If an error occurs, return NA for both columns
    data.frame(p.value = NA, W = NA)
  })
  return(result)
}



fidelity_initial <- fid %>%
  group_by(metric) %>%
  summarise(
    initial_control = t_test_p(.data$control[Apo_Sym == "APO"], .data$control[Apo_Sym == "SYM"]),
    initial_shade = t_test_p(.data$shade[Apo_Sym == "APO"], .data$shade[Apo_Sym == "SYM"]),
    initial_deep = t_test_p (.data$deep[Apo_Sym == "APO"], .data$deep[Apo_Sym == "SYM"])
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4))) 

write.csv(fidelity_initial, file = "STATS/TLAP_STATS_fidelity.csv")

fidelity_final <- fid %>%
  group_by(metric) %>%
  summarise(
    final_shade = t_test_p(.data$shade[.data$new == "APO"],.data$shade[.data$new == "SYM"]),
    final_control = t_test_p(.data$control[.data$new == "APO"],.data$control[.data$new == "SYM"] ),
    final_deep = t_test_p(.data$deep[.data$new == "APO"],.data$deep[.data$new == "SYM"]) ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) 

