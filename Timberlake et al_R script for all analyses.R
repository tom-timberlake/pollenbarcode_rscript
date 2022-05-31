#######################################################################################################################################################
################################## TIMBERLAKE ET AL. - R SCRIPT FOR ALL ANALYSES AND FIGURE PLOTTING  #################################################
#######################################################################################################################################################



#######################################################################################################################
############# Floral resource analysis - modeling smooth, non-linear trend in floral abundance over time  #############
#######################################################################################################################

floral_abun_data <- read.csv(file.choose()) 

#### Running the GAM model 
library(mgcv)

GAM_model <- gam(floral_units + 0.0001 ~s(day, fx = T, k= 8), #k = knots - this value was chosen through careful model validation and inspection of AIC values
                     family = Gamma(link = "log"),
                     data = floral_abun_data)


##Model checking and validation
summary(GAM_model) 
gam.check(GAM_model)
AIC(GAM_model)

plot(jitter(floral_units, amount = 0.00001) ~ day, # plotting sugar versus week with a minor jitter adjustment
     data = floral_abun_data,
     ylab = "floral_units",    #labels y axis      
     xlab = "date", # labels x axis
     cex.lab = 1.8, # size of y and x axis labels
     cex.axis = 1.5,
     font.lab=2, #axis labels in bold
     las = 1,  # Ensure that all the axis labels are horizontal
     col = "black", 
     pch = 19, # changes the type of dot for data
     bty = "L") # changes the box type to an L rather than square


##Using pred function to plot model predictions and generate estimates of floral abundance on each day of the year

pdat <- expand.grid(day = seq(0,260,1)) #predicts values for hedge from the model from 0 to 32 in increments of exactly one day
pred <- predict (GAM_model, newdata = pdat, na.rm = T, type= "response", se.fit = TRUE)
predframe <- data.frame (pdat,level=0, preds = pred$fit, se = pred$se.fit)
lines(predframe$preds~predframe$day, lwd=3, col="black")

write.csv(pred, file="Chosen_file_name.csv") #write daily predictions of floral units to csv file



###########################################################################################################################
############# Obj. 1: Generating null networks and testing fit of each null model to the observed network #################
###########################################################################################################################

#install.packages("econullnetr")
library("econullnetr")

##Set working directory
setwd("C:/Users/thoma/OneDrive - University of Bristol/PhD/Ch.4_Pollen metabarcoding/EcoNullNetR results/Temporary results")

####Loading in data frames
#Consumer data
pollinators_early_pooled <- read.csv(file.choose()) 

#Resource data
plants_early_FU <- read.csv(file.choose()) #FU data subdivided by study site
plants_early_nectar <- read.csv(file.choose()) #nectar data subdivided by study site
plants_early_pollen <- read.csv(file.choose()) #pollen data subdivided by study site

#Plotting order data
plotting_order_early <- read.csv(file.choose()) # choose: 'Plotting order_Early_P-A' in Barcoding data/csv files/plotting order

####Generating the null distribution

set.seed(1234)

##Floral unit data
null.net_early_pooled_FU<-generate_null_net(pollinators_early_pooled[, 3:20], plants_early_FU[, 2:18], sims=10000,
                                         c.samples=pollinators_early_pooled[, 2],
                                         r.samples=plants_early_FU[, 1])

##Nectar data
null.net_early_pooled_nectar<-generate_null_net(pollinators_early_pooled[, 3:20], plants_early_nectar[, 2:18], sims=10000,
                                            c.samples=pollinators_early_pooled[, 2],
                                            r.samples=plants_early_nectar[, 1])
##Pollen data
null.net_early_pooled_pollen<-generate_null_net(pollinators_early_pooled[, 3:20], plants_early_pollen[, 2:18], sims=10000,
                                                c.samples=pollinators_early_pooled[, 2],
                                                r.samples=plants_early_pollen[, 1])

####Testing significance of interactions

#FU data
Early_interactions_test_pooled_FU <- test_interactions(null.net_early_pooled_FU, 0.98)
write.csv(Early_interactions_test_pooled_FU, "Interaction test_Early_P-A_Pooled species_FU.csv")

#Nectar data
Early_interactions_test_pooled_nectar <- test_interactions(null.net_early_pooled_nectar, 0.98)
write.csv(Early_interactions_test_pooled_nectar, "Interaction test_Early_P-A_Pooled species_nectar.csv")

#Pollen data
Early_interactions_test_pooled_pollen <- test_interactions(null.net_early_pooled_pollen, 0.98)
write.csv(Early_interactions_test_pooled_pollen, "Interaction test_Early_P-A_Pooled species_pollen.csv")


#####Analyzing fit of each null model to the data

##Install packages
library(nlme)
library(MuMIn)

model_deviation <- read.csv(file.choose()) #dataset showing the deviation between the interactions predicted by each null model and the observed interaction values for each individual plant species

#Running the GLMM
GLMM_deviation <- lme(absolute_deviation ~ null_model, random=~1|plant_species, method="REML", data = model_deviation)
summary(GLMM_deviation)

#Testing the fit of each null model
r.squaredGLMM(GLMM_deviation)
res.aov <- aov(absolute_deviation ~ null_model, data = proportional_deviation)
summary(res.aov)




#################################################################################################################
###### Obj.2: Measure diet breadth, interaction diversity & floral preferences and compare to null models #######
#################################################################################################################


###### Comparing diet breadth between sampling period, species and sites 

diet_div_data <- read.csv(file.choose()) 

#ANOVA for comparing diet breadth between sampling periods
diet_breadth_ANOVA1 <- aov(diet_diversity_observed ~ season, data = diet_div_data)
summary(diet_breadth_ANOVA1) # Summary of the analysis
TukeyHSD(diet_breadth_ANOVA1) ##Performs a parwise comparison on each independent variable category with adjustments for multiple testing 

#ANOVA for comparing diet breadth between species
diet_breadth_ANOVA2 <- aov(diet_diversity_observed ~ species, data = diet_div_data)
summary(diet_breadth_ANOVA2) # Summary of the analysis
TukeyHSD(diet_breadth_ANOVA2) ##Performs a parwise comparison on each independent variable category with adjustments for multiple testing 

#ANOVA for comparing diet breadth between sites
diet_breadth_ANOVA3 <- aov(diet_diversity_observed ~ farm, data = diet_div_data)
summary(diet_breadth_ANOVA3) # Summary of the analysis
TukeyHSD(diet_breadth_ANOVA3) ##Performs a parwise comparison on each independent variable category with adjustments for multiple testing 



#### Calculating Shannon diversity for the observed network as well as 95% CIs for the pollen and nectar-based null networks

net.stats_early_nectar <- bipartite_stats(null.net_early_pooled_nectar, index.type = "networklevel",
                                          indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")
write.csv(net.stats_early_nectar, "Network-metrics_early-season_P-A_nectar.csv")

net.stats_early_pollen <- bipartite_stats(null.net_early_pooled_pollen, index.type = "networklevel",
                                          indices=c("interaction evenness", "Shannon diversity"), intereven = "sum")
write.csv(net.stats_early_pollen, "Network-metrics_early-season_P-A_pollen.csv")



##### Calculating and plotting pollinator foraging preferences

svg("Preference_plot_Early-season_P-A_Nectar&pollen_0.98.svg", bg="transparent", width=6,height=6)

##Nectar data
plot_preferences(null.net_early_pooled_nectar, "Bombus.sp.", signif.level = 0.98, style = "dots",
                 type = "counts", res.order=plotting_order_early,  res.col = c("#67A9CF", "#F7F7F7", "#EF8A62"), 
                 l.cex = 1, p.cex = 1.5, font = 3,
                 xlab = "Number of visiting individuals")

par(new=TRUE)

##Pollen data
plot_preferences(null.net_early_pooled_pollen, "Bombus.sp.", signif.level = 0.98, style = "dots",
                 type = "counts", res.order=plotting_order_early,  res.col = c("#67A9CF", "#F7F7F7", "#EF8A62"), 
                 l.cex = 1, p.cex = 1.5, font = 3,
                 xlab = "Number of visiting individuals")
dev.off()



##########################################################################################################################
#######  Obj.3: Document phenological changes in bumblebee visitation patterns through the year ###########################
##########################################################################################################################


#### Running and Plotting NMDS on Barcoding Data 

##Load required packages
install.packages("vegan")
library(vegan)

##Set working directory
setwd("C:/Users/tt15117/OneDrive - University of Bristol/AA_Chapter 2_Pollen metabarcoding/EcoNullNetR results/Temporary results/NMDS")

##Load in required data
barcoding_data <- read.csv(file.choose()) ##Choose 'All seasons_Full barcoding results_P-A' from Documents\Pollen Barcoding\CSV files for R

##Setting data out as a matrix
barcoding_matrix <- barcoding_data[,5:178]
rownames(barcoding_matrix) <- barcoding_data[,1]
season <- barcoding_data[,2]
farm <- barcoding_data[,3]
insect <- barcoding_data[,4]

####Performing the NMDS 
Barcode_NMDS=metaMDS(barcoding_matrix,k=2,trymax=100) # K=The number of reduced dimensions, try max = number of iterations

Barcode_NMDS # Shows results of NMDS. 

stressplot(Barcode_NMDS) # Shepard plot, which shows scatter around the regression between the interpoint distances in the final configuration (distances between each pair of communities) against their original dissimilarities
# Large scatter around the line suggests that original dissimilarities are not well preserved in the reduced number of dimensions


###Plotting the NMDS 

par(mfrow=c(1,3))
svg("NMDS plot_All barcoding data_all 3 plots.svg", bg="transparent", width = 16, height=6)

##Plotting by season
ordiplot(Barcode_NMDS, type="n", xlim=c(-2.7,1.5), ylim=c(-1.2,1.5), bty="L", las = 1)
ordihull(Barcode_NMDS,groups=season,draw="polygon",col=c("darkorange", "deepskyblue4", "firebrick" ), alpha = 75,
         label=FALSE)
points(Barcode_NMDS, display = "sites", pch=c(15,16,17)[season], col=c("darkorange", "deepskyblue4", "firebrick")[season])

legend(0,2.35,
       legend = c("Early Season", "Mid Season", "Late Season"), #sets the legend labels to show
       title= "", title.adj=0.5,   #sets the legend title and shifts its position slightly
       pch=c(15,16,17), #sets the points style
       col= c("darkorange", "deepskyblue4", "firebrick" ), #sets the colours in the order they are listed above
       cex=0.8,
       bty = "n") #removes the box around the legend


##Plotting by farm
ordiplot(Barcode_NMDS, type="n", xlim=c(-2.7,1.5), ylim=c(-1.2,1.5), bty="L", las = 1)
ordihull(Barcode_NMDS,groups=farm,draw="polygon",col=c("orangered1", "turquoise4",  "red4"), alpha = 75,
         label=FALSE)
points(Barcode_NMDS, display = "sites", pch=c(15,16,17)[farm], col=c("orangered1", "turquoise4",  "red4")[farm])

legend(0,2.35,
       legend = c("Birches Farm", "Eastwood Farm", "Elmtree Farm"), #sets the legend labels to show
       title= "", title.adj=0.5,   #sets the legend title and shifts its position slightly
       pch=c(15,16,17), #sets the points style
       col= c("orangered1", "turquoise4",  "red4"), #sets the colours in the order they are listed above
       cex=0.8,
       bty = "n") #removes the box around the legend


### Plotting by pollinator species
ordiplot(Barcode_NMDS, type="n", xlim=c(-2.7,1.5), ylim=c(-1.2,1.5), bty="L", las = 1)
ordihull(Barcode_NMDS,groups=insect,draw="polygon",col=c("darkgreen", "darkorange3",  "red3", "deeppink3", "deepskyblue4"), alpha = 75,
         label=FALSE)
points(Barcode_NMDS, display = "sites", pch=c(15,16,17,18,19)[insect], col=c("darkgreen", "darkorange3",  "red3", "deeppink3", "deepskyblue4")[insect], cex=0.8)

legend(0,2.35,
       legend = c("B.hortorum", "B.hypnorum", "B.lapidarius", "B.pascuorum", "B.terrestris"), #sets the legend labels to show
       title= "", title.adj=0.5,   #sets the legend title and shifts its position slightly
       pch=c(15,16,17,18,19), #sets the points style
       col= c("darkgreen", "darkorange3",  "red3", "deeppink3", "deepskyblue4"), #sets the colours in the order they are listed above
       cex=0.8,
       
       bty = "n") #removes the box around the legend

dev.off()



##### Pollen-nectar phenology analysis

##Install packages
library(ggplot2)
library(nlme)
library(MuMIn)
library(emmeans)
library(multcomp)
library(contrast)
library(dplyr)
library(plotrix)

####Loading in data frames
pollen_nectar_phenology <- read.csv(file.choose()) #select 'Nectar_pollen rank_by season' 
pollen_nectar_ratios <- read.csv(file.choose()) #select 'Individual_bees_pollen-nectar ratios' 

pollen_nectar_phenology$season <- factor(pollen_nectar_phenology$season, levels = c("early","mid","late"))
pollen_nectar_ratios$season <- factor(pollen_nectar_ratios$season, levels = c("early","mid","late"))

###### ###### Summarising pollen and nectar real values for individual bees  ###### ######
pollen_nectar_values_grouped <- group_by(pollen_nectar_ratios, season) # group data by treatment 
pollen_nectar_values_summary <- summarize(pollen_nectar_values_grouped,
                                          mean_nectar = mean(mean_nectar_value_mg), # calculate the mean for each group
                                          sd_nectar = sd(mean_nectar_value_mg),
                                          se_nectar = std.error(mean_nectar_value_mg),
                                          mean_pollen = mean(mean_pollen_value_mm3), # calculate the mean for each group
                                          sd_pollen = sd(mean_pollen_value_mm3),
                                          se_pollen = std.error(mean_pollen_value_mm3)) # calculate the standard deviation for each group


##GLMM of pollen actual values
GLMM_pollen_value <- lme(mean_pollen_value_mm3 ~ season, random = list( ~1|farm,  ~1|insect), method="REML", data = pollen_nectar_ratios)
summary(GLMM_pollen_value)
r.squaredGLMM(GLMM_pollen_value)
#Post-hoc test
summary(glht(GLMM_pollen_value, emm(pairwise ~ season)), test=adjusted(type="bonferroni"))

##GLMM of nectar actual values
GLMM_nectar_value <- lme(mean_nectar_value_mg ~ season, random = list( ~1|farm,  ~1|insect), method="REML", data = pollen_nectar_ratios)
summary(GLMM_nectar_value)
r.squaredGLMM(GLMM_nectar_value)
#Post-hoc test
summary(glht(GLMM_nectar_value, emm(pairwise ~ season)), test=adjusted(type="bonferroni"))


###### ###### Plotting pollen vs nectar actual values for individual bees  ###### ######
pollen-nectar_plot <- ggplot(pollen_nectar_ratios, aes(x=mean_pollen_value_mm3, y=mean_nectar_value_mg, color=season, shape=season)) + 
  geom_point(size = 1  ) +  #Adds datapoints to plot
  geom_point(x=1.662, y=1.646, color="darkorange", shape=0, size=5)+ #mean FU values - taken from summary calculations above
  geom_point(x=2.159, y=0.88, color="deepskyblue4", shape=1,size=5)+ #mean FU values - taken from summary calculations above
  geom_point(x=0.603, y=0.374, color="firebrick", shape=2, size=5)+ #mean FU values - taken from summary calculations above
  geom_point(x=1.81, y=2.25, color="darkorange", shape=15, size=5)+ #mean bee values - taken from summary calculations above
  geom_point(x=1.81, y=2.06, color="deepskyblue4", shape=16,size=5)+ #mean bee values - taken from summary calculations above
  geom_point(x=1.92, y=2.92, color="firebrick", shape=17, size=5)+ #mean bee values - taken from summary calculations above
  theme_classic(base_size = 22, ) +  #sets plot theme 
  scale_colour_manual(values = c("darkorange", "deepskyblue4", "firebrick" )) +
  scale_shape_manual(values=c(15, 16, 17))+
  labs(x=bquote("Mean pollen value (mm" ^3* "/FU)"),
       y=bquote("Mean nectar value (mg/FU/day)"), size=0.5) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))
#ylim(18,50)+
#xlim(18,50)

svg("Polen_nectar_values_by season_individual bees.svg", bg="transparent", width=5.8,height=4)
pollen-nectar_plot
dev.off()


########################################################################################################################################
#######################################################   END OF SCRIPT  ###############################################################
########################################################################################################################################