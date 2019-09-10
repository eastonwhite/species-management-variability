# Main manuscript figures for:
# Ecological management success and failure is highly stochastic: an experimental test
# Created by Easton R. White
# Last edited: 2-Sep-2019


###########################################
###########################################


# Figure 1: Main simulation results
model_output <- read.csv('../model_outputs/simulation_results2.csv',header=T)

# Build figure 1 with model_output
require(dplyr)
require(ggplot2)
require(ggpubr)
require(ggthemes)
require(viridis)


model_totals = model_output %>%
  group_by(gen,landscape,hRate,dTime) %>%
  summarize(total_pop = sum(count))

# Keeping pop less than 500 beetles by gen 6
model_totals_success = model_totals %>%
  filter(gen==6) %>%
  group_by(landscape,hRate,dTime) %>%
  summarize(success = sum(total_pop<300)) %>%
  group_by(hRate,dTime) %>%
  summarize(percent_success = sum(success)/max(landscape))

model_totals_success$abs_value = abs(model_totals_success$percent_success-0.5)
model_totals_success$abs_value = (model_totals_success$abs_value- min(model_totals_success$abs_value))/(max(model_totals_success$abs_value)-min(model_totals_success$abs_value))


fig1 <- ggplot(aes(x=hRate,y=percent_success,group=dTime,color=dTime),data=model_totals_success) + geom_line() + scale_color_viridis(begin = 0,end=0.9) + theme_tufte() +theme(axis.text = element_text(size=12),text=element_text(size=16)) +
  ylab(' ') +
  xlab(" ") +
  labs(color="Dispersal time (hours)") +
  theme_classic(base_size = 14, base_family = "") +
  ylim(0,1) +
  annotate("text",x=0.8,y=0.02,label='(a)',col='black')


# Plot 2

prob_prevent_spread = model_output %>%
  filter(gen==4,patch %in% c(4,5,6)) %>%
  group_by(landscape,hRate,dTime) %>%
  summarize(success = sum(count)) %>%
  group_by(hRate,dTime) %>%
  summarize(num_success = sum(success==0))
#BC_totals_failed$num_replicates = BC_totals_failed_tally$n
prob_prevent_spread$percent_success = prob_prevent_spread$num_success/max(prob_prevent_spread$num_success)


fig2 <- ggplot(aes(x=hRate,y=percent_success,group=dTime,color=dTime),data=prob_prevent_spread) +   geom_line() + scale_color_viridis(begin = 0,end=0.9) +
  ylab(' ') +
  xlab(" ") +
  labs(color="Dispersal time (hours)") +
  theme_classic(base_size = 14, base_family = "")+
  ylim(0,1) +
  annotate("text",x=0.8,y=0.02,label='(b)',col='black')

pdf(file='figures/stochastic_model_outputs.pdf',width=8,height=4)

# Combine plots
figure=ggarrange(fig1, fig2, ncol=2, nrow=1, common.legend = TRUE, legend="top")
annotate_figure(figure,
                bottom = text_grob("Harvest rate", color = "black",
                                   hjust = 0.9,vjust=-1, x = 0.6, size = 16),
                left = text_grob("Probability of success", color = "black", rot = 90,size=16,vjust=2,hjust=0.4)
)


dev.off()




###########################################
###########################################



# Build figure 2

pdf(file = 'figures/spread_rate.pdf',width = 5,height=5)

model_output <- read.csv('../model_outputs/simulation_results_spatial_spread2.csv',header=T)

model_output$trial <- ceiling(model_output$landscape/10)

model_output$dTime = as.factor(model_output$dTime)
levels(model_output$dTime)= c('Low','Med','High')

# Pull in simulation results and build plot

# Need mean and sd of each bunch
spread_speed = model_output %>%
  filter(gen==4) %>%
  group_by(hRate,dTime,landscape,trial) %>%
  summarise(max_patch = ifelse(length(which(count>0)),patch[tail(which(count>0),1)],0.001)/4 ) %>% # Artifically put in a value of 0.1 if entire landscape goes extinct...
  group_by(hRate,dTime,trial) %>% 
  summarise(avg_dist_spread = mean(max_patch),sd_dist_spread = sd(max_patch)) %>%
  group_by(hRate,dTime) %>%
  summarise(avg_dist_spread = mean(avg_dist_spread),avg_sd_dist_spread = mean(sd_dist_spread))

require(ggplot2)
require(viridis)
p <- ggplot(data=spread_speed) +
  geom_ribbon(aes(x=hRate,ymin=avg_dist_spread-2*avg_sd_dist_spread,ymax=avg_dist_spread+2*avg_sd_dist_spread,fill=as.factor(dTime))) + scale_fill_viridis(discrete=T,alpha=0.1,guide=F) +
  geom_line(aes(x=hRate,y=avg_dist_spread,color=as.factor(dTime)),size=1.5) + scale_color_viridis(discrete = T,alpha=0.5) +
  labs(color="Dispersal") +
  xlab("Harvest rate") +
  ylab("Spread rate") +
  theme_classic(base_size = 14, base_family = "")

# Add experimental data to the plot
BC <- read.csv('../clean_data/BeetleControl_CompleteData_clean.csv',stringsAsFactors = FALSE)

BC$dTime = as.factor(BC$dTime)
levels(BC$dTime)= c('Low','Med','High')

require(dplyr)

spread_speed_BC = BC %>%
  filter(gen==4) %>%
  group_by(hRate,dTime,landscape) %>%
  summarise(max_patch = ifelse(length(which(count>0)),patch[tail(which(count>0),1)],0.001)/4 ) 
# %>%
#   group_by(hRate,dTime) %>%
#   summarise(avg_dist_spread = mean(max_patch),sd_dist_spread = sd(max_patch))

#
p + geom_point(aes(x=jitter(hRate,0.5),y=jitter(max_patch,0.5),col=as.factor(dTime)),data=spread_speed_BC,shape=12)  +   theme(legend.position="top")

#scale_y_continuous(limits = c(-0.05,1.6)) + scale_x_continuous(limits = c(-0.05,0.79)) +

# +
#   geom_errorbar(aes(x=hRate,ymin= avg_dist_spread - sd_dist_spread,ymax= avg_dist_spread + sd_dist_spread),data=spread_speed_BC)

dev.off()






###########################################
###########################################

# Figure 3

# Plot of total population size over time for different treatments, also report percent of successful treatments...
pdf(file = 'figures/variability_in_controlling_pop_size.pdf',width = 7,height=5)

BC <- read.csv('../clean_data/BeetleControl_CompleteData_clean.csv',stringsAsFactors = F)

BC$dTime = as.factor(BC$dTime)
BC$hRate = as.factor(BC$hRate)
levels(BC$dTime)= c('low','med','high')
levels(BC$hRate)= c('none','low','med','high')

BC=na.omit(BC) # Why na.omit?
require(dplyr)
require(ggplot2)
require(viridis)

threshold = 300

BC_totals = BC %>%
  group_by(gen,landscape,hRate,dTime) %>%
  summarize(total_pop = sum(count))

BC_totals_gen6 = BC_totals %>%
  filter(gen==6) %>%
  group_by(hRate,dTime) %>%
  summarize(mean_total_pop = mean(total_pop))
# This result shows that on average a control strategy might meet the target, but there still might only be a 70% of the strategy actually meeting the target

BC_totals_failed = BC_totals %>%
  filter(gen==6) %>%
  group_by(landscape,hRate,dTime) %>%
  summarize(success = sum(total_pop<threshold)) %>%
  group_by(hRate,dTime) %>%
  summarize(percent_success = sum(success)/8)
BC_totals_failed$x=0
BC_totals_failed$y=1500
BC_totals_failed$landscape=1


p=ggplot(data = BC_totals,aes(x=gen,y=total_pop)) + geom_point(aes(color=as.factor(landscape),by=as.factor(landscape))) +geom_line(aes(color=as.factor(landscape),by=as.factor(landscape))) +facet_grid(hRate~dTime) +  theme_classic(base_size = 14, base_family = "") +
  scale_color_viridis(discrete = TRUE,begin = 0,end=0.7) +
  labs(color='Replicate') +
  xlab("Time") +
  ylab("Number of beetles") +
  annotate("rect", xmin = -0.5, xmax = 6.6, ymin = threshold, ymax = 2500,alpha = .1,col='lightgrey')+   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background = element_blank(),plot.margin = unit(c(1,3,1,1), "lines"))+
  coord_cartesian(ylim=c(0, 2000),xlim=c(0,6))

p+geom_text(mapping = aes(x,y,label=paste('Pr(Success)=',percent_success)),hjust=0,vjust=0,data=BC_totals_failed,size=2.5)  


dev.off()

###########################################
###########################################

# Statistical model to see how probability of success depends on dispersal and harvest rates
BC_totals$success = 0
BC_totals$success[BC_totals$total_pop<threshold]=1

ml=glm(formula = BC_totals$success~BC_totals$hRate+BC_totals$dTime +BC_totals$hRate*BC_totals$dTime + BC_totals$gen, family = binomial(link='logit'))
summary(ml)
exp(coef(ml))
# Use modern dive package and function for pretty table

explained_deviance = 100*(ml$null.deviance-ml$deviance)/ml$null.deviance

require(stargazer)
stargazer(ml, title="Regression Results", dep.var.labels=c("Probability Success"),covariate.labels=c("Harvest Rate","Dispersal Time","Generation","Harvest*Dispersal","Intercept"),omit.stat=c("LL","ser","f"), no.space=TRUE,ci=TRUE,model.numbers = FALSE,add.lines = list(c("Explained Deviance", paste(round(explained_deviance,2),'\\%',sep = ''))),notes = "I make this look good!")









###########################################
###########################################




# Figure 4


pdf(file = 'figures/variability_in_controlling_spread.pdf',width = 7,height=5)

BC <- read.csv('../clean_data/BeetleControl_CompleteData_clean.csv',stringsAsFactors = F)

BC$dTime = as.factor(BC$dTime)
BC$hRate = as.factor(BC$hRate)
levels(BC$dTime)= c('low','med','high')
levels(BC$hRate)= c('none','low','med','high')

require(ggplot2)
require(viridis)
BC_totals_failed_tally = BC %>%
  filter(gen==4,patch %in% c(4,5,6)) %>%
  group_by(landscape,hRate,dTime) %>%
  summarize(success = sum(count)) %>%
  group_by(hRate,dTime) %>%
  tally

BC_totals_failed = BC %>%
  filter(gen==4,patch %in% c(4,5,6)) %>%
  group_by(landscape,hRate,dTime) %>%
  summarize(success = sum(count)) %>%
  group_by(hRate,dTime) %>%
  summarize(num_success = sum(success>0))
BC_totals_failed$x=2.1
BC_totals_failed$y=400
BC_totals_failed$landscape=1
BC_totals_failed$num_replicates = BC_totals_failed_tally$n
BC_totals_failed$percent_success = BC_totals_failed$num_success/BC_totals_failed$num_replicates

BC_gen4 = subset(BC,BC$gen==4)
p=ggplot(data = BC_gen4,aes(x=patch,y=count)) + geom_point(aes(color=as.factor(landscape))) +geom_line(aes(color=as.factor(landscape))) +facet_grid(hRate~dTime)  +  theme_classic(base_size = 14, base_family = "") +
  scale_color_viridis(discrete = TRUE,begin = 0,end=0.7) +
  labs(color='Replicate') +
  xlab("Patch") +
  ylab("Number of beetles") +
  annotate("rect", xmin = 3.8, xmax = 6.6, ymin = -10, ymax = 600,alpha = .1,col='lightgrey')+   
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),strip.background = element_blank(),plot.margin = unit(c(1,3,1,1), "lines"))+
  coord_cartesian(ylim=c(0,500),xlim=c(1,6))

p+geom_text(mapping = aes(x,y,label=paste('Pr(Success)=',1-percent_success)),hjust=0,vjust=0,data=BC_totals_failed,size=2.5)

dev.off()




###############


# Statistical outputs for fig 4
BC_stopped_dispersal = BC %>%
  filter(gen==4,patch %in% c(4,5,6)) %>%
  group_by(landscape,hRate,dTime) %>%
  summarize(total_pop = sum(count))

BC_stopped_dispersal$success = 0
BC_stopped_dispersal$success[BC_stopped_dispersal$total_pop<1]=1

# Statistical model
ml=glm(formula = BC_stopped_dispersal$success~BC_stopped_dispersal$hRate*BC_stopped_dispersal$dTime, family = binomial(link='logit'))
#ml1=glm(formula = BC_stopped_dispersal$success~BC_stopped_dispersal$hRate*BC_stopped_dispersal$dTime, family = binomial(link='probit'))
summary(ml)
exp(coef(ml))


explained_deviance = 100*(ml$null.deviance-ml$deviance)/ml$null.deviance
