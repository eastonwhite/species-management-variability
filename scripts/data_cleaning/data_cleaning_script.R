# Created by Easton R. White (eastonrwhite@gmail.com)
# Last edited: 2-Sep-2019
# BeetleControl project


# Cleaning of BeetleControl experimental data

BC <- read.csv('../raw_data/BeetleControl_CompleteData.csv',stringsAsFactors = FALSE)

# Add a generation 0 with same initial starting conditions
BC_gen0 = subset(BC,BC$gen==1)
BC_gen0$gen = 0
BC_gen0[,c('total_count','total_w50','total_weight','harvest_count','harvest_weight')] =0
BC_gen0$total_count[BC_gen0$patch==1] = 20

BC = rbind(BC_gen0,BC)

BC = subset(BC,BC$landscape %in% 1:10)
BC$count = BC$total_count
BC$count[is.na(BC$count)] = 50*BC$total_weight[is.na(BC$count)]/BC$total_w50[is.na(BC$count)]

# Add harvest data for generation 6 (was not written down as there was no generation 7)
BC$harvest_count[BC$harvest_weight>0] = NA
BC$harvest_count[BC$gen==6 & BC$landscape > 6] = NA
BC$harvest_count[BC$gen==6 & BC$landscape > 6] = ceiling(BC$hRate[BC$gen==6 & BC$landscape > 6]*BC$total_count[BC$gen==6 & BC$landscape > 6])
BC$harvest_weight[BC$gen==6 & BC$landscape > 6] = NA
BC$harvest_weight[BC$gen==6 & BC$landscape > 6] = BC$hRate[BC$gen==6 & BC$landscape > 6]*BC$total_weight[BC$gen==6 & BC$landscape > 6]
BC$harvest_weight[BC$gen==6 & BC$landscape > 6 & is.na(BC$harvest_count)==F]=0

# Create harvest coluom
BC$harvest = BC$harvest_count
BC$harvest[is.na(BC$harvest)] = ceiling(BC$count[is.na(BC$harvest)]*BC$harvest_weight[is.na(BC$harvest)]/BC$total_weight[is.na(BC$harvest)])

#BC$harvest[is.na(BC$count)] = NA


#BC[BC$gen>1 & BC$hRate==0 & BC$dTime==6 & BC$landscape==1,c('total_count','total_w50','total_weight','harvest_count','harvest_weight','count','issue','harvest')] = NA
# Replace data point where "zero harvest" treatment was recorded as harvested. The next generation data support that nothing was harvesting. Therefore we just need to set the harvest to 0 for that treatment
BC[BC$gen>1 & BC$hRate==0 & BC$dTime==6 & BC$landscape==1,c('harvest_count','harvest_weight','harvest')] = 0

# Remove replicates 7 and 8 after generation 4 because of mistake with counting and harvesting
BC[BC$gen %in% c(5,6) & BC$landscape %in% c(7:8),c('total_count','total_w50','total_weight','harvest_count','harvest_weight','count','issue')] = NA


# Landscape 8 had a weird large value in patch 3, gen 4. I have left that value in for now as there was no explanation as to why.
# Look at this with BC[BC$gen %in% c(4) & BC$dTime %in% c(1) & BC$hRate==0.5 & BC$patch==3,]


# Optional code to remove NAs from the dataset
#BC = subset(BC,is.na(BC$issue)==F)


# Write clean data file

write.csv(BC[,c('hRate','dTime','landscape','patch','gen','count','harvest')],file = "../clean_data/BeetleControl_CompleteData_clean.csv",quote = F,row.names = F)







#############################################
#############################################
# Tests to make sure the data looks clean 

#plot(BC$harvest,BC$hRate*50*BC$total_weight/BC$total_w50) 
# Can see rogue point at BC[BC$gen==2 & BC$hRate==0 & BC$dTime==6 & BC$landscape==1 & BC$patch==1,]

#plot(BC$count,(1-BC$hRate)*50*BC$total_weight/BC$total_w50)

#BC[BC$hRate==0 & BC$harvest>0,]


