


model = lm(sum_pop_size ~ H + D + H:D,data=model_output)
require(car)
Anova(model,type='II')


interaction.plot(x.factor     = model_output$H,
                 trace.factor = model_output$D,
                 response     = model_output$sum_pop_size,
                 fun = mean)


# ANOVA assumptions
hist(model_output$sum_pop_size[model_output$H==0 & model_output$D==0.01])


model_output$H = as.factor(model_output$H)
model_output$D = as.factor(model_output$D)
mymodel = anova(lm(sum_pop_size ~ H*D ,data=model_output))

# Tests for homogeneity of variances. It shows that variances are different between groups, violating ANOVA assumption
leveneTest(lm(sum_pop_size ~ H*D,data=model_output))

bartlett.test(sum_pop_size ~ H + D,data=model_output)
#Could use repeated measures ANOVA to look a time series of population data - it is like sampling the same human over their life to see how their health changes


pairwise.t.test(model_output$sum_pop_size, model_output$H + model_output$D,p.adjust.method = 'bonferroni')

