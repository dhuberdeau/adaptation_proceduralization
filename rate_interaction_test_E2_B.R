library(lme4)
library(boot)
library(effects)
## load data, rotation B (08/25/2017)
E1_rate <- read.csv("~/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_proceduralization/rotationB_output_matrix", header = FALSE, na.strings = 'NaN')
E1_cat <- read.csv("~/OneDrive/Documents/JHU/BLAM_lab/Matlab/adaptation_proceduralization/rotationB_factor_matrix", header = FALSE, na.strings = 'NaN')

dat <- data.frame(cycle = E1_cat[,1],
                  type = E1_cat[,2], 
                  subject = E1_cat[,3],
                  Y1 = E1_rate[,1])

dat$cycle <- as.numeric(dat$cycle)
dat$type <- factor(dat$type)
dat$subject <- factor(dat$subject)
dat$Y1 <- as.numeric(dat$Y1)
lm1 <- lmer(Y1 ~ type + cycle + (1|subject), data = dat, REML = FALSE)
lm2 <- update(lm1, formula = . ~ . + type:cycle)
lm3 <- lmer(Y1 ~ type + (1|subject), data = dat, REML = FALSE)
lm4 <- lmer(Y1 ~ cycle + (1|subject), data = dat, REML = FALSE)

anova(lm1, lm2) # LRT
anova(lm1, lm3)
anova(lm1, lm4)
plot(allEffects(lm2))

