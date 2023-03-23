library(minpack.lm)
library(ggplot2)
library(ggsignif)
library(gridExtra)
library(cowplot)
setwd("C:/Users/belov/Desktop/Ecophys PhD/Whole Plant Conductance (Ch.1)/Retention curve")
dat <- read.csv('VanGenuchten.csv')
dat$NegMPa <- -dat$MPa
dat$h <- - 10197.44 * dat$MPa
mod.fit <- nlsLM(GWC ~ theta_r + (theta_s - theta_r) / (1 + (alpha * h) ^ n) ^ (1 - 1/n),
                 data = dat,
                 start = list(theta_s = 0.28, theta_r = 0.005, alpha = 0.05, n = 2.8),
                 control = nls.lm.control(maxiter = 200))
dat$GWC.pred <- predict(mod.fit)
plot(dat$GWC, log10(dat$h))
coef <- as.matrix(coef(mod.fit), nrow = 2, ncol = 2)
theta_s <- coef[1]
theta_r <- coef[2]
alpha <- coef[3]
n <- coef[4]
m <- 1 - (1/n)

setwd("C:/Users/belov/Desktop/Ecophys PhD/Whole Plant Conductance 2019")
trans <- read.csv('Compiled_hourly_transpiration_WP_weather_data.csv')
trans.day <- subset(trans, Hour > 7 & Hour < 17)

figS4 <- ggplot(data = dat, aes(x = GWC, y = MPa)) +
  geom_point() +
  geom_line(data = trans.day, aes(x = gwc, y = WP)) +
  labs(x = 'GWC(%)', 
       y =  expression(paste(-Psi[s],' (MPa)'))) +
  ylim(-8, 0) +
  theme_bw() +
  annotate("text", x = 0.15, y = -2.1, label = "y = 0.00462 + (0.298 - 0.00462)/(1 + (0.00033 * h)^ 2.98)^(0.664)", 
           parse = FALSE, size = 4) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
        axis.title.x = element_text(vjust = 0, size = rel(1)),
        axis.text = element_text(size = rel(1.1)),
        axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.title = element_text(size = rel(0.8)),
        legend.text = element_text(size = rel(0.8)))

setwd("C:/Users/belov/Desktop/Ecophys PhD/Whole Plant Conductance (Ch.1)/Retention curve")
plotpath <- file.path("C:/Users/belov/Desktop/Ecophys PhD/Whole Plant Conductance (Ch.1)/Retention curve")
pdf(file = plotpath)
figS4
dev.off()