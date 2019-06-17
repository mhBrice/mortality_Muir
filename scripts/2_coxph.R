### PACKAGES ####
library(gtools)
library(dplyr)

library(survival)
library(survminer)
library(sjPlot)
library(rms)

library(RColorBrewer)
library(colourlovers)
library(scales)

library(graphicsutils)

### FUNCTIONS 

source('functions/risk_fun.R')
source('functions/mozaicplot.R')

### DATA ####

surv_sp <- readRDS("data/surv_sp.RDS")
data_muir2 <- readRDS("data/data_muir_clean.RDS")

species <- names(surv_sp)
spnames <- c('Acer saccharum', 'Fagus grandifolia', 'Ostrya virginiana', 'Tilia americana', 'Tsuga canadensis')


### FIGURE 1. MOSAIC PLOT ####

spnames2 <- c('Acer\nsaccharum', 'Fagus\ngrandifolia', 'Ostrya\nvirginiana', 'Tilia\namericana', 'Tsuga\ncanadensis')

col_surv <- rev(swatch(clpalette('4598774'))[[1]])[-1]

col_surv <- alpha(col_surv,.6)

dat <- subset(data_muir, Year_death >= 1998)
tab <- table(droplevels(dat$Species), droplevels(dat$Ice_storm), as.factor(dat$State))
dimnames(tab) <- list(c(spnames2), 
                      c("Ice storm 1", "Ice storm 2", "Ice storm 3", "Ice storm 4"),
                      c("Alive", "Dead"))

pdf("results/fig1_mosaic.pdf",  
    width = 6.8, height = 4.6)
# quartz(width = 6.8, height = 4.6)
par(mar = c(.5,4,2.5,3.8))
mozaicplot(tab = tab, col = col_surv)

dev.off()


# Percentage by damage class
# tr <- addmargins(table(droplevels(dat$Species), droplevels(dat$Ice_storm)))
# (tr[,4])/tr[,5]
# (tr[,3]+tr[,4])/tr[,5]


### FIGURE A1. MAP OF INDIVIDUAL TREES BY DAMAGE ####

xy_muir <- st_as_sf(dat, coords = c("X", "Y"), crs = 2145)
zone_muir <- st_read(dsn = 'data/zone_etude.shp')


title_ice <- c("Ice storm 1", "Ice storm 2", "Ice storm 3", "Ice storm 4")

pdf("results/figA1_map.pdf", 
    width = 6.2, height = 5)
# quartz(width = 6.2, height = 5)
par(mfrow=c(2,2), mar = c(.2,.2,.2,.2))
for(i in 1:4){
  plot(st_geometry(zone_muir), border = "grey", lwd = 2, xpd = NA)
  tmp <- subset(xy_muir, Ice_storm == i)
  plot(st_geometry(tmp), add=T, col = col_surv[i], pch = 19, xpd = NA)
  mtext(title_ice[i], 3, cex = .9, line = -1.5, adj = .4)
}
dev.off()


### FIGURE 2. MORTALITY RATES & SURIVAL CURVES ####


# Mean survival

col_surv <- swatch(clpalette('694737'))[[1]]

y_step <- c("1987-1996", "1996-1998","1998-1999","1999-2000", 
            "2000-2004", "2004-2011", "2011-2018")
fit <- surv_fit(Surv(Year_death, State) ~ Species,
                data = data_muir)

pdf("results/fig2_mortality.pdf",
    height = 3.5, width = 7)
# quartz(height = 3.5,width = 7)
par(mfrow =c(1,2), mar=c(4.5,2.7,1,.7))

# MORTALITY RATES
plot0(xlim=c(1995,2018), ylim=c(0,5))
axis(2,las=1, cex.axis = .7)
axis(1, at=cut_y, labels=F)
text(cut_y, rep(-1.2,7), y_step, srt=90, xpd=NA, cex=.7)
mtext("Mortality rate (%)", 2, cex = .9, line = 1.7)
mtext("Time interval (years)", 1, cex = .9, line = 3.5)
box2(1:2, lwd=1.5)
for(sp in species) {
  
  surv_tmp <- subset(SURV_sp, Species==sp)
  
  surv_tab  <- table(surv_tmp$time0, surv_tmp$State)
  surv_tab[,1] <- surv_tab[,1] + surv_tab[,2]
  
  lines((1-((surv_tab[,1]-surv_tab[,2])/surv_tab[,1])^(1/diff(c(1987,cut_y))))*100~cut_y, 
        type = "b",
        col = col_surv[which(species==sp)], lty=which(species==sp), 
        lwd = 1.5, pch = 19, cex = .7)
  
}
abline(v=1998, col = alpha("grey35",.8))
legend("topright", legend=spnames[c(3,2,1,4,5)], lty=c(3,2,1,4,5), lwd = 1.5, 
       col = col_surv[c(3,2,1,4,5)], cex = .7, bty="n", xpd = NA, inset = c(0,-.05))
mtext(letters[1], 3, adj = -.17, line = 0)

# SURVIVAL CURVES
plot(fit, xlim=c(1990,2020), ylim=c(.55,1), col = col_surv, 
     lty=1:5, lwd = 1.5, las = 1, axes = F)
box2(1:2, lwd = 1.5)
abline(v=1998, col = alpha("grey35",.8))
axis(1, cex.axis = .7)
axis(2, las = 1, cex.axis = .7)
mtext("Survival probability", 2, cex = .9, line = 2)
mtext("Time (years)", 1, cex = .9, line = 2)
mtext(letters[2], 3, adj = -.2)

dev.off()


### Get mortality rate by year ####

cut_y <- c(1996,1998,1999,2000, 2004, 2011, 2018) 

SURV_sp <- survSplit(data = data_muir, 
                     cut = cut_y, 
                     end = "Year_death", start = 'time0', zero = 1987,
                     event = "State")

# All species together
surv_tab  <- table(SURV_sp$time0, SURV_sp$State)
surv_tab[,1] <- surv_tab[,1] + surv_tab[,2]
(1-((surv_tab[,1] - surv_tab[,2])/surv_tab[,1])^(1/diff(c(1987,cut_y))))*100

# By species
x <- list()
for(sp in species) {
  surv_tmp <- subset(SURV_sp, Species==sp)
  surv_tab  <- table(surv_tmp$time0, surv_tmp$State)
  (surv_tab[,1] <- surv_tab[,1] + surv_tab[,2])
  x[[sp]] <- (1-((surv_tab[,1] - surv_tab[,2])/surv_tab[,1])^(1/diff(c(1987,cut_y))))*100
}
x

### Get overall mortality rates ####

SURV_sp <- survSplit(data = data_muir, 
                     cut = 2018, 
                     end = "Year_death", start = 'time0', zero = 1987,
                     event = "State")

# All species together
surv_tab  <- table(SURV_sp$time0, SURV_sp$State)
surv_tab[,1] <- surv_tab[,1] + surv_tab[,2]
(1-((surv_tab[,1] - surv_tab[,2])/surv_tab[,1])^(1/31))*100


# By species
x <- list()
for(sp in species) {
  surv_tmp <- subset(SURV_sp, Species==sp)
  surv_tab  <- table(surv_tmp$time0, surv_tmp$State)
  (surv_tab[,1] <- surv_tab[,1] + surv_tab[,2])
  x[[sp]] <- (1-((surv_tab[,1] - surv_tab[,2])/surv_tab[,1])^(1/31))*100
}
x


############################################
### SURVIVAL COX MODEL ####
############################################


### SCALE VARIABLES + CREATE SPECIES DF ####

var2scale <- c("Initial_DBH", "Hegyi_Acer", "Hegyi_Fagus", "Hegyi_Ostrya",  
               "Hegyi_Tilia",  "Hegyi_Tsuga",  "Tree_density", "Hegyi")

cut_y <- c(1998, 1999, 2000, 2004, 2011, 2018, 2020) 

for(sp in species) {
  surv_tmp <- surv_sp2[[sp]]
  surv_tmp[,var2scale] <- scale(surv_tmp[,var2scale])
  surv_tmp <- subset(surv_tmp, Year_death >= 1998 & !is.na(Ice_storm))
  surv_tmp$Ice_storm <- droplevels(surv_tmp$Ice_storm)
  # pour coxph
  surv_tmp$surv_object <- Surv(time = surv_tmp$Year_death, event = surv_tmp$State)
  
  assign(sp, surv_tmp)
}


### VARIABLES ####

var2mod <- c('Ice_storm',
             "Initial_DBH", "Tree_density",
             "Hegyi_Acer", "Hegyi_Fagus", "Hegyi_Ostrya",
             "Hegyi_Tilia",  "Hegyi_Tsuga")

var2fgr <- c('Ice_storm',
             "Initial_DBH", "Tree_density",
             "Hegyi_Acer", "Hegyi_Fagus", "Hegyi_Ostrya" ,
             "Hegyi_Tilia",  "Hegyi_Tsuga",
             "BBD")

var2plot <- c('Ice_storm2', 'Ice_storm3', 'Ice_storm4',
              "Initial_DBH", "Tree_density", 
              "Hegyi_Acer", "Hegyi_Fagus", "Hegyi_Ostrya" ,
              "Hegyi_Tilia",  "Hegyi_Tsuga",
              "BBD")

### FORMULAS COX-PH #### 

form_cox <- as.formula(paste0("surv_object ~ ", 
                              paste(var2mod, collapse = "+"))) 

form_cox2 <- as.formula(paste0("surv_object ~ ", 
                              paste(var2fgr, collapse = "+"))) 

# Test interaction
# form_cox <- surv_object ~ Verglas*Initial_DBH + Verglas*Tree_density + Hegyi_Acer +
#   Hegyi_Fagus + Hegyi_Ostrya + Hegyi_Tilia + Hegyi_Tsuga
# 
# form_cox2 <- surv_object ~ Verglas*Initial_DBH + Verglas*Tree_density + Hegyi_Acer +
#   Hegyi_Fagus + Hegyi_Ostrya + Hegyi_Tilia + Hegyi_Tsuga +
#   Verglas*BBD

# form_cox2 <-surv_object ~ Verglas + Initial_DBH + Tree_density + Hegyi + BBD

### RUN MODELS ####

(ph_asa = coxph(form_cox, data = ASA))
cox.zph(ph_asa)
(ph_fgr = coxph(form_cox2, data = FGR))
cox.zph(ph_fgr)
(ph_ovi = coxph(form_cox, data = OVI))
cox.zph(ph_ovi)
(ph_tam = coxph(form_cox, data = TAM))
(cox.zph(ph_tam))
(ph_tca = coxph(form_cox, data = TCA))
cox.zph(ph_tca)

# Check that Variance Inflation Factor (VIF) are below 10
vif(ph_asa)
vif(ph_fgr)
vif(ph_ovi)
vif(ph_tam)
vif(ph_tca)

# Get risk ratio + CI + pvalue stars for Cox PH

risk_ph_asa <- risk_ph(ph_asa, var2plot)
risk_ph_fgr <- risk_ph(ph_fgr, var2plot)
risk_ph_ovi <- risk_ph(ph_ovi, var2plot)
risk_ph_tam <- risk_ph(ph_tam, var2plot)
risk_ph_tca <- risk_ph(ph_tca, var2plot)

### TABLE A1-A5. MODEL SUMMARY ####

tab_model(ph_asa, title = "",  
          dv.labels = "Survival of Acer saccharum", 
          string.ci = "CI (95%)", string.p = "P-value", 
          file = "results/model_asa.html")

tab_model(ph_fgr, title = "",  
          dv.labels = "Survival of Fagus grandifolia", 
          string.ci = "CI (95%)", string.p = "P-value", 
          file = "results/model_fgr.html")

tab_model(ph_ovi, title = "",  
          dv.labels = "Survival of Ostrya virginiana", 
          string.ci = "CI (95%)", string.p = "P-value", 
          file = "results/model_ovi.html")

tab_model(ph_tam, title = "",  
          dv.labels = "Survival of Tilia americana", 
          string.ci = "CI (95%)", string.p = "P-value", 
          file = "results/model_tam.html")

tab_model(ph_tca, title = "",  
          dv.labels = "Survival of Tsuga canadensis", 
          string.ci = "CI (95%)", string.p = "P-value", 
          file = "results/model_tca.html")


### FIGURE 3. PLOT RISK RATIO ####

# Graphical param

labels.sig <- c('Ice storm 2','Ice storm 3', 'Ice storm 4',
                'Initial DBH', "Tree density", 
                'Hegyi Acer', 'Hegyi Fagus', 'Hegyi Ostrya', 
                'Hegyi Tilia', 'Hegyi Tsuga', 
                'BBD')


# Layout
m <- matrix(c(1,2,3,4,5), 5, 1, byrow = F)

# Plot

pdf("results/fig3_cox_risk_ratio.pdf", 
    width = 4.5, height = 7)
# quartz(width = 4.5, height = 7)
layout(m, heights = c(1,1,1,1,1))

par(mar=c(.4,3.5,.4,0.5), oma = c(5.1,0,0,0))

plot_risk(mod = ph_asa, var = var2plot, sp = "Acer saccharum", ylab = F)
plot_risk(mod = ph_fgr, var = var2plot, sp = "Fagus grandifolia", ylab = F)
plot_risk(mod = ph_ovi, var = var2plot, sp = "Ostrya virginiana")
plot_risk(mod = ph_tam, var = var2plot, sp = "Tilia americana", ylab = F)
plot_risk(mod = ph_tca, var = var2plot, lab = labels.sig, sp = "Tsuga canadensis", ylab = F)

dev.off()

### FIGURE 4. SURVIVAL CURVES FROM FITTED COX MODELS ####

mymodel <- list(ph_asa, ph_fgr, ph_ovi, ph_tam, ph_tca)

# Graphical parameters
col_surv <- rev(swatch(clpalette('4598774'))[[1]])[-1]

m <- matrix(c(1:10), 2, 5, byrow = F)

# Plot
pdf("results/fig4_survfit_col.pdf",
    height = 3.5, width = 7.48)
# quartz(height = 3.5, width = 7.48)
layout(m, widths = c(1,1,1,1,.5))
par(mar = c(1.2,1.2,1,.5), oma = c(2.1,2.1,1,0))

for(sp in 1:4) {
  dat <- dat_unsc[[sp]]
  
  mod <- coxph(form_cox, data = dat)
  
  # Verglas
  newdata <- data.frame(expand.grid(Ice_storm = levels(dat$Ice_storm),
                                    Initial_DBH = mean(dat$Initial_DBH),
                                    Tree_density = mean(dat$Tree_density),
                                    Hegyi_Acer = mean(dat$Hegyi_Acer),
                                    Hegyi_Ostrya = mean(dat$Hegyi_Ostrya),
                                    Hegyi_Fagus = mean(dat$Hegyi_Fagus),
                                    Hegyi_Tilia = mean(dat$Hegyi_Tilia),
                                    Hegyi_Tsuga = mean(dat$Hegyi_Tsuga),
                                    BBD = 0))
  
  plot(survfit(mod, newdata = newdata), 
       xlim = c(1998, 2018.2), ylim = c(0, 1), 
       col = col_surv, lty = 4:1,lwd = 1.3, 
       axes = F, xpd = NA)
  box2(1:2, lwd = 1.5)
  axis(1, labels = F)
  axis(2, labels = ifelse(sp==1, T, F), las = 1, xpd = NA, cex.axis = .9)
  mtext(spnames[sp], 3, line = .5, cex = .8, font = 3)
  
  
  # DBH
  
  mindbh <- quantile(dat$Initial_DBH, .1)
  maxdbh <- quantile(dat$Initial_DBH, .9)
  meddbh <- quantile(dat$Initial_DBH, .5)
  newdata <- data.frame(expand.grid(Ice_storm = "1",
                                    Initial_DBH = c(mindbh, meddbh, maxdbh),
                                    Tree_density = mean(dat$Tree_density),
                                    Hegyi_Acer = mean(dat$Hegyi_Acer),
                                    Hegyi_Ostrya = mean(dat$Hegyi_Ostrya),
                                    Hegyi_Fagus = mean(dat$Hegyi_Fagus),
                                    Hegyi_Tilia = mean(dat$Hegyi_Tilia),
                                    Hegyi_Tsuga = mean(dat$Hegyi_Tsuga),
                                    BBD = 0))
  
  plot(survfit(mod, newdata = newdata), 
       xlim = c(1998, 2018.2), ylim = c(0, 1), 
       col = col_surv, lty = 4:2, lwd = 1.5, 
       axes = F, xpd = NA)
  box2(1:2, lwd = 1.5)
  axis(1, cex.axis = .9)
  axis(2, labels = ifelse(sp==1, T, F), las = 1, xpd = NA, cex.axis = .9)
  
}
par(mar=c(0,0,0,0))
# Legend
plot0()
legend("center", legend = 1:4, col = col_surv, lty = 4:1, lwd = 1.5, bty = "n", 
       title = "Ice storm", xpd = NA)

plot0()
legend("center", legend = c("Small", "Medium", "Large"), 
       col = col_surv, lty = 4:2, lwd = 1.5, bty = "n", 
       title = "Initial DBH", xpd = NA)

mtext("Survival probability", 2, line = 1, outer = T, cex = .85)
mtext("Time (years)", 1, line = 1, outer = T, cex = .85)

dev.off()


#### FIGURE 5. BOXPLOT DAMAGE VS TREE SIZE #####

# Use unscale variable
for(sp in species) {
  surv_tmp <- surv_sp[[sp]]
  surv_tmp <- subset(surv_tmp, Year_death >= 1998 & !is.na(Ice_storm))
  surv_tmp$Ice_storm <- droplevels(surv_tmp$Ice_storm)
  surv_tmp$surv_object <- Surv(time = surv_tmp$Year_death, event = surv_tmp$State)
  
  assign(paste0(sp, "_unsc"), surv_tmp)
}

dat_unsc <- list(ASA_unsc, FGR_unsc, OVI_unsc, TAM_unsc, TCA_unsc)


# Layout
m <- matrix(c(1:5), 1, 5, byrow = F)


pdf("results/fig5_damageVSsize.pdf", 
    height = 2.3, width = 8.48)
# quartz(height = 2.3, width = 8.48)
layout(m, widths = c(1,1,1,1,1))
par(mar = c(1.2,1.2,1.2,.5), oma = c(2.1,2.1,1,0))
for(sp in 1:5) {
  dat <- dat_unsc[[sp]]
  
  aov_sp <- summary(aov(Initial_DBH ~ Ice_storm, data = dat))
  pval <- aov_sp[[1]][1,5]
  star <- stars.pval(pval)
  pval <- ifelse(pval < 0.001, "< 0.001", signif(pval, 1))
  
  boxplot2(Initial_DBH ~ Ice_storm, data = dat)
  box2(1:2, lwd = 1.5)
  axis(1, at = as.numeric(levels(dat[,"Ice_storm"])), cex.axis = .85)
  axis(2, las = 1, cex.axis = .85)
  mtext(spnames[sp], 3, line = 1, cex = .8, font = 3)
  mtext(paste0("p-value ", pval, star), 3, line = -.2, cex = .7)
}

mtext("Initial DBH (cm)", 2, line = 1, outer = T, cex = .85)
mtext("Levels of canopy loss from the ice storm", 1, line = 1, outer = T, cex = .85)

dev.off()


#### FIGURE A2. BOXPLOT DAMAGE VS TREE DENSITY #####


pdf("results/figA2_damageVSdensity.pdf", 
    height = 2.3, width = 8.48)
# quartz(height = 2.3, width = 8.48)
layout(m, widths = c(1,1,1,1,1))
par(mar = c(1.2,1.2,1.2,.5), oma = c(2.1,2.1,1,0))
for(sp in 1:5) {
  dat <- dat_unsc[[sp]]
  
  aov_sp <- summary(aov(Tree_density ~ Ice_storm, data = dat))
  pval <- aov_sp[[1]][1,5]
  star <- stars.pval(pval)
  pval <- ifelse(pval < 0.001, "< 0.001", signif(pval, 1))
  
  boxplot2(Tree_density ~ Ice_storm, data = dat)
  box2(1:2, lwd = 1.5)
  axis(1, at = as.numeric(levels(dat[,"Ice_storm"])), cex.axis = .85)
  axis(2, las = 1, cex.axis = .85)
  mtext(spnames[sp], 3, line = 1, cex = .8, font = 3)
  mtext(paste0("p-value ", pval, star), 3, line = -.2, cex = .7)
}

mtext("Tree density", 2, line = 1, outer = T, cex = .85)
mtext("Levels of canopy loss from the ice storm", 1, line = 1, outer = T, cex = .85)

dev.off()


