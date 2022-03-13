d <- read.csv("data/raw_data.csv", stringsAsFactors = FALSE)

norm_quantile_map <- function(v) {
  v_quant <- rank(v, na.last = "keep") / (length(v) + 1)
  quant_out <- qnorm(v_quant)
  quant_out
}

table(d$Sample)
d$sample_name <- d$Sample
d$sample_name <- gsub("New ", "", d$sample_name)
d$sample_name <- gsub("Old ", "", d$sample_name)
d$log_conc <- log2(d$Calc..Concentration)

d$Assay <- gsub("IL-1\\?", "IL-1B", d$Assay)
d$Assay <- gsub("TNF-\\?", "TNF-A", d$Assay)
assay_unique <- unique(d$Assay)
n_assay <- length(assay_unique)
d$geno <- sapply(strsplit(d$sample_name, split = " "), function(x) x[1])
d$bio <- sapply(strsplit(d$sample_name, split = " "), function(x) x[2])
d$id <- sapply(strsplit(d$sample_name, split = " "), function(x) x[3])
d$geno_bio <- paste0(d$geno, "_", d$bio)
geno_bio_unique <- c("WT_blood", "Dp5_blood", "Dp5_fluid")
d$geno_bio <- factor(d$geno_bio, levels = geno_bio_unique)
d$litter <- substr(d$id, 1, 4)
d$Plate <- C(as.factor(d$Plate), contr = contr.treatment)
# d$Plate <- as.factor(d$Plate)
# d <- d[!is.na(d$log_conc), ]
mean(is.na(d$log_conc))
assay_curr <- assay_unique[1]
anova_pval <- c()
t_out_names <- c("Value", "Std.Error", "DF", "t-value", "p-value")
t_out_names <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
t_list <- list()
for (t_out_name in t_out_names) {
  t_list[[t_out_name]] <- matrix(NA, n_assay, 6)
  dimnames(t_list[[t_out_name]]) <- list(assay_unique, 
                                         c("Dp5_blood_v_WT_blood", "Dp5_fluid_v_WT_blood", "Dp5_fluid_v_Dp5_blood", 
                                           "WT_blood", "Dp5_blood", "Dp5_fluid"))
}  


range_of_detection <- list('IL-17A' = c(0.30, 2150),
                           'IL-21' =  c(6.5, 40600),
                           'IL-22' = c(1.2, 1850),
                           'IL-6' = c(4.8, 16000),
                           'IL-10' = c(3.8, 22800), 
                           'IL-1B' = c(3.1, 13000),
                           'TNF-A' = c(1.3, 6200),
                           'VEGF-A' = c(0.77, 12100))


paste(paste0(names(range_of_detection), " = ", sapply(range_of_detection, function(x) x[1]), ""), collapse = ", ")

d$log_conc_lod <- d$log_conc
d$conc_lod <- d$Calc..Concentration
d$conc_log_fc <- NA
for (assay_curr in assay_unique) {
  lod_lower <- (range_of_detection[[assay_curr]][1])#min(d$log_conc[d$Assay == assay_curr], na.rm = TRUE)
  lod_upper <- (range_of_detection[[assay_curr]][2])
  which_replace_lower <- which(d$Assay == assay_curr & (is.na(d$conc_lod) | d$conc_lod < lod_lower))
  d$conc_lod[which_replace_lower] <- lod_lower * .5
  which_replace_upper <- which(d$Assay == assay_curr & d$conc_lod > lod_upper)
  which_this_assay <- which(d$Assay == assay_curr) 
  d$conc_log_fc[which_this_assay] <- log2(d$conc_lod[which_this_assay])# / lod_lower)
}

plate_unique <- unique(d$Plate)
plate_mns <- sapply(plate_unique, function(x) mean(d[d$Plate == x, "conc_log_fc"]))
d$conc_log_fc_plate_corr <- d$conc_log_fc - plate_mns[d$Plate] + mean(plate_mns)



d$sample_plate <- paste0(d$sample_name, "_", d$Plate)
sample_plate_unique <- unique(d$sample_plate)
n_sample_plate <- length(sample_plate_unique)
d$sample_plate_assay <- paste0(d$sample_name, "_", d$Plate, "_", d$Assay)
sample_plate_assay_unique <- unique(d$sample_plate_assay)
n_sample_plate_assay <- length(sample_plate_assay_unique)
sample_name_unique <- unique(d$sample_name)
n_sample <- length(sample_name_unique)
d$sample_assay <- paste0(d$sample_name, "_", d$Assay)
sample_assay_unique <- unique(d$sample_assay)
n_sample_assay <- length(sample_assay_unique)
d_sample <- d[match(sample_assay_unique, d$sample_assay), ]
d_sample$mean_log_conc <- sapply(d_sample$sample_assay, function(x) mean(d[d$sample_assay == x, "conc_log_fc_plate_corr"], na.rm = T))
d_sample_plate <- d[match(sample_plate_unique, d$sample_plate), ]
d_sample_plate_assay <- d[match(sample_plate_assay_unique, d$sample_plate_assay), ]
d_sample_plate$mean_log_conc <- sapply(d_sample_plate$sample_plate, function(x) mean(d[d$sample_plate == x, "conc_log_fc_plate_corr"], na.rm = T))
d_sample_plate_assay$mean_log_conc <- sapply(d_sample_plate_assay$sample_plate, function(x) mean(d[d$sample_plate == x, "conc_log_fc_plate_corr"], na.rm = T))

na_prop <- matrix(NA, n_sample_plate, n_assay, dimnames = list(sample_plate_unique, assay_unique))
for (sample_plate in sample_plate_unique) {
  for (assay in assay_unique) {
    na_prop[sample_plate, assay] <- sum(is.na(d[d$sample_plate == sample_plate & d$Assay == assay, "log_conc"]))
  }
} 

graphics.off()
dir.create("figures", showWarnings = FALSE)
pdf(file = "figures/missingness.pdf", 12, 6)
par(oma = c(10, 8, 2, 2))
colv <- grey(c(1, .5, 0))
image(x = 1:n_sample_plate, y = 1:n_assay, z = na_prop, xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = colv, bty = "n")
axis(side = 1, at = 1:n_sample_plate, labels = sample_plate_unique, las = 2, cex.axis = .7)
axis(side = 2, at = 1:n_assay, labels = assay_unique, las = 2, cex.axis = .7)
par(xpd = NA)
legend(x = -n_sample_plate / 6, y = n_assay, legend = c(0, 1, 2), col = colv, title = "# NaN", pch = 19, pt.cex = 1.3)
par(xpd = F)
dev.off()

assay_unique
paired_test_pval <- 
for (assay_curr in assay_unique) {
  for (geno_bio in geno_bio_unique) {
    
  }
}


table(table(d$sample_name) / 8)

test_df <- data.frame(test = c("Dp5_blood_vs_WT_blood", "Dp5_fluid_vs_WT_blood", "Dp5_fluid_vs_Dp5_blood"))
test_df$geno_bio_1 <- sapply(strsplit(test_df$test, split = "_vs_"), function(x) x[1])
test_df$geno_bio_2 <- sapply(strsplit(test_df$test, split = "_vs_"), function(x) x[2])
nonpar_pval <- matrix(NA, n_assay, 3, dimnames = list(assay_unique, test_df$test))
for (assay_curr in assay_unique) {
  d_sample_curr <- d_sample[d_sample$Assay == assay_curr, ]
  for (j in 1:nrow(test_df)) {
    y_curr_1 <- d_sample_curr[d_sample_curr$geno_bio == test_df$geno_bio_1[j], "mean_log_conc"]
    y_curr_2 <- d_sample_curr[d_sample_curr$geno_bio == test_df$geno_bio_2[j], "mean_log_conc"]
    nonpar_pval[assay_curr, test_df$test[j]] <- wilcox.test(y_curr_1, y_curr_2)$p.value
  }
}
nonpar_qval <- nonpar_pval
nonpar_qval[] <- p.adjust(c(nonpar_pval), method = "BH")

graphics.off()
RColorBrewer::display.brewer.all(n=NULL, type="all", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=TRUE)

conc_range <- range(d_sample$mean_log_conc)
ylim_log <- range(d_sample$mean_log_conc, na.rm = T) + c(-1, 3)
colv <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
names(colv) <- assay_unique
jit_eps1 <- .15
jit_eps <- .4
# pdf("figures/barplot.pdf", 12, 7)
output_plot <- TRUE
if (output_plot) {
  jpeg("figures/barplot.jpeg", 12, 7, units = "in",  res = 500)
}
par(mfrow = c(1, 8), mar = c(0, 0, 0, 0), oma = c(10, 6, 4, 4))
assay_curr <- assay_unique[2]
for (assay_curr in assay_unique) {
  d_sample_curr <- d_sample[d_sample$Assay == assay_curr, ]
  plot(0, ty = "n", xlim = c(.5, 3.5), ylim = ylim_log, xaxt = "n", yaxt = "n", yaxs = "i")
  mtext(side = 3, text = assay_curr, col = colv[assay_curr], line = 1)
  for (geno_bio in geno_bio_unique) {
    y_curr <- d_sample_curr[d_sample_curr$geno_bio == geno_bio, "mean_log_conc"]
    mn_curr <- mean(y_curr)
    se_curr <- sd(y_curr) / sqrt(length(y_curr))
    ind_curr <- match(geno_bio, geno_bio_unique)
    polygon(x = ind_curr + c(-1, 1, 1, -1, -1) * jit_eps, 
            y = c(ylim_log[1], ylim_log[1], mn_curr, mn_curr, ylim_log[1]), 
            col = colv[assay_curr])
    points(x = ind_curr + runif(n = length(y_curr), -jit_eps1, jit_eps1), 
           y = y_curr, bg = 1, col = 1, pch = 21)
    lines(x = rep(ind_curr, 2), y = mn_curr + c(-2, 2) * se_curr)
    lines(x = rep(ind_curr, 2) + c(-1, 1) * jit_eps1, y = mn_curr + c(2, 2) * se_curr)
    lines(x = rep(ind_curr, 2) + c(-1, 1) * jit_eps1, y = mn_curr + c(-2, -2) * se_curr)
  }
  if (match(assay_curr, assay_unique) == 1) {
    ats <- c(.1, .2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
    axis(side = 2, at = log2(ats), labels = ats, las = 1, cex.axis = 1.1)
  }
  n_in_grp <- sapply(geno_bio_unique, function(x) sum(d_sample_curr$geno_bio == x))
  # mtext(side = 1, at = 1:3, text = paste0(gsub("_", " ", geno_bio_unique), " (n = ", n_in_grp, ")"), las = 2, line = 2)
  mtext(side = 1, at = 1:3, text = paste0(gsub("_", " ", geno_bio_unique)), las = 2, line = 2)
  for (j in 1:nrow(test_df)) {
    if (nonpar_qval[assay_curr, test_df$test[j]] < .05) {
      ind1 <- match(test_df$geno_bio_1[j], geno_bio_unique)
      ind2 <- match(test_df$geno_bio_2[j], geno_bio_unique)
      yat <- ylim_log[2] - j * .75
      pvalc <- nonpar_pval[assay_curr, test_df$test[j]]
      qvalc <- nonpar_qval[assay_curr, test_df$test[j]]
      starlab <- ""
      if (pvalc < .0001) {
        starlab <- "****"
      } else if (pvalc < .001) {
        starlab <- "***"
      } else if (pvalc < .01) {
        starlab <- "**"
      } else if (pvalc < .05) {
        starlab <- "*"
      } 
      if (starlab != "") {
        lines(x = c(ind1, ind2), y = rep(yat, 2))
        lines(x = rep(ind1, 2), y = yat + c(0, -.1))
        lines(x = rep(ind2, 2), y = yat + c(0, -.1))
        text(x = mean(c(ind1, ind2)), y = yat + .25, labels = starlab, cex = 1.5, col = ifelse(qvalc < .05, 1, "grey"))
      }
    }
  }
  # abline(h = log2(range_of_detection[[assay_curr]][1]), lty = 2, lwd = 2)
}
mtext(outer = TRUE, side = 2, text = "pg/ml (log2 scale)", line = 4)
pretty(ylim_log)
if (output_plot) {
  dev.off()
}

cbind(c(nonpar_pval), c(nonpar_qval))



null_t_list <- list()
par(mfrow = c(3, 3))
n_perms <- 1000
d_sample_plate_assay$y <- d_sample_plate_assay$mean_log_conc
for (perm in 0:n_perms) {
  try({
    print(perm)
    d_sample_perm <- d_sample
    if (perm >= 1) {
      set.seed(perm)
      d_sample_perm$geno_bio <- sample(d_sample$geno_bio, replace = F)
    }
    for (assay_curr in assay_unique) {
      d_curr <- d_sample_plate_assay[d_sample_plate_assay$Assay == assay_curr, ]
      d_curr$geno_bio <- d_sample_perm[match(d_curr$sample_plate, d_sample_perm$sample_plate), "geno_bio"]
      if (qn) {
        d_curr$log_conc <- norm_quantile_map(d_curr$log_conc)
      }
      lm_out <- lm(y ~ geno_bio + Plate, data = d_curr)
      lm_out_no_intercept <- lm(y ~ -1 + geno_bio + Plate, data = d_curr)
      # qqnorm(lm_out$residuals)
      # qqline(lm_out$residuals)
      d_curr_releveled <- d_curr
      d_curr_releveled$geno_bio <- relevel(d_curr_releveled$geno_bio, ref = "Dp5_blood")
      lm_out_dp5_blood_intercept <- lm(y ~ geno_bio + Plate, data = d_curr_releveled)
      # nlme::lme(fixed = y ~ geno_bio + Plate, 
      #                      random = list(~ 1 | sample_name, ~ 1 | litter),
      #                      data = d_curr_releveled)
      # lme_out_no_intercept <- nlme::lme(fixed = y ~ -1 + geno_bio + Plate, 
      #                                          random = list(~ 1 | sample_name, ~ 1 | litter),
      #                                          data = d_curr_releveled)
      ttab_curr <- summary(lm_out)$coef
      ttab_curr_dp5_blood_intercept <- summary(lm_out_dp5_blood_intercept)$coef
      ttab_curr_no_intercept <- summary(lm_out_no_intercept)$coef
      for (t_out_name in t_out_names) {
        t_list[[t_out_name]][assay_curr, "Dp5_blood_v_WT_blood"] <- ttab_curr["geno_bioDp5_blood", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_fluid_v_WT_blood"] <- ttab_curr["geno_bioDp5_fluid", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_fluid_v_Dp5_blood"] <- ttab_curr_dp5_blood_intercept["geno_bioDp5_fluid", t_out_name]
        t_list[[t_out_name]][assay_curr, "WT_blood"] <- ttab_curr_no_intercept["geno_bioWT_blood", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_blood"] <- ttab_curr_no_intercept["geno_bioDp5_blood", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_fluid"] <- ttab_curr_no_intercept["geno_bioDp5_fluid", t_out_name]
      }
      if (perm >= 1) {
        null_t_list[[perm]] <- t_list
      } else {
        true_t <- t_list
      }
      # anova_curr <- anova(lme_out, Terms = "geno_bio")
      # anova_pval[assay_curr] <- anova_curr$`p-value`
    }
  })
}




sapply(null_t_list, function(x) x$'Pr(>|t|)')


t_list$'q-value' <- t_list$'p-value'
t_list$'q-value'[] <- NA
t_list$'q-value'[, 1:3] <- p.adjust(t_list$'p-value'[, 1:3], method = "BH")

t_list$'q-value'
t_list$'p-value'











####


null_t_list <- list()
par(mfrow = c(3, 3))
n_perms <- 100
for (perm in 0:n_perms) {
  try({
    print(perm)
    d_sample_perm <- d_sample
    if (perm >= 1) {
      set.seed(perm)
      d_sample_perm$geno_bio <- sample(d_sample$geno_bio, replace = F)
    }
    for (assay_curr in assay_unique) {
      d_curr <- d[d$Assay == assay_curr, ]
      d_curr$geno_bio <- d_sample_perm[match(d_curr$sample_name, d_sample_perm$sample_name), "geno_bio"]
      if (qn) {
        d_curr$log_conc <- norm_quantile_map(d_curr$log_conc)
      }
      lme_out <- nlme::lme(fixed = y ~ geno_bio + Plate,
                           random = list(~ 1 | sample_name, ~ 1 | litter),
                           data = d_curr)
      # qqnorm(lme_out$residuals)
      # qqline(lme_out$residuals)
      d_curr_releveled <- d_curr
      d_curr_releveled$geno_bio <- relevel(d_curr_releveled$geno_bio, ref = "Dp5_blood")
      lme_out_dp5_blood_intercept <- nlme::lme(fixed = y ~ geno_bio + Plate, 
                                               random = list(~ 1 | sample_name, ~ 1 | litter),
                                               data = d_curr_releveled)
      lme_out_no_intercept <- nlme::lme(fixed = y ~ -1 + geno_bio + Plate, 
                                        random = list(~ 1 | sample_name, ~ 1 | litter),
                                        data = d_curr_releveled)
      ttab_curr <- summary(lme_out)$tTable
      ttab_curr_dp5_blood_intercept <- summary(lme_out_dp5_blood_intercept)$tTable
      ttab_curr_no_intercept <- summary(lme_out_no_intercept)$tTable
      for (t_out_name in t_out_names) {
        t_list[[t_out_name]][assay_curr, "Dp5_blood_v_WT_blood"] <- ttab_curr["geno_bioDp5_blood", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_fluid_v_WT_blood"] <- ttab_curr["geno_bioDp5_fluid", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_fluid_v_Dp5_blood"] <- ttab_curr_dp5_blood_intercept["geno_bioDp5_fluid", t_out_name]
        t_list[[t_out_name]][assay_curr, "WT_blood"] <- ttab_curr_no_intercept["geno_bioWT_blood", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_blood"] <- ttab_curr_no_intercept["geno_bioDp5_blood", t_out_name]
        t_list[[t_out_name]][assay_curr, "Dp5_fluid"] <- ttab_curr_no_intercept["geno_bioDp5_fluid", t_out_name]
      }
      if (perm >= 1) {
        null_t_list[[perm]] <- t_list
      } else {
        true_t <- t_list
      }
      anova_curr <- anova(lme_out, Terms = "geno_bio")
      anova_pval[assay_curr] <- anova_curr$`p-value`
    }
  })
}

