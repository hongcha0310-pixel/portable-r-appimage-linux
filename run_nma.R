options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(meta)
  library(netmeta)
})

outdir <- file.path('analysis', 'results')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

capture.output(sessionInfo(), file = file.path(outdir, 'sessionInfo.txt'))
writeLines(c(
  paste0('R version: ', R.version.string),
  paste0('meta version: ', as.character(packageVersion('meta'))),
  paste0('netmeta version: ', as.character(packageVersion('netmeta')))
), file.path(outdir, 'versions.txt'))

# Higher positive SMD values consistently indicate a more favorable outcome.
# Arm-based continuous outcomes are transformed with meta::pairwise(sm = "SMD"),
# which uses the multi-arm-consistent Cohen d formulation.

convert_direct <- function(df) {
  if (nrow(df) == 0) {
    return(data.frame(TE=numeric(), seTE=numeric(), studlab=character(),
                      treat1=character(), treat2=character(), source=character()))
  }
  out <- vector('list', nrow(df))
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    if (r$type == 'binary_or_to_smd') {
      a <- as.numeric(r$a); b <- as.numeric(r$b); c <- as.numeric(r$c); d <- as.numeric(r$d)
      logOR <- log((a * d) / (b * c))
      seLogOR <- sqrt(1/a + 1/b + 1/c + 1/d)
      conv <- sqrt(3) / pi
      TE <- logOR * conv
      seTE <- seLogOR * conv
    } else if (r$type == 'md_ci_to_smd') {
      n1 <- as.numeric(r$n1); n2 <- as.numeric(r$n2)
      md <- as.numeric(r$md); lower <- as.numeric(r$lower); upper <- as.numeric(r$upper)
      dfree <- n1 + n2 - 2
      tcrit <- qt(0.975, df = dfree)
      se_md <- (upper - lower) / (2 * tcrit)
      sd_pooled <- se_md / sqrt(1/n1 + 1/n2)
      dval <- md / sd_pooled
      if (r$direction == 'lower_better') dval <- -dval
      TE <- dval
      seTE <- sqrt((n1+n2)/(n1*n2) + (dval^2)/(2*dfree))
    } else {
      stop('Unknown direct contrast type: ', r$type)
    }
    out[[i]] <- data.frame(TE=TE, seTE=seTE, studlab=r$study,
                           treat1=r$treat1, treat2=r$treat2,
                           source=r$source, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

keep_ic_component <- function(dat) {
  trts <- unique(c(dat$treat1, dat$treat2))
  if (!'IC' %in% trts) return(dat[0, , drop=FALSE])
  seen <- 'IC'
  repeat {
    before <- length(seen)
    hit <- dat$treat1 %in% seen | dat$treat2 %in% seen
    seen <- unique(c(seen, dat$treat1[hit], dat$treat2[hit]))
    if (length(seen) == before) break
  }
  dat[dat$treat1 %in% seen & dat$treat2 %in% seen, , drop=FALSE]
}

fit_nma <- function(dat) {
  dat <- keep_ic_component(dat)
  netmeta(TE, seTE, treat1, treat2, studlab,
          data = dat,
          sm = 'SMD',
          common = FALSE,
          random = TRUE,
          method.tau = 'REML',
          reference.group = 'IC',
          baseline.reference = TRUE,
          small.values = 'undesirable',
          all.treatments = TRUE,
          warn = TRUE)
}

extract_vs_ic <- function(nma, scenario='Primary') {
  trts <- setdiff(rownames(nma$TE.random), 'IC')
  if (length(trts) == 0) return(data.frame())
  # Preserve both orientations so the extraction can be audited directly.
  data.frame(
    scenario = scenario,
    treatment = trts,
    TE_treatment_vs_IC = nma$TE.random[trts, 'IC'],
    lower_treatment_vs_IC = nma$lower.random[trts, 'IC'],
    upper_treatment_vs_IC = nma$upper.random[trts, 'IC'],
    p_treatment_vs_IC = nma$pval.random[trts, 'IC'],
    TE_IC_vs_treatment = nma$TE.random['IC', trts],
    lower_IC_vs_treatment = nma$lower.random['IC', trts],
    upper_IC_vs_treatment = nma$upper.random['IC', trts],
    p_IC_vs_treatment = nma$pval.random['IC', trts],
    row.names = NULL,
    check.names = FALSE
  )
}

make_league <- function(nma) {
  trts <- rownames(nma$TE.random)
  m <- matrix('', nrow=length(trts), ncol=length(trts), dimnames=list(trts,trts))
  for (i in seq_along(trts)) for (j in seq_along(trts)) {
    if (i == j) m[i,j] <- trts[i]
    else m[i,j] <- sprintf('%.3f [%.3f, %.3f]', nma$TE.random[i,j], nma$lower.random[i,j], nma$upper.random[i,j])
  }
  as.data.frame(m, check.names=FALSE)
}

save_plots <- function(nma, outcome, dir_out) {
  png(file.path(dir_out, paste0(outcome, '_network.png')), width=2400, height=1800, res=300, bg='white')
  netgraph(nma, seq='optimal', points=TRUE, cex.points=4.2, pch.points=21,
           col.points='black', bg.points='white', cex=1.25,
           thickness='number.of.studies', number.of.studies=TRUE,
           cex.number.of.studies=1.0, plastic=FALSE, iterate=TRUE,
           start.layout='circle', lwd=3)
  dev.off()
  pdf(file.path(dir_out, paste0(outcome, '_network.pdf')), width=8, height=6)
  netgraph(nma, seq='optimal', points=TRUE, cex.points=4.2, pch.points=21,
           col.points='black', bg.points='white', cex=1.25,
           thickness='number.of.studies', number.of.studies=TRUE,
           cex.number.of.studies=1.0, plastic=FALSE, iterate=TRUE,
           start.layout='circle', lwd=3)
  dev.off()

  png(file.path(dir_out, paste0(outcome, '_forest.png')), width=2400, height=1800, res=300, bg='white')
  forest(nma, reference='IC', common=FALSE, random=TRUE,
         rightcols=c('effect','ci','Pscore'),
         xlab='Standardized mean difference (SMD)',
         smlab='Random-effects network meta-analysis')
  dev.off()
  pdf(file.path(dir_out, paste0(outcome, '_forest.pdf')), width=8.5, height=6.5)
  forest(nma, reference='IC', common=FALSE, random=TRUE,
         rightcols=c('effect','ci','Pscore'),
         xlab='Standardized mean difference (SMD)',
         smlab='Random-effects network meta-analysis')
  dev.off()

  rk <- netrank(nma, method='P-score', small.values='undesirable')
  ps <- sort(rk$ranking.random, decreasing=TRUE)
  png(file.path(dir_out, paste0(outcome, '_Pscore.png')), width=2200, height=1500, res=300, bg='white')
  par(mar=c(8,5,2,1))
  barplot(ps, las=2, ylim=c(0,1), ylab='P-score', main=paste(outcome, 'P-score ranking'))
  abline(h=c(0,0.5,1), lty=c(1,3,1), col='grey70')
  dev.off()
  pdf(file.path(dir_out, paste0(outcome, '_Pscore.pdf')), width=8, height=5.5)
  par(mar=c(8,5,2,1))
  barplot(ps, las=2, ylim=c(0,1), ylab='P-score', main=paste(outcome, 'P-score ranking'))
  abline(h=c(0,0.5,1), lty=c(1,3,1), col='grey70')
  dev.off()
}

analyse_outcome <- function(outcome, arms, direct_raw, high_rob_studies) {
  odir <- file.path(outdir, outcome)
  dir.create(odir, recursive=TRUE, showWarnings=FALSE)

  pw <- meta::pairwise(studlab = study, treat = treatment,
                       n = n, mean = mean, sd = sd,
                       data = arms, sm = 'SMD',
                       keep.all.comparisons = TRUE,
                       append = FALSE)
  arm_contr <- as.data.frame(pw[, c('TE','seTE','studlab','treat1','treat2')])
  arm_contr$source <- 'Arm-based; meta::pairwise(sm="SMD")'
  direct <- convert_direct(direct_raw)
  dat <- rbind(arm_contr, direct)
  dat <- dat[is.finite(dat$TE) & is.finite(dat$seTE) & dat$seTE > 0, , drop=FALSE]
  write.csv(dat, file.path(odir, paste0(outcome, '_contrasts_used.csv')), row.names=FALSE)

  nma <- fit_nma(dat)
  saveRDS(nma, file.path(odir, paste0(outcome, '_netmeta_model.rds')))
  capture.output(print(nma, digits=4), file=file.path(odir, paste0(outcome, '_model_print.txt')))
  capture.output(summary(nma), file=file.path(odir, paste0(outcome, '_model_summary.txt')))

  # Core statistics
  I2raw <- as.numeric(nma$I2)
  I2pct <- ifelse(is.na(I2raw), NA_real_, ifelse(I2raw <= 1, I2raw*100, I2raw))
  core <- data.frame(
    outcome=outcome,
    studies=length(unique(dat$studlab)),
    treatments=length(unique(c(dat$treat1, dat$treat2))),
    tau=as.numeric(nma$tau),
    tau2=as.numeric(nma$tau)^2,
    Q=as.numeric(nma$Q),
    df_Q=as.numeric(nma$df.Q),
    p_Q=as.numeric(nma$pval.Q),
    I2_percent=I2pct
  )
  write.csv(core, file.path(odir, paste0(outcome, '_core_stats.csv')), row.names=FALSE)

  # Reference estimates and matrices
  vic <- extract_vs_ic(nma, 'Primary')
  write.csv(vic, file.path(odir, paste0(outcome, '_vs_IC.csv')), row.names=FALSE)
  write.csv(data.frame(treatment=names(netrank(nma, method='P-score', small.values='undesirable')$ranking.random),
                       Pscore=as.numeric(netrank(nma, method='P-score', small.values='undesirable')$ranking.random)),
            file.path(odir, paste0(outcome, '_Pscore.csv')), row.names=FALSE)
  write.csv(cbind(treatment=rownames(nma$TE.random), as.data.frame(nma$TE.random, check.names=FALSE)),
            file.path(odir, paste0(outcome, '_TE_random_matrix.csv')), row.names=FALSE)
  write.csv(cbind(treatment=rownames(nma$lower.random), as.data.frame(nma$lower.random, check.names=FALSE)),
            file.path(odir, paste0(outcome, '_lower_random_matrix.csv')), row.names=FALSE)
  write.csv(cbind(treatment=rownames(nma$upper.random), as.data.frame(nma$upper.random, check.names=FALSE)),
            file.path(odir, paste0(outcome, '_upper_random_matrix.csv')), row.names=FALSE)
  league <- make_league(nma)
  write.csv(cbind(treatment=rownames(league), league), file.path(odir, paste0(outcome, '_league.csv')), row.names=FALSE)

  # Network geometry data from netgraph
  ng <- netgraph(nma, figure=FALSE, number.of.studies=TRUE, seq='optimal')
  write.csv(ng$nodes, file.path(odir, paste0(outcome, '_network_nodes.csv')), row.names=FALSE)
  write.csv(ng$edges, file.path(odir, paste0(outcome, '_network_edges.csv')), row.names=FALSE)

  # Global and local inconsistency
  dd <- tryCatch(decomp.design(nma), error=function(e) e)
  capture.output(print(dd), file=file.path(odir, paste0(outcome, '_design_decomposition.txt')))
  if (is.list(dd) && !is.null(dd$Q.decomp)) write.csv(dd$Q.decomp, file.path(odir, paste0(outcome, '_design_Q.csv')), row.names=FALSE)
  if (is.list(dd) && !is.null(dd$Q.inc.random)) write.csv(dd$Q.inc.random, file.path(odir, paste0(outcome, '_design_Q_random.csv')), row.names=FALSE)
  ns <- tryCatch(netsplit(nma, method='Back-calculation', common=FALSE, random=TRUE, ci=TRUE, test=TRUE), error=function(e) e)
  capture.output(print(ns), file=file.path(odir, paste0(outcome, '_netsplit.txt')))
  capture.output(str(ns), file=file.path(odir, paste0(outcome, '_netsplit_structure.txt')))

  # Sensitivity 1: exclude high risk of bias studies
  dat_hr <- dat[!dat$studlab %in% high_rob_studies, , drop=FALSE]
  dat_hr <- keep_ic_component(dat_hr)
  if (nrow(dat_hr) > 0 && length(unique(c(dat_hr$treat1,dat_hr$treat2))) >= 2) {
    nma_hr <- fit_nma(dat_hr)
    write.csv(extract_vs_ic(nma_hr, 'Exclude high RoB'), file.path(odir, paste0(outcome, '_sensitivity_highRoB.csv')), row.names=FALSE)
    saveRDS(nma_hr, file.path(odir, paste0(outcome, '_sensitivity_highRoB.rds')))
  }

  # Sensitivity 2: exclude crossover trial
  dat_x <- dat[dat$studlab != 'Truque-Díaz 2025, crossover', , drop=FALSE]
  dat_x <- keep_ic_component(dat_x)
  if (nrow(dat_x) > 0 && length(unique(c(dat_x$treat1,dat_x$treat2))) >= 2) {
    nma_x <- fit_nma(dat_x)
    write.csv(extract_vs_ic(nma_x, 'Exclude crossover'), file.path(odir, paste0(outcome, '_sensitivity_crossover.csv')), row.names=FALSE)
    saveRDS(nma_x, file.path(odir, paste0(outcome, '_sensitivity_crossover.rds')))
  }

  # Leave-one-study-out analysis for all reference comparisons
  loo <- list()
  for (st in sort(unique(dat$studlab))) {
    d <- dat[dat$studlab != st, , drop=FALSE]
    d <- keep_ic_component(d)
    if (nrow(d) == 0 || length(unique(c(d$treat1,d$treat2))) < 2) next
    fit <- tryCatch(fit_nma(d), error=function(e) NULL)
    if (is.null(fit)) next
    z <- extract_vs_ic(fit, paste0('Exclude: ', st))
    z$excluded_study <- st
    loo[[length(loo)+1]] <- z
  }
  if (length(loo)) write.csv(do.call(rbind, loo), file.path(odir, paste0(outcome, '_leave_one_out.csv')), row.names=FALSE)

  save_plots(nma, outcome, odir)
  invisible(list(nma=nma, core=core, vs_ic=vic))
}


# Embedded analysis input data derived from the locked extraction workbook.
strength_arms <- read.csv(textConnection(r"x(study,treatment,n,mean,sd,source
Calatayud 2020 / Cruz-Montecinos 2023,EX,10,12.72,7.143287437850708,"Calatayud 2020, Table 3"
Calatayud 2020 / Cruz-Montecinos 2023,IC,10,-5.08,7.597606110512448,"Calatayud 2020, Table 3"
Cuesta-Barriuso 2014a,MAN,11,-0.05,0.153,"Cuesta-Barriuso 2014, Table 3"
Cuesta-Barriuso 2014a,EX,10,-0.053000000000000005,0.157,"Cuesta-Barriuso 2014, Table 3"
Cuesta-Barriuso 2014a,IC,10,0,0,"Cuesta-Barriuso 2014, Table 3"
Donoso-Úbeda 2023,MAN,26,280.83,77.2,"Donoso-Úbeda 2023, Table II"
Donoso-Úbeda 2023,IC,26,271.33,26.26,"Donoso-Úbeda 2023, Table II"
Donoso-Úbeda 2025,MAN,29,26.29,34.39,"Donoso-Úbeda 2024/2025, Table 2"
Donoso-Úbeda 2025,IC,29,-2.29,4.78,"Donoso-Úbeda 2024/2025, Table 2"
Goto 2014,DIG,16,1.44,0.7,"Goto 2014, Table 2"
Goto 2014,EX,16,1.37,0.7,"Goto 2014, Table 2"
Parhampour 2022,EX,28,7.299999999999999,0.958772363204345,"Parhampour 2022, Table 3"
Parhampour 2022,IC,14,5.69,0.5,"Parhampour 2022, Table 3"
Runkel 2016 / Runkel 2017,EX,24,5.05,1.9,"Runkel 2016, Figure 2b"
Runkel 2016 / Runkel 2017,IC,28,4.8,2,"Runkel 2016, Figure 2b"
Tat 2021,MULTI,8,17.1,2.5925925925925926,"Tat 2021, Table 3"
Tat 2021,EX,9,16.8,7.037037037037036,"Tat 2021, Table 3"
Villalón-González 2025/2026,EX,16,484.63,197.93,"Villalón-González 2025/2026, Table 3"
Villalón-González 2025/2026,IC,17,430.42,152.19,"Villalón-González 2025/2026, Table 3")x"), check.names=FALSE, stringsAsFactors=FALSE)

strength_direct <- read.csv(textConnection(r"x(study,type,treat1,treat2,a,b,c,d,n1,n2,md,lower,upper,direction,source)x"), check.names=FALSE, stringsAsFactors=FALSE)

pain_arms <- read.csv(textConnection(r"x(study,treatment,n,mean,sd,source
Cuesta-Barriuso 2014a,MAN,11,-1.1,1.818,"Cuesta-Barriuso 2014, Table 3"
Cuesta-Barriuso 2014a,EX,10,-0.44699999999999995,0.831,"Cuesta-Barriuso 2014, Table 3"
Cuesta-Barriuso 2014a,IC,10,-0.441,0.682,"Cuesta-Barriuso 2014, Table 3"
Cuesta-Barriuso 2017,EX,10,-0.32,0.79,"Cuesta-Barriuso 2017, Table 3"
Cuesta-Barriuso 2017,IC,10,-0.3,0.56,"Cuesta-Barriuso 2017, Table 3"
Cuesta-Barriuso 2020 / Donoso-Úbeda 2020,MAN,33,-1.74,1.94,"Donoso-Úbeda 2020, Table 3"
Cuesta-Barriuso 2020 / Donoso-Úbeda 2020,IC,32,-4.08,2.92,"Donoso-Úbeda 2020, Table 3"
"Cuesta-Barriuso 2021, knee",MAN,14,-2,1.2600000000000002,"Cuesta-Barriuso 2021 knee, Table 3"
"Cuesta-Barriuso 2021, knee",IC,14,-4.86,2.74,"Cuesta-Barriuso 2021 knee, Table 3"
Cuesta-Barriuso elbow family,MAN,35,-0.6,1.3,"Cuesta-Barriuso 2022, Table 2"
Cuesta-Barriuso elbow family,IC,34,-2.8,2.2,"Cuesta-Barriuso 2022, Table 2"
Donoso-Úbeda 2023,MAN,26,-1.8599999999999999,2.05,"Donoso-Úbeda 2023, Table II"
Donoso-Úbeda 2023,IC,26,-2.09,1.32,"Donoso-Úbeda 2023, Table II"
Donoso-Úbeda 2025,MAN,29,0.68,0.88,"Donoso-Úbeda 2024/2025, Table 2"
Donoso-Úbeda 2025,IC,29,0.4,1.8,"Donoso-Úbeda 2024/2025, Table 2"
Goto 2014,DIG,16,-16.7,15.3,"Goto 2014, Table 2"
Goto 2014,EX,16,-5.1,8.8,"Goto 2014, Table 2"
Khami 2020,EPM,20,-1.33,1.44,"Khami 2020, Table 3"
Khami 2020,IC,16,-5.19,2.71,"Khami 2020, Table 3"
Pires Oliveira 2020,ACU,15,-2.23,2.54,"Pires Oliveira 2020, Table 2"
Pires Oliveira 2020,EPM,13,-5.1899999999999995,2.48,"Pires Oliveira 2020, Table 2"
Tat 2021,MULTI,8,0,1.4814814814814814,"Tat 2021, Table 3"
Tat 2021,EX,9,-4,2.222222222222222,"Tat 2021, Table 3"
Tat 2026,MULTI,18,-2.31,2.12,"Tat 2026, Table 5"
Tat 2026,IC,17,-3.8799999999999994,2.23,"Tat 2026, Tables 4–5"
"Truque-Díaz 2025, pilot",MAN,11,0,3.5555555555555554,"Truque-Díaz 2025 pilot, Table 3"
"Truque-Díaz 2025, pilot",IC,13,-1.5,2.9999999999999996,"Truque-Díaz 2025 pilot, Table 3"
Ucero-Lozano 2024,IMV,24,-2.34,1.474578793421815,"Ucero-Lozano 2024, Table III"
Ucero-Lozano 2024,IC,24,-3.43,1.6054585679858813,"Ucero-Lozano 2024, Table III"
Ucero-Lozano 2025,IMV,14,-1.08,0.8757061434319141,"Ucero-Lozano 2025, Table 2"
Ucero-Lozano 2025,IC,14,-1.85,0.8260206884854225,"Ucero-Lozano 2025, Table 2")x"), check.names=FALSE, stringsAsFactors=FALSE)

pain_direct <- read.csv(textConnection(r"x(study,type,treat1,treat2,a,b,c,d,n1,n2,md,lower,upper,direction,source
Calatayud 2020 / Cruz-Montecinos 2023,binary_or_to_smd,EX,IC,7,3,2,8,10,10,,,,higher_better,Calatayud 2020 clinically important pain improvement
"Truque-Díaz 2025, crossover",md_ci_to_smd,MAN,IC,,,,,6,7,-0.85,-2.01,0.29,lower_better,Truque-Díaz 2025 crossover phase 1)x"), check.names=FALSE, stringsAsFactors=FALSE)

hjhs_arms <- read.csv(textConnection(r"x(study,treatment,n,mean,sd,source
Calatayud 2020 / Cruz-Montecinos 2023,EX,10,2.89,4.13779468024229,"Calatayud 2020, Table 3"
Calatayud 2020 / Cruz-Montecinos 2023,IC,10,-0.99,4.13779468024229,"Calatayud 2020, Table 3"
Cuesta-Barriuso 2020 / Donoso-Úbeda 2020,MAN,33,-5.45,3.07,"Donoso-Úbeda 2020, Table 3"
Cuesta-Barriuso 2020 / Donoso-Úbeda 2020,IC,32,-7.03,3.27,"Donoso-Úbeda 2020, Table 3"
"Cuesta-Barriuso 2021, knee",MAN,14,-8,2.18,"Cuesta-Barriuso 2021 knee, Table 3"
"Cuesta-Barriuso 2021, knee",IC,14,-9.429999999999998,2.13,"Cuesta-Barriuso 2021 knee, Table 3"
Cuesta-Barriuso elbow family,MAN,35,-8.4,3.6,"Cuesta-Barriuso 2022, Table 2"
Cuesta-Barriuso elbow family,IC,34,-9.9,4.4,"Cuesta-Barriuso 2022, Table 2"
Feng 2026,DIG,8,-31.63,13.14,Feng 2026
Feng 2026,EX,9,-22.11,8.73,Feng 2026
Khami 2020,EPM,20,-9.65,3.15,"Khami 2020, Table 3"
Khami 2020,IC,16,-14.19,4.52,"Khami 2020, Table 3"
Makkar 2025,MULTI,16,-22.5,7.5,"Makkar 2025, Table 3"
Makkar 2025,EX,15,-21.76,7.9,"Makkar 2025, Table 3"
Parhampour 2019,EX,36,-4.760000000000001,1.1191872050733962,"Parhampour 2019, Table 3"
Parhampour 2019,IC,12,-9.08,1.67,"Parhampour 2019, Table 3"
Pires Oliveira 2020,ACU,15,-64.07,17.3,"Pires Oliveira 2020, Table 2"
Pires Oliveira 2020,EPM,13,-60.92,17.28,"Pires Oliveira 2020, Table 2"
Tat 2021,MULTI,8,7,3.7037037037037033,"Tat 2021, Table 5"
Tat 2021,EX,9,0,2.222222222222222,"Tat 2021, Table 5"
Tat 2026,MULTI,18,-8.17,9.48,"Tat 2026, Table 5"
Tat 2026,IC,17,-14.41,11.3,"Tat 2026, Table 4"
"Truque-Díaz 2025, pilot",MAN,11,-7,5.185185185185185,"Truque-Díaz 2025 pilot, Table 3"
"Truque-Díaz 2025, pilot",IC,13,-7,3.333333333333333,"Truque-Díaz 2025 pilot, Table 3")x"), check.names=FALSE, stringsAsFactors=FALSE)

hjhs_direct <- read.csv(textConnection(r"x(study,type,treat1,treat2,a,b,c,d,n1,n2,md,lower,upper,direction,source
"Truque-Díaz 2025, crossover",md_ci_to_smd,MAN,IC,,,,,6,7,-2.39,-2.84,-1.94,lower_better,Truque-Díaz 2025 crossover phase 1)x"), check.names=FALSE, stringsAsFactors=FALSE)


high_rob <- c(
  'Goto 2014',
  'Khami 2020',
  'Pires Oliveira 2020',
  'Runkel 2016 / Runkel 2017',
  'Tat 2021'
)

res_strength <- analyse_outcome('Strength', strength_arms, strength_direct, high_rob)
res_pain     <- analyse_outcome('Pain', pain_arms, pain_direct, high_rob)
res_hjhs     <- analyse_outcome('HJHS', hjhs_arms, hjhs_direct, high_rob)

summary_core <- rbind(res_strength$core, res_pain$core, res_hjhs$core)
write.csv(summary_core, file.path(outdir, 'summary_core_stats.csv'), row.names=FALSE)

all_vs_ic <- rbind(
  transform(res_strength$vs_ic, outcome='Strength'),
  transform(res_pain$vs_ic, outcome='Pain'),
  transform(res_hjhs$vs_ic, outcome='HJHS')
)
write.csv(all_vs_ic, file.path(outdir, 'all_primary_vs_IC.csv'), row.names=FALSE)

cat('R NMA completed successfully.\n')
print(summary_core)
print(all_vs_ic)
