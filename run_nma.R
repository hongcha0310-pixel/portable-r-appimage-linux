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
           cex.number.of.studies=1.0, plastic=FALSE,
           start.layout='circle', lwd=3)
  dev.off()
  pdf(file.path(dir_out, paste0(outcome, '_network.pdf')), width=8, height=6)
  netgraph(nma, seq='optimal', points=TRUE, cex.points=4.2, pch.points=21,
           col.points='black', bg.points='white', cex=1.25,
           thickness='number.of.studies', number.of.studies=TRUE,
           cex.number.of.studies=1.0, plastic=FALSE,
           start.layout='circle', lwd=3)
  dev.off()

  png(file.path(dir_out, paste0(outcome, '_forest.png')), width=2400, height=1800, res=300, bg='white')
  forest(nma, reference.group='IC', pooled='random',
         rightcols=c('effect','ci'),
         xlab='Standardized mean difference (SMD)',
         smlab='Random-effects network meta-analysis')
  dev.off()
  pdf(file.path(dir_out, paste0(outcome, '_forest.pdf')), width=8.5, height=6.5)
  forest(nma, reference.group='IC', pooled='random',
         rightcols=c('effect','ci'),
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
