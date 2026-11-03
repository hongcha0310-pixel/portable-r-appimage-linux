options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(meta)
  library(netmeta)
})

outdir <- file.path("analysis", "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

writeLines(c(
  paste0("R version: ", R.version.string),
  paste0("meta version: ", as.character(packageVersion("meta"))),
  paste0("netmeta version: ", as.character(packageVersion("netmeta")))
), file.path(outdir, "versions.txt"))
capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))

# All arm means are harmonized before analysis so higher values indicate better outcomes.
# Therefore, a positive SMD always favors treat1 over treat2.

read_embedded_csv <- function(txt) {
  read.csv(text = txt, check.names = FALSE, stringsAsFactors = FALSE,
           na.strings = c("", "NA"))
}

convert_direct <- function(df) {
  if (nrow(df) == 0) {
    return(data.frame(TE=numeric(), seTE=numeric(), studlab=character(),
                      treat1=character(), treat2=character(), source=character(),
                      stringsAsFactors=FALSE))
  }
  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    rr <- df[i, , drop=FALSE]
    if (rr$type == "binary_or_to_smd") {
      a <- as.numeric(rr$a); b <- as.numeric(rr$b)
      c0 <- as.numeric(rr$c); d <- as.numeric(rr$d)
      log_or <- log((a * d) / (b * c0))
      se_log_or <- sqrt(1/a + 1/b + 1/c0 + 1/d)
      hh <- sqrt(3) / pi
      TE <- log_or * hh
      seTE <- se_log_or * hh
    } else if (rr$type == "md_ci_to_smd") {
      n1 <- as.numeric(rr$n1); n2 <- as.numeric(rr$n2)
      md <- as.numeric(rr$md); lower <- as.numeric(rr$lower); upper <- as.numeric(rr$upper)
      dfree <- n1 + n2 - 2
      tcrit <- qt(0.975, df = dfree)
      se_md <- (upper - lower) / (2 * tcrit)
      sd_pooled <- se_md / sqrt(1/n1 + 1/n2)
      dval <- md / sd_pooled
      if (rr$direction == "lower_better") dval <- -dval
      TE <- dval
      seTE <- sqrt((n1+n2)/(n1*n2) + (dval^2)/(2*dfree))
    } else {
      stop("Unknown direct contrast type: ", rr$type)
    }
    out[[i]] <- data.frame(
      TE=TE, seTE=seTE, studlab=as.character(rr$study),
      treat1=as.character(rr$treat1), treat2=as.character(rr$treat2),
      source=as.character(rr$source), stringsAsFactors=FALSE
    )
  }
  do.call(rbind, out)
}

keep_ic_component <- function(dat) {
  if (nrow(dat) == 0) return(dat)
  trts <- unique(c(dat$treat1, dat$treat2))
  if (!"IC" %in% trts) return(dat[0, , drop=FALSE])
  seen <- "IC"
  repeat {
    hit <- dat$treat1 %in% seen | dat$treat2 %in% seen
    new_seen <- unique(c(seen, dat$treat1[hit], dat$treat2[hit]))
    if (length(new_seen) == length(seen)) break
    seen <- new_seen
  }
  dat[dat$treat1 %in% seen & dat$treat2 %in% seen, , drop=FALSE]
}

fit_nma <- function(dat) {
  dat <- keep_ic_component(dat)
  netmeta(
    TE, seTE, treat1, treat2, studlab,
    data = dat, sm = "SMD",
    common = FALSE, random = TRUE,
    method.tau = "REML",
    reference.group = "IC",
    small.values = "undesirable",
    warn = TRUE
  )
}

extract_vs_ic <- function(nma, scenario="Primary") {
  trts <- setdiff(rownames(nma$TE.random), "IC")
  if (!length(trts)) return(data.frame())
  data.frame(
    scenario=scenario,
    treatment=trts,
    SMD=nma$TE.random[trts, "IC"],
    lower=nma$lower.random[trts, "IC"],
    upper=nma$upper.random[trts, "IC"],
    p=nma$pval.random[trts, "IC"],
    row.names=NULL,
    check.names=FALSE
  )
}

core_stats <- function(nma, outcome, dat) {
  I2raw <- suppressWarnings(as.numeric(nma$I2))
  I2pct <- ifelse(is.na(I2raw), NA_real_, ifelse(I2raw <= 1, 100*I2raw, I2raw))
  tau <- suppressWarnings(as.numeric(nma$tau))
  data.frame(
    outcome=outcome,
    studies=length(unique(dat$studlab)),
    treatments=length(unique(c(dat$treat1, dat$treat2))),
    tau=tau,
    tau2=tau^2,
    Q=suppressWarnings(as.numeric(nma$Q)),
    df_Q=suppressWarnings(as.numeric(nma$df.Q)),
    p_Q=suppressWarnings(as.numeric(nma$pval.Q)),
    I2_percent=I2pct
  )
}

make_league <- function(nma) {
  trts <- rownames(nma$TE.random)
  m <- matrix("", nrow=length(trts), ncol=length(trts), dimnames=list(trts,trts))
  for (i in seq_along(trts)) {
    for (j in seq_along(trts)) {
      if (i == j) {
        m[i,j] <- trts[i]
      } else {
        m[i,j] <- sprintf("%.3f [%.3f, %.3f]", nma$TE.random[i,j], nma$lower.random[i,j], nma$upper.random[i,j])
      }
    }
  }
  as.data.frame(m, check.names=FALSE)
}

network_tables <- function(dat) {
  edge_key <- ifelse(dat$treat1 < dat$treat2,
                     paste(dat$treat1, dat$treat2, sep="|||"),
                     paste(dat$treat2, dat$treat1, sep="|||"))
  tmp <- data.frame(edge=edge_key, study=dat$studlab, stringsAsFactors=FALSE)
  split_edges <- strsplit(unique(tmp$edge), "\\|\\|\\|")
  edges <- do.call(rbind, lapply(seq_along(split_edges), function(i) {
    e <- unique(tmp$edge)[i]
    data.frame(treat1=split_edges[[i]][1], treat2=split_edges[[i]][2],
               studies=length(unique(tmp$study[tmp$edge == e])), stringsAsFactors=FALSE)
  }))
  trts <- sort(unique(c(dat$treat1, dat$treat2)))
  nodes <- do.call(rbind, lapply(trts, function(t) {
    st <- unique(c(dat$studlab[dat$treat1==t], dat$studlab[dat$treat2==t]))
    data.frame(treatment=t, studies=length(st), stringsAsFactors=FALSE)
  }))
  list(nodes=nodes, edges=edges)
}

safe_capture <- function(expr, file) {
  tryCatch({
    capture.output(force(expr), file=file)
    TRUE
  }, error=function(e) {
    writeLines(paste("ERROR:", conditionMessage(e)), file)
    FALSE
  })
}

save_optional_outputs <- function(nma, outcome, odir) {
  # Network plots
  for (ext in c("png", "pdf")) {
    outfile <- file.path(odir, paste0(outcome, "_network.", ext))
    tryCatch({
      if (ext == "png") png(outfile, width=2400, height=1800, res=300, bg="white") else pdf(outfile, width=8, height=6)
      netgraph(nma, seq="optimal", points=TRUE, cex.points=4.2, pch.points=21,
               col.points="black", bg.points="white", cex=1.25,
               thickness="number.of.studies", number.of.studies=TRUE,
               cex.number.of.studies=1.0, plastic=FALSE, start.layout="circle", lwd=3)
      dev.off()
    }, error=function(e) {
      try(dev.off(), silent=TRUE)
      writeLines(paste("Network plot error:", conditionMessage(e)), file.path(odir, paste0(outcome, "_network_plot_error.txt")))
    })
  }
  # Forest plots
  for (ext in c("png", "pdf")) {
    outfile <- file.path(odir, paste0(outcome, "_forest.", ext))
    tryCatch({
      if (ext == "png") png(outfile, width=2400, height=1800, res=300, bg="white") else pdf(outfile, width=8.5, height=6.5)
      forest(nma, reference.group="IC", pooled="random",
             rightcols=c("effect","ci"),
             xlab="Standardized mean difference (SMD)",
             smlab="Random-effects network meta-analysis")
      dev.off()
    }, error=function(e) {
      try(dev.off(), silent=TRUE)
      writeLines(paste("Forest plot error:", conditionMessage(e)), file.path(odir, paste0(outcome, "_forest_plot_error.txt")))
    })
  }
}

analyse_outcome <- function(outcome, arms, direct_raw, high_rob_studies) {
  message("--- ", outcome, " ---")
  odir <- file.path(outdir, outcome)
  dir.create(odir, recursive=TRUE, showWarnings=FALSE)

  pw <- meta::pairwise(
    studlab=study, treat=treatment,
    n=n, mean=mean, sd=sd,
    data=arms, sm="SMD"
  )
  arm_contr <- as.data.frame(pw[, c("TE","seTE","studlab","treat1","treat2")])
  arm_contr$source <- "Arm-based; meta::pairwise(sm='SMD')"
  direct <- convert_direct(direct_raw)
  dat <- rbind(arm_contr, direct)
  dat <- dat[is.finite(dat$TE) & is.finite(dat$seTE) & dat$seTE > 0, , drop=FALSE]
  write.csv(dat, file.path(odir, paste0(outcome, "_contrasts_used.csv")), row.names=FALSE)

  nma <- fit_nma(dat)
  saveRDS(nma, file.path(odir, paste0(outcome, "_netmeta_model.rds")))
  capture.output(print(nma, digits=4), file=file.path(odir, paste0(outcome, "_model_print.txt")))
  capture.output(summary(nma), file=file.path(odir, paste0(outcome, "_model_summary.txt")))

  core <- core_stats(nma, outcome, dat)
  write.csv(core, file.path(odir, paste0(outcome, "_core_stats.csv")), row.names=FALSE)
  vic <- extract_vs_ic(nma, "Primary")
  write.csv(vic, file.path(odir, paste0(outcome, "_vs_IC.csv")), row.names=FALSE)

  # Ranking
  rank_obj <- tryCatch(netrank(nma, small.values="undesirable"), error=function(e) NULL)
  if (!is.null(rank_obj) && !is.null(rank_obj$ranking.random)) {
    ranks <- data.frame(treatment=names(rank_obj$ranking.random), Pscore=as.numeric(rank_obj$ranking.random))
    ranks <- ranks[order(-ranks$Pscore), , drop=FALSE]
    write.csv(ranks, file.path(odir, paste0(outcome, "_Pscore.csv")), row.names=FALSE)
    tryCatch({
      png(file.path(odir, paste0(outcome, "_Pscore.png")), width=2200, height=1500, res=300, bg="white")
      par(mar=c(8,5,2,1))
      barplot(setNames(ranks$Pscore, ranks$treatment), las=2, ylim=c(0,1), ylab="P-score", main=paste(outcome, "P-score ranking"))
      dev.off()
    }, error=function(e) try(dev.off(), silent=TRUE))
  }

  # League data and network geometry
  league <- make_league(nma)
  write.csv(cbind(treatment=rownames(league), league), file.path(odir, paste0(outcome, "_league.csv")), row.names=FALSE)
  nt <- network_tables(dat)
  write.csv(nt$nodes, file.path(odir, paste0(outcome, "_network_nodes.csv")), row.names=FALSE)
  write.csv(nt$edges, file.path(odir, paste0(outcome, "_network_edges.csv")), row.names=FALSE)

  # Global inconsistency (protected so it can never stop primary NMA)
  dd <- tryCatch(decomp.design(nma), error=function(e) e)
  capture.output(print(dd), file=file.path(odir, paste0(outcome, "_design_decomposition.txt")))
  # Local inconsistency (protected)
  ns <- tryCatch(netsplit(nma), error=function(e) e)
  capture.output(print(ns), file=file.path(odir, paste0(outcome, "_netsplit.txt")))

  # Sensitivity: exclude high RoB
  dat_hr <- keep_ic_component(dat[!dat$studlab %in% high_rob_studies, , drop=FALSE])
  if (nrow(dat_hr) > 0 && length(unique(c(dat_hr$treat1, dat_hr$treat2))) >= 2) {
    fit_hr <- tryCatch(fit_nma(dat_hr), error=function(e) NULL)
    if (!is.null(fit_hr)) {
      write.csv(extract_vs_ic(fit_hr, "Exclude high RoB"), file.path(odir, paste0(outcome, "_sensitivity_highRoB.csv")), row.names=FALSE)
      saveRDS(fit_hr, file.path(odir, paste0(outcome, "_sensitivity_highRoB.rds")))
    }
  }

  # Sensitivity: exclude crossover
  dat_x <- keep_ic_component(dat[dat$studlab != "Truque-Díaz 2025, crossover", , drop=FALSE])
  if (nrow(dat_x) > 0 && length(unique(c(dat_x$treat1, dat_x$treat2))) >= 2) {
    fit_x <- tryCatch(fit_nma(dat_x), error=function(e) NULL)
    if (!is.null(fit_x)) {
      write.csv(extract_vs_ic(fit_x, "Exclude crossover"), file.path(odir, paste0(outcome, "_sensitivity_crossover.csv")), row.names=FALSE)
      saveRDS(fit_x, file.path(odir, paste0(outcome, "_sensitivity_crossover.rds")))
    }
  }

  # Leave-one-study-out
  loo <- list()
  for (st in sort(unique(dat$studlab))) {
    d <- keep_ic_component(dat[dat$studlab != st, , drop=FALSE])
    if (nrow(d) == 0 || length(unique(c(d$treat1,d$treat2))) < 2) next
    fit <- tryCatch(fit_nma(d), error=function(e) NULL)
    if (is.null(fit)) next
    z <- extract_vs_ic(fit, paste0("Exclude: ", st))
    z$excluded_study <- st
    loo[[length(loo)+1]] <- z
  }
  if (length(loo)) write.csv(do.call(rbind, loo), file.path(odir, paste0(outcome, "_leave_one_out.csv")), row.names=FALSE)

  save_optional_outputs(nma, outcome, odir)
  list(nma=nma, core=core, vs_ic=vic)
}

strength_arms <- read_embedded_csv("study,treatment,n,mean,sd,source\nCalatayud 2020 / Cruz-Montecinos 2023,EX,10,12.72,7.143287437850708,\"Calatayud 2020, Table 3\"\nCalatayud 2020 / Cruz-Montecinos 2023,IC,10,-5.08,7.597606110512448,\"Calatayud 2020, Table 3\"\nCuesta-Barriuso 2014a,MAN,11,-0.05,0.153,\"Cuesta-Barriuso 2014, Table 3\"\nCuesta-Barriuso 2014a,EX,10,-0.053000000000000005,0.157,\"Cuesta-Barriuso 2014, Table 3\"\nCuesta-Barriuso 2014a,IC,10,0,0,\"Cuesta-Barriuso 2014, Table 3\"\nDonoso-Úbeda 2023,MAN,26,280.83,77.2,\"Donoso-Úbeda 2023, Table II\"\nDonoso-Úbeda 2023,IC,26,271.33,26.26,\"Donoso-Úbeda 2023, Table II\"\nDonoso-Úbeda 2025,MAN,29,26.29,34.39,\"Donoso-Úbeda 2024/2025, Table 2\"\nDonoso-Úbeda 2025,IC,29,-2.29,4.78,\"Donoso-Úbeda 2024/2025, Table 2\"\nGoto 2014,DIG,16,1.44,0.7,\"Goto 2014, Table 2\"\nGoto 2014,EX,16,1.37,0.7,\"Goto 2014, Table 2\"\nParhampour 2022,EX,28,7.299999999999999,0.958772363204345,\"Parhampour 2022, Table 3\"\nParhampour 2022,IC,14,5.69,0.5,\"Parhampour 2022, Table 3\"\nRunkel 2016 / Runkel 2017,EX,24,5.05,1.9,\"Runkel 2016, Figure 2b\"\nRunkel 2016 / Runkel 2017,IC,28,4.8,2,\"Runkel 2016, Figure 2b\"\nTat 2021,MULTI,8,17.1,2.5925925925925926,\"Tat 2021, Table 3\"\nTat 2021,EX,9,16.8,7.037037037037036,\"Tat 2021, Table 3\"\nVillalón-González 2025/2026,EX,16,484.63,197.93,\"Villalón-González 2025/2026, Table 3\"\nVillalón-González 2025/2026,IC,17,430.42,152.19,\"Villalón-González 2025/2026, Table 3\"\n")

strength_direct <- read_embedded_csv("study,type,treat1,treat2,a,b,c,d,n1,n2,md,lower,upper,direction,source\n")

pain_arms <- read_embedded_csv("study,treatment,n,mean,sd,source\nCuesta-Barriuso 2014a,MAN,11,-1.1,1.818,\"Cuesta-Barriuso 2014, Table 3\"\nCuesta-Barriuso 2014a,EX,10,-0.44699999999999995,0.831,\"Cuesta-Barriuso 2014, Table 3\"\nCuesta-Barriuso 2014a,IC,10,-0.441,0.682,\"Cuesta-Barriuso 2014, Table 3\"\nCuesta-Barriuso 2017,EX,10,-0.32,0.79,\"Cuesta-Barriuso 2017, Table 3\"\nCuesta-Barriuso 2017,IC,10,-0.3,0.56,\"Cuesta-Barriuso 2017, Table 3\"\nCuesta-Barriuso 2020 / Donoso-Úbeda 2020,MAN,33,-1.74,1.94,\"Donoso-Úbeda 2020, Table 3\"\nCuesta-Barriuso 2020 / Donoso-Úbeda 2020,IC,32,-4.08,2.92,\"Donoso-Úbeda 2020, Table 3\"\n\"Cuesta-Barriuso 2021, knee\",MAN,14,-2,1.2600000000000002,\"Cuesta-Barriuso 2021 knee, Table 3\"\n\"Cuesta-Barriuso 2021, knee\",IC,14,-4.86,2.74,\"Cuesta-Barriuso 2021 knee, Table 3\"\nCuesta-Barriuso elbow family,MAN,35,-0.6,1.3,\"Cuesta-Barriuso 2022, Table 2\"\nCuesta-Barriuso elbow family,IC,34,-2.8,2.2,\"Cuesta-Barriuso 2022, Table 2\"\nDonoso-Úbeda 2023,MAN,26,-1.8599999999999999,2.05,\"Donoso-Úbeda 2023, Table II\"\nDonoso-Úbeda 2023,IC,26,-2.09,1.32,\"Donoso-Úbeda 2023, Table II\"\nDonoso-Úbeda 2025,MAN,29,0.68,0.88,\"Donoso-Úbeda 2024/2025, Table 2\"\nDonoso-Úbeda 2025,IC,29,0.4,1.8,\"Donoso-Úbeda 2024/2025, Table 2\"\nGoto 2014,DIG,16,-16.7,15.3,\"Goto 2014, Table 2\"\nGoto 2014,EX,16,-5.1,8.8,\"Goto 2014, Table 2\"\nKhami 2020,EPM,20,-1.33,1.44,\"Khami 2020, Table 3\"\nKhami 2020,IC,16,-5.19,2.71,\"Khami 2020, Table 3\"\nPires Oliveira 2020,ACU,15,-2.23,2.54,\"Pires Oliveira 2020, Table 2\"\nPires Oliveira 2020,EPM,13,-5.1899999999999995,2.48,\"Pires Oliveira 2020, Table 2\"\nTat 2021,MULTI,8,0,1.4814814814814814,\"Tat 2021, Table 3\"\nTat 2021,EX,9,-4,2.222222222222222,\"Tat 2021, Table 3\"\nTat 2026,MULTI,18,-2.31,2.12,\"Tat 2026, Table 5\"\nTat 2026,IC,17,-3.8799999999999994,2.23,\"Tat 2026, Tables 4–5\"\n\"Truque-Díaz 2025, pilot\",MAN,11,0,3.5555555555555554,\"Truque-Díaz 2025 pilot, Table 3\"\n\"Truque-Díaz 2025, pilot\",IC,13,-1.5,2.9999999999999996,\"Truque-Díaz 2025 pilot, Table 3\"\nUcero-Lozano 2024,IMV,24,-2.34,1.474578793421815,\"Ucero-Lozano 2024, Table III\"\nUcero-Lozano 2024,IC,24,-3.43,1.6054585679858813,\"Ucero-Lozano 2024, Table III\"\nUcero-Lozano 2025,IMV,14,-1.08,0.8757061434319141,\"Ucero-Lozano 2025, Table 2\"\nUcero-Lozano 2025,IC,14,-1.85,0.8260206884854225,\"Ucero-Lozano 2025, Table 2\"\n")

pain_direct <- read_embedded_csv("study,type,treat1,treat2,a,b,c,d,n1,n2,md,lower,upper,direction,source\nCalatayud 2020 / Cruz-Montecinos 2023,binary_or_to_smd,EX,IC,7,3,2,8,10,10,,,,higher_better,Calatayud 2020 clinically important pain improvement\n\"Truque-Díaz 2025, crossover\",md_ci_to_smd,MAN,IC,,,,,6,7,-0.85,-2.01,0.29,lower_better,Truque-Díaz 2025 crossover phase 1\n")

hjhs_arms <- read_embedded_csv("study,treatment,n,mean,sd,source\nCalatayud 2020 / Cruz-Montecinos 2023,EX,10,2.89,4.13779468024229,\"Calatayud 2020, Table 3\"\nCalatayud 2020 / Cruz-Montecinos 2023,IC,10,-0.99,4.13779468024229,\"Calatayud 2020, Table 3\"\nCuesta-Barriuso 2020 / Donoso-Úbeda 2020,MAN,33,-5.45,3.07,\"Donoso-Úbeda 2020, Table 3\"\nCuesta-Barriuso 2020 / Donoso-Úbeda 2020,IC,32,-7.03,3.27,\"Donoso-Úbeda 2020, Table 3\"\n\"Cuesta-Barriuso 2021, knee\",MAN,14,-8,2.18,\"Cuesta-Barriuso 2021 knee, Table 3\"\n\"Cuesta-Barriuso 2021, knee\",IC,14,-9.429999999999998,2.13,\"Cuesta-Barriuso 2021 knee, Table 3\"\nCuesta-Barriuso elbow family,MAN,35,-8.4,3.6,\"Cuesta-Barriuso 2022, Table 2\"\nCuesta-Barriuso elbow family,IC,34,-9.9,4.4,\"Cuesta-Barriuso 2022, Table 2\"\nFeng 2026,DIG,8,-31.63,13.14,Feng 2026\nFeng 2026,EX,9,-22.11,8.73,Feng 2026\nKhami 2020,EPM,20,-9.65,3.15,\"Khami 2020, Table 3\"\nKhami 2020,IC,16,-14.19,4.52,\"Khami 2020, Table 3\"\nMakkar 2025,MULTI,16,-22.5,7.5,\"Makkar 2025, Table 3\"\nMakkar 2025,EX,15,-21.76,7.9,\"Makkar 2025, Table 3\"\nParhampour 2019,EX,36,-4.760000000000001,1.1191872050733962,\"Parhampour 2019, Table 3\"\nParhampour 2019,IC,12,-9.08,1.67,\"Parhampour 2019, Table 3\"\nPires Oliveira 2020,ACU,15,-64.07,17.3,\"Pires Oliveira 2020, Table 2\"\nPires Oliveira 2020,EPM,13,-60.92,17.28,\"Pires Oliveira 2020, Table 2\"\nTat 2021,MULTI,8,7,3.7037037037037033,\"Tat 2021, Table 5\"\nTat 2021,EX,9,0,2.222222222222222,\"Tat 2021, Table 5\"\nTat 2026,MULTI,18,-8.17,9.48,\"Tat 2026, Table 5\"\nTat 2026,IC,17,-14.41,11.3,\"Tat 2026, Table 4\"\n\"Truque-Díaz 2025, pilot\",MAN,11,-7,5.185185185185185,\"Truque-Díaz 2025 pilot, Table 3\"\n\"Truque-Díaz 2025, pilot\",IC,13,-7,3.333333333333333,\"Truque-Díaz 2025 pilot, Table 3\"\n")

hjhs_direct <- read_embedded_csv("study,type,treat1,treat2,a,b,c,d,n1,n2,md,lower,upper,direction,source\n\"Truque-Díaz 2025, crossover\",md_ci_to_smd,MAN,IC,,,,,6,7,-2.39,-2.84,-1.94,lower_better,Truque-Díaz 2025 crossover phase 1\n")


high_rob <- c(
  "Goto 2014",
  "Khami 2020",
  "Pires Oliveira 2020",
  "Runkel 2016 / Runkel 2017",
  "Tat 2021"
)

res_strength <- analyse_outcome("Strength", strength_arms, strength_direct, high_rob)
res_pain     <- analyse_outcome("Pain", pain_arms, pain_direct, high_rob)
res_hjhs     <- analyse_outcome("HJHS", hjhs_arms, hjhs_direct, high_rob)

summary_core <- rbind(res_strength$core, res_pain$core, res_hjhs$core)
write.csv(summary_core, file.path(outdir, "summary_core_stats.csv"), row.names=FALSE)

all_vs_ic <- rbind(
  transform(res_strength$vs_ic, outcome="Strength"),
  transform(res_pain$vs_ic, outcome="Pain"),
  transform(res_hjhs$vs_ic, outcome="HJHS")
)
write.csv(all_vs_ic, file.path(outdir, "all_primary_vs_IC.csv"), row.names=FALSE)

writeLines("R NMA completed successfully.", file.path(outdir, "SUCCESS.txt"))
cat("R NMA completed successfully.\n")
print(summary_core)
print(all_vs_ic)
