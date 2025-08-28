###################################################
## HIV latency modeling - Data Loading (Sample Version)
## Note: Actual experimental data are not included.
##       This script demonstrates data handling flow.
###################################################

# 必要なライブラリ
library(FME)
library(tidyverse)
library(gridExtra)
library(grid)
library(coda)

###################################################
## Execution options
###################################################
read_data     <- TRUE
exec_fit      <- TRUE
plot_fit      <- TRUE
exec_MCMC     <- TRUE
plot_MCMC     <- TRUE
exec_coda     <- TRUE
exec_globsens <- TRUE
plot_globsens <- TRUE
code_version  <- ""
# プロットの基本テーマ
pltsetd <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "bold"),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5)
  )

###################################################
## Basic configuration
###################################################
mcmc_iter   <- 5000
mcmc_burnin <- mcmc_iter / 2
Tmin        <- 0
Tmax        <- 6
step_size   <- 0.05
stime       <- seq(Tmin, Tmax, step_size)

lbt <- 0.01
ubt <- 100

# 初期パラメータ（サンプル）
pars <- c(
  w = 0.5,
  condition1_u2 = 0.1,
  condition1_u3 = 0.1,
  condition1_u4 = 0.1,
  condition1_u5 = 0.1,
  condition1_u6 = 0.1,
  condition1_g  = 1,
  condition2_u2 = 0.1,
  condition2_u3 = 0.1,
  condition2_u4 = 0.1,
  condition2_u5 = 0.1,
  condition2_u6 = 0.1,
  condition2_g  = 1,
  condition3_u2 = 0.1,
  condition3_u3 = 0.1,
  condition3_u4 = 0.1,
  condition3_u5 = 0.1,
  condition3_u6 = 0.1,
  condition3_g  = 1,
  abp = 3.3, arn = 0.1, fb = 1, fp = 1, fr = 1
)

estimated_pars <- c(k = 1.0e7)  # サンプル値

###################################################
## Data reading (サンプルデータを読み込む)
###################################################
# 実データは非公開。本コードでは data/ 内のサンプルを読み込む。

get_infected_data <- function(folder_name) {
  result <- vector("list", length = 3)
  for (i in 1:3) {
    file_path <- paste0("data/", folder_name, "/summary_R", i, ".txt")
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    raw <- read.csv(file_path, sep = "\t", header = TRUE)
    
    # 整形（色ごとの強度と細胞数）
    df <- raw %>%
      group_by(Time, Color) %>%
      summarise(
        Intensity   = sum(Intensity, na.rm = TRUE),
        Cell_Count  = sum(Cell_Count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from  = Color,
        values_from = c(Intensity, Cell_Count),
        values_fill = list(Intensity = 0, Cell_Count = 0)
      ) %>%
      mutate(
        time = Time,
        Fb = Intensity_Blue,
        Fp = Intensity_Purple,
        Fr = Intensity_Red,
        xib = Cell_Count_Blue,
        xip = Cell_Count_Purple,
        xlr = Cell_Count_Red,
        xln = Cell_Count_NoColor
      ) %>%
      select(time, xib, xip, xlr, xln, Fb, Fp, Fr)
    
    result[[i]] <- df
  }
  return(result)
}

# 各条件のデータ読み込み
condition1 <- get_infected_data("condition1")
condition2 <- get_infected_data("condition2")
condition3 <- get_infected_data("condition3")

###################################################
## Data merge & averaging
###################################################
all_data <- lapply(1:3, function(i) {
  df_list <- list(
    condition1 = condition1[[i]],
    condition2 = condition2[[i]],
    condition3 = condition3[[i]]
  )
  # 列名に条件名を付与して結合
  df_list_renamed <- lapply(names(df_list), function(nm) {
    df <- df_list[[nm]]
    colnames(df)[-1] <- paste0(colnames(df)[-1], "_", nm)
    df
  })
  Reduce(function(x, y) merge(x, y, by = "time"), df_list_renamed)
})

# 平均値を計算
calculate_mean_dataframe <- function(dfs) {
  combined <- do.call(cbind, dfs)
  mean_df <- data.frame(lapply(1:ncol(dfs[[1]]), function(i) {
    rowMeans(combined[, seq(i, ncol(combined), ncol(dfs[[1]]))], na.rm = TRUE)
  }))
  colnames(mean_df) <- colnames(dfs[[1]])
  mean_df
}

condition1_mean <- calculate_mean_dataframe(condition1)
condition2_mean <- calculate_mean_dataframe(condition2)
condition3_mean <- calculate_mean_dataframe(condition3)


###################################################
## ODE model
###################################################
system(paste("R CMD SHLIB sub_analysis", code_version,  ".cpp", sep = ""))
dyn.load(paste("sub_analysis", code_version, .Platform$dynlib.ext, sep = ""))

ODEs_devided <- function(pars) {
  Tmin <- 1
  segments <- c(2, 3, 4, 5, 6)
  
  # 初期値を condition1,2,3 の平均データから取得
  output <- data.frame(
    time = Tmin,
    xib_condition1 = condition1_mean$xib[condition1_mean$time == Tmin],
    xip_condition1 = condition1_mean$xip[condition1_mean$time == Tmin],
    xlr_condition1 = condition1_mean$xlr[condition1_mean$time == Tmin],
    xln_condition1 = condition1_mean$xln[condition1_mean$time == Tmin],
    xib_condition2 = condition2_mean$xib[condition2_mean$time == Tmin],
    xip_condition2 = condition2_mean$xip[condition2_mean$time == Tmin],
    xlr_condition2 = condition2_mean$xlr[condition2_mean$time == Tmin],
    xln_condition2 = condition2_mean$xln[condition2_mean$time == Tmin],
    xib_condition3 = condition3_mean$xib[condition3_mean$time == Tmin],
    xip_condition3 = condition3_mean$xip[condition3_mean$time == Tmin],
    xlr_condition3 = condition3_mean$xlr[condition3_mean$time == Tmin],
    xln_condition3 = condition3_mean$xln[condition3_mean$time == Tmin],
    Fb_condition1  = condition1_mean$Fb[condition1_mean$time == Tmin],
    Fp_condition1  = condition1_mean$Fp[condition1_mean$time == Tmin],
    Fr_condition1  = condition1_mean$Fr[condition1_mean$time == Tmin],
    Fb_condition2  = condition2_mean$Fb[condition2_mean$time == Tmin],
    Fp_condition2  = condition2_mean$Fp[condition2_mean$time == Tmin],
    Fr_condition2  = condition2_mean$Fr[condition2_mean$time == Tmin],
    Fb_condition3  = condition3_mean$Fb[condition3_mean$time == Tmin],
    Fp_condition3  = condition3_mean$Fp[condition3_mean$time == Tmin],
    Fr_condition3  = condition3_mean$Fr[condition3_mean$time == Tmin],
    segment = paste0("Segment ", Tmin)
  )
  
  # ODE用初期条件
  rhs <- as.numeric(output[1, -c(1, ncol(output))])
  names(rhs) <- names(output)[-c(1, ncol(output))]
  
  # パラメータ設定
  parms <- c(
    condition1_w = pars["w"],
    condition1_u2 = pars["condition1_u2"],
    condition1_u3 = pars["condition1_u3"],
    condition1_u4 = pars["condition1_u4"],
    condition1_u5 = pars["condition1_u5"],
    condition1_u6 = pars["condition1_u6"],
    condition1_g  = pars["condition1_g"],
    
    condition2_w = pars["w"],
    condition2_u2 = pars["condition2_u2"],
    condition2_u3 = pars["condition2_u3"],
    condition2_u4 = pars["condition2_u4"],
    condition2_u5 = pars["condition2_u5"],
    condition2_u6 = pars["condition2_u6"],
    condition2_g  = pars["condition2_g"],
    
    condition3_w = pars["w"],
    condition3_u2 = pars["condition3_u2"],
    condition3_u3 = pars["condition3_u3"],
    condition3_u4 = pars["condition3_u4"],
    condition3_u5 = pars["condition3_u5"],
    condition3_u6 = pars["condition3_u6"],
    condition3_g  = pars["condition3_g"],
    
    k   = estimated_pars["k"],
    abp = pars["abp"], arn = pars["arn"],
    fb  = pars["fb"],  fp  = pars["fp"], fr  = pars["fr"],
    segment = Tmin
  )
  
  # 各セグメントでODEを解く
  for (seg in segments) {
    times <- seq(Tmin, seg, by = 0.05)
    odes_times <- seq(0, seg - Tmin, by = 0.05)
    
    parms["segment"] <- seg
    
    # 外部C++コードでODEを計算
    out <- ode(
      y = rhs, parms = parms, times = odes_times,
      func = "derivs", initfunc = "initparms",
      nout = 9,
      outnames = c(
        "Fb_condition1", "Fp_condition1", "Fr_condition1",
        "Fb_condition2", "Fp_condition2", "Fr_condition2",
        "Fb_condition3", "Fp_condition3", "Fr_condition3"
      ),
      dllname = paste0("sub_analysis", code_version)
    )
    out <- as.data.frame(out)
    out$time <- times
    out$segment <- paste0("Segment ", seg)
    
    output <- rbind(output, out)
    rhs <- tail(out, 1)[, names(rhs)] * 0.5  # passageで細胞数を半減
    Tmin <- seg
  }
  
  return(output)
}

###################################################
## Sum of Squared Residuals (SSR)
## - Compare model output with experimental data
###################################################
SSR <- function(pars) {
  val <- 0.0
  out <- ODEs_devided(pars)
  compare_times <- c(2, 3, 4, 5, 6)
  cell_take_out <- 0.5
  
  selected_columns <- grep("condition", colnames(out), value = TRUE)
  
  for (i in 1:3) {
    datavalues <- log10(all_data[[i]][, -1] + 1e-20)
    datavalues$time <- all_data[[i]][, "time"]
    datavalues[datavalues < -1] <- 0  # ノイズのクリップ
    
    for (j in compare_times) {
      row_model <- out[out$time == j & out$segment == paste0("Segment ", j), selected_columns]
      row_data  <- datavalues[datavalues$time == j, selected_columns]
      diff <- (log10(row_model * cell_take_out + 1e-20) - row_data)^2
      val <- val + sum(diff, na.rm = TRUE)
    }
  }
  return(val)
}

###################################################
## Log-SSR wrapper
###################################################
SSRlog <- function(lpars) {
  SSR(10^lpars)
}

###################################################
## Constraint setting (for parameter fitting)
###################################################
set_constr <- function(low_vec, up_vec) {
  n <- length(low_vec)
  stopifnot(length(up_vec) == n)
  AL <- diag(n); AU <- -diag(n)
  Amat <- rbind(AL, AU)
  bvec <- c(low_vec, -up_vec)
  return(list(low.val = low_vec, up.val = up_vec, Amat = Amat, bvec = bvec))
}

###################################################
## Sensitivity wrapper
###################################################
sensrfun <- function(pars) {
  ODEs_devided(pars)
}

###################################################
## Parameter fitting (Nelder-Mead optimization)
###################################################
if (exec_fit == TRUE) {
  lbtvec <- lbt * pars
  ubtvec <- ubt * pars
  cbox   <- set_constr(log10(lbtvec), log10(ubtvec))
  
  # Constrained optimization (log-scale)
  est <- constrOptim(
    theta   = log10(pars),
    f       = SSRlog,
    method  = "Nelder-Mead",
    ui      = cbox$Amat,
    ci      = cbox$bvec,
    control = list(maxit = 500, reltol = 1e-16)
  )
  fit_mixture <- 10^(est$par)
  save(fit_mixture, file = paste0("fit_mixture", code_version, ".Rdata"))
  write.table(fit_mixture, "fit_mixture.txt", sep = "\t")
}

###################################################
## MCMC sampling
###################################################
load(paste0("fit_mixture", code_version, ".Rdata"))
chain <- 3

if (exec_MCMC == TRUE) {
  mcmc_pars <- fit_mixture
  results <- mclapply(1:chain, function(c) {
    MCMC_mixture <- modMCMC(
      f     = SSRlog,
      p     = log10(mcmc_pars),
      lower = log10(lbt * mcmc_pars),
      upper = log10(ubt * mcmc_pars),
      niter = mcmc_iter,
      burninlength = mcmc_burnin,
      jump  = 0.075 * log10(fit_mixture),
      wvar0 = 0.2,
      updatecov = 50
    )
    MCMC_mixture$pars <- 10^(MCMC_mixture$pars)
    save(MCMC_mixture, file = paste0("MCMC_mixture", code_version, c, ".Rdata"))
    return(MCMC_mixture)
  }, mc.cores = detectCores() - 1)
  names(results) <- paste0("MCMC_mixture", 1:chain)
}

###################################################
## Convergence diagnostics (Gelman-Rubin)
###################################################
mcmc_list    <- list()
mcmclog_list <- list()
for (c in 1:chain) {
  load(paste0("MCMC_mixture", code_version, c, ".Rdata"))
  cobj    <- as.mcmc(MCMC_mixture$pars)
  cobjlog <- log10(cobj)
  mcmc_list[[c]]    <- cobj
  mcmclog_list[[c]] <- cobjlog
}
mcmc_list    <- mcmc.list(mcmc_list)
mcmclog_list <- mcmc.list(mcmclog_list)

# Gelman-Rubin PSRF
gelman_result    <- as.data.frame(gelman.diag(mcmc_list)$psrf)
gelmanlog_result <- as.data.frame(gelman.diag(mcmclog_list)$psrf)
write.table(gelman_result,    "gelman_rubin_results.txt", sep = "\t")
write.table(gelmanlog_result, "gelman_rubin_results_log10.txt", sep = "\t")

###################################################
## Credible intervals (95% CI)
###################################################
load(paste0("MCMC_mixture", code_version, "1.Rdata"))
csv <- MCMC_mixture$pars[order(MCMC_mixture$SS),][1:MCMC_mixture$naccepted,]
quantiles <- apply(csv, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
write.csv(quantiles, "95credible_intervals.csv")

###################################################
## Global sensitivity analysis
###################################################
sR_pars <- fit_mixture
load(paste0("MCMC_mixture", code_version, "1.Rdata"))
sR <- sensRange(func = sensrfun, parms = sR_pars, parInput = MCMC_mixture$pars)

# 列名整理
convert_col_name <- function(name) {
  num_dots <- nchar(gsub("[^.]", "", name))
  if (num_dots > 1) sub("\\.(?!.*\\.)", "00", name, perl = TRUE) else name
}
colnames(sR) <- sapply(colnames(sR), convert_col_name)
sR <- sR[, !grepl("segment", names(sR))]
sR[] <- lapply(sR, function(x) as.numeric(as.character(x)))


###################################################
## Visualization
## - Compare model fit with data
## - Plot estimated reactivation parameters
###################################################

# Basic plot theme
pltsetd <- theme_classic() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 20, hjust = 0.5)
  )

xlabel <- "Days"
ylabel <- "Cells (log10)"

# Run model with fitted parameters
fitted <- ODEs_devided(fit_mixture)
fitted_log10 <- fitted %>%
  mutate(across(-c(time, segment), ~ log10(. + 1e-20)))

# Experimental data (log10変換済み)
data_list <- lapply(all_data, function(df) {
  data.frame(time = df$time, log10(df[, -1] + 1e-20))
})

# === Plot model vs data for selected variables ===
plot_vars <- c("xib_condition1", "xip_condition1", "xlr_condition1",
               "xib_condition2", "xip_condition2", "xlr_condition2",
               "xib_condition3", "xip_condition3", "xlr_condition3")

color_map <- c("xib"="#72C7F3","xip"="#A753A6","xlr"="#F566AF")

plots <- list()
for (var in plot_vars) {
  prefix <- strsplit(var, "_")[[1]][1]
  col <- color_map[[prefix]]
  
  # 実験データ平均＋SE
  ted_combined <- bind_rows(data_list) %>%
    group_by(time) %>%
    summarise(
      mean_value = mean(.data[[var]], na.rm = TRUE),
      se = sd(.data[[var]], na.rm = TRUE)/sqrt(n()),
      ymin = mean_value - se,
      ymax = mean_value + se
    )
  
  p <- ggplot() +
    geom_line(data = fitted_log10, aes(x=time, y=.data[[var]]),
              color=col, size=1.2, alpha=0.7) +
    geom_point(data = ted_combined, aes(x=time, y=mean_value),
               color=col, size=3) +
    geom_errorbar(data = ted_combined, aes(x=time, ymin=ymin, ymax=ymax),
                  width=0.2, color=col) +
    labs(title=var, x=xlabel, y=ylabel) +
    scale_y_continuous(limits=c(0,6)) +
    pltsetd
  plots[[var]] <- p
}

# Save combined figure
combined_plot <- grid.arrange(grobs = plots, ncol=3)
ggsave("combined_plot.png", combined_plot, width=12, height=8)

###################################################
## Reactivation parameter plot (median + 95% CI)
###################################################
load(paste0("MCMC_mixture", code_version, "1.Rdata"))
accepted_samples <- MCMC_mixture$pars[order(MCMC_mixture$SS),][1:MCMC_mixture$naccepted,]

medians <- apply(accepted_samples, 2, median)
ci_lower <- apply(accepted_samples, 2, function(x) quantile(x,0.025))
ci_upper <- apply(accepted_samples, 2, function(x) quantile(x,0.975))

# 条件ごとのラベル
u_params <- grep("u[2-6]", names(medians), value=TRUE)
groups <- ifelse(grepl("condition1", u_params),"Condition1",
                 ifelse(grepl("condition2", u_params),"Condition2","Condition3"))

df <- data.frame(
  param = u_params,
  group = groups,
  median = medians[u_params],
  ci_lower = ci_lower[u_params],
  ci_upper = ci_upper[u_params]
)

reactivation_plot <- ggplot(df, aes(x=param, y=median, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), color="black") +
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper),
                width=0.2, position=position_dodge(0.9)) +
  scale_fill_manual(values=c("Condition1"="#72C7F3",
                             "Condition2"="#A753A6",
                             "Condition3"="#F566AF")) +
  labs(x="Parameter", y="Reactivation rate (/day)") +
  pltsetd +
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggsave("reactivation_pars_with_CI.png", reactivation_plot, width=8, height=6)
