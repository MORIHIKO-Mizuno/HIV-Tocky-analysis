rm(list = ls(all = TRUE))

# 文字列結合のための関数, 意味わからん
"+" <- function(e1, e2) {
  if (is.character(c(e1, e2))) {
    paste(e1, e2, sep = "")
  } else {
    base::"+"(e1, e2)
  }
}

library(FME)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
library(coda)
library(dplyr)
###################################################
## Execution option
###################################################

read_data <- "yes"
exec_fit <- "yes"
plot_fit <- "yes"

exec_MCMC <- "yes"
plot_MCMC <- "yes"
exec_coda <- "yes"
exec_globsens <- "yes"
plot_globsens <- "yes"

datatypes <- c(1, 2, 3)
pltsetd <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "bold"),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 17, color = "black", face = "bold"),
    axis.text.y = element_text(size = 17, color = "black", face = "bold"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
  )
###################################################
## Basic configuration
###################################################

# コードのバージョンを指定

mcmc_iter <-  500000 # 500000
mcmc_burnin <- mcmc_iter / 2
Tmin <- 0
Tmax <- 6
step_size <- 0.05
stime <- seq(Tmin, Tmax, step_size)

lbt <- 1 / 100 # lower bound threshold
ubt <- 100 # upper bound threshold

# ラムダが大きいなら遅延の項を消してもいい
pars <- c(
  w = 0.5,

  JQ1_u2_infected = 0.1,
  JQ1_u3_infected = 0.1,
  JQ1_u4_infected = 0.1,
  JQ1_u5_infected = 0.1,
  JQ1_u6_infected = 0.1,
  JQ1_g_infected = 1,
  # notreat_w_infected = 4.00417647,

  notreat_u2_infected = 0.1,
  notreat_u3_infected = 0.1,
  notreat_u4_infected = 0.1,
  notreat_u5_infected = 0.1,
  notreat_u6_infected = 0.1,
  notreat_g_infected = 1,
  # stim_w_infected = 2.24420237,

  stim_u2_infected = 0.1,
  stim_u3_infected = 0.1,
  stim_u4_infected = 0.1,
  stim_u5_infected = 0.1,
  stim_u6_infected = 0.1,
  stim_g_infected = 1,
  abp = 3.3,
  arn = 0.1,
  fb = 1,
  fp = 1,
  fr = 1
)
estimated_pars <- c(
  JQ1_g_infected = 1.296807, notreat_g_infected = 1.206803, stim_g_infected = 8.018361e-01,
  k = 3.207954e+08
)
code_version <- ""

setwd("/Users/2m/Library/ode" + code_version)

###################################################
## Data reading
###################################################
blue_line <- 10**2.91812
red_line <- 10**2.93463

get_infected_data <- function(folder_name) {
  int_data <- vector(mode = "list", length = 3)
  for (i in c(1:3)) {
    ### koko!!
    rfn <- paste("data/", folder_name, "/summary_R", as.character(i), ".txt", sep = "")
    print(rfn)
    dsf <- read.csv(rfn, sep = "\t", header = T)
    
    
    # まず、Total_Blue と Total_Red の値を Color ごとに展開する
    data_wide <- dsf %>%
      group_by(Time, Color) %>%
      summarise(
        Intensity = sum(Intensity, na.rm = TRUE),
        Cell_Count = sum(Cell_Count, na.rm = TRUE), .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Color,
        values_from = c(Intensity, Cell_Count),
        values_fill = list(Intensity = 0, Cell_Count = 0)
      )
    
    # 新しい列 (Fbb, Fpb, Fpr, Frr) の計算
    final_data <- data_wide %>%
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
    dtime <- final_data[[1]] # または final_data$time
    tn <- length(dtime) # nrow() の代わりに length() を使用
    
    int_data[[i]] <- final_data
  }
  
  # 各データフレームからtimeが0の行を抽出
  
  return(int_data)
}

JQ1_infected <- get_infected_data("JQ1_infected")
notreat_infected <- get_infected_data("notreat_infected")
stim_infected <- get_infected_data("stim_infected")

# SSRのために全てのデータを結合
all_data <- vector(mode = "list", length = 3)
for (i in c(1, 2, 3)) {
  # 例として、データフレームが格納されたリストを定義
  df_list <- list(
    JQ1_infected = JQ1_infected[[i]],
    notreat_infected = notreat_infected[[i]],
    stim_infected = stim_infected[[i]]
  )
  
  # 列名を変更する関数
  rename_columns <- function(df, df_name) {
    colnames(df)[-1] <- paste(colnames(df)[-1], df_name, sep = "_")
    return(df)
  }
  
  # リスト内の各データフレームの列名を変更
  df_list_renamed <- lapply(seq_along(df_list), function(i) {
    rename_columns(df_list[[i]], names(df_list)[i])
  })
  
  # リスト内のデータフレームをtime列で結合
  merged_df <- Reduce(function(x, y) merge(x, y, by = "time"), df_list_renamed)
  
  all_data[[i]] <- merged_df
}
# データフレームの平均を計算する関数
calculate_mean_dataframe <- function(dataframe_list) {
  # データフレームの数を取得
  n <- length(dataframe_list)
  
  # 各データフレームの行数と列数が同じであることを確認
  if (!all(sapply(dataframe_list, nrow) == nrow(dataframe_list[[1]])) ||
      !all(sapply(dataframe_list, ncol) == ncol(dataframe_list[[1]]))) {
    stop("All dataframes must have the same dimensions.")
  }
  
  # データフレームの結合
  combined_df <- do.call(cbind, dataframe_list)
  
  # 平均を計算する
  mean_df <- data.frame(lapply(1:ncol(dataframe_list[[1]]), function(i) {
    rowMeans(combined_df[, seq(i, ncol(combined_df), ncol(dataframe_list[[1]]))])
  }))
  
  # 列名を設定
  colnames(mean_df) <- colnames(dataframe_list[[1]])
  
  return(mean_df)
}
# 初期値を設定するためにデータフレームを作成
JQ1_mean <- calculate_mean_dataframe(JQ1_infected)
notreat_mean <- calculate_mean_dataframe(notreat_infected)
stim_mean <- calculate_mean_dataframe(stim_infected)
###################################################
## functions
###################################################

# cでの処理
## fatal ##
system(paste("R CMD SHLIB LogM", code_version, ".cpp", sep = ""))
dyn.load(paste("LogM", code_version, .Platform$dynlib.ext, sep = ""))

## ODE model ##

ODEs_devided <- function(pars) {
  # Calculate differential equations for each time period
  Tmin <- 1
  segments <- c(2, 3, 4, 5, 6)
  output <- data.frame(
    time = Tmin,
    xib_JQ1_infected = JQ1_mean$xib[JQ1_mean$time == Tmin],
    xip_JQ1_infected = JQ1_mean$xip[JQ1_mean$time == Tmin],
    xlr_JQ1_infected = JQ1_mean$xlr[JQ1_mean$time == Tmin],
    xln_JQ1_infected = JQ1_mean$xln[JQ1_mean$time == Tmin],
    xib_notreat_infected = notreat_mean$xib[notreat_mean$time == Tmin],
    xip_notreat_infected = notreat_mean$xip[notreat_mean$time == Tmin],
    xlr_notreat_infected = notreat_mean$xlr[notreat_mean$time == Tmin],
    xln_notreat_infected = notreat_mean$xln[notreat_mean$time == Tmin],
    xib_stim_infected = stim_mean$xib[stim_mean$time == Tmin],
    xip_stim_infected = stim_mean$xip[stim_mean$time == Tmin],
    xlr_stim_infected = stim_mean$xlr[stim_mean$time == Tmin],
    xln_stim_infected = stim_mean$xln[stim_mean$time == Tmin],
    Fb_JQ1_infected = JQ1_mean$Fb[JQ1_mean$time == Tmin],
    Fp_JQ1_infected = JQ1_mean$Fp[JQ1_mean$time == Tmin],
    Fr_JQ1_infected = JQ1_mean$Fr[JQ1_mean$time == Tmin],
    Fb_notreat_infected = notreat_mean$Fb[notreat_mean$time == Tmin],
    Fp_notreat_infected = notreat_mean$Fp[notreat_mean$time == Tmin],
    Fr_notreat_infected = notreat_mean$Fr[notreat_mean$time == Tmin],
    Fb_stim_infected = stim_mean$Fb[stim_mean$time == Tmin],
    Fp_stim_infected = stim_mean$Fp[stim_mean$time == Tmin],
    Fr_stim_infected = stim_mean$Fr[stim_mean$time == Tmin],
    segment = paste0("Segment ", Tmin)
  )
  # initial conditions
  rhs <- c(
    xib_JQ1_infected = JQ1_mean$xib[JQ1_mean$time == Tmin],
    xip_JQ1_infected = JQ1_mean$xip[JQ1_mean$time == Tmin],
    xlr_JQ1_infected = JQ1_mean$xlr[JQ1_mean$time == Tmin],
    xln_JQ1_infected = JQ1_mean$xln[JQ1_mean$time == Tmin],
    xib_notreat_infected = notreat_mean$xib[notreat_mean$time == Tmin],
    xip_notreat_infected = notreat_mean$xip[notreat_mean$time == Tmin],
    xlr_notreat_infected = notreat_mean$xlr[notreat_mean$time == Tmin],
    xln_notreat_infected = notreat_mean$xln[notreat_mean$time == Tmin],
    xib_stim_infected = stim_mean$xib[stim_mean$time == Tmin],
    xip_stim_infected = stim_mean$xip[stim_mean$time == Tmin],
    xlr_stim_infected = stim_mean$xlr[stim_mean$time == Tmin],
    xln_stim_infected = stim_mean$xln[stim_mean$time == Tmin]
  )
  parms <- c(
    JQ1_w_infected = as.numeric(pars["w"]),
  
    JQ1_u2_infected = as.numeric(pars["JQ1_u2_infected"]),
    JQ1_u3_infected = as.numeric(pars["JQ1_u3_infected"]),
    JQ1_u4_infected = as.numeric(pars["JQ1_u4_infected"]),
    JQ1_u5_infected = as.numeric(pars["JQ1_u5_infected"]),
    JQ1_u6_infected = as.numeric(pars["JQ1_u6_infected"]),
    JQ1_g_infected = as.numeric(pars["JQ1_g_infected"]),
    
    notreat_w_infected = as.numeric(pars["w"]),

    notreat_u2_infected = as.numeric(pars["notreat_u2_infected"]),
    notreat_u3_infected = as.numeric(pars["notreat_u3_infected"]),
    notreat_u4_infected = as.numeric(pars["notreat_u4_infected"]),
    notreat_u5_infected = as.numeric(pars["notreat_u5_infected"]),
    notreat_u6_infected = as.numeric(pars["notreat_u6_infected"]),
    notreat_g_infected = as.numeric(pars["notreat_g_infected"]),
    
    stim_w_infected = as.numeric(pars["w"]),

    stim_u2_infected = as.numeric(pars["stim_u2_infected"]),
    stim_u3_infected = as.numeric(pars["stim_u3_infected"]),
    stim_u4_infected = as.numeric(pars["stim_u4_infected"]),
    stim_u5_infected = as.numeric(pars["stim_u5_infected"]),
    stim_u6_infected = as.numeric(pars["stim_u6_infected"]),
    stim_g_infected = as.numeric(pars["stim_g_infected"]),
    k = as.numeric(estimated_pars["k"]),
    abp = as.numeric(pars["abp"]),
    arn = as.numeric(pars["arn"]),
    fb = as.numeric(pars["fb"]),
    fp = as.numeric(pars["fp"]),
    fr = as.numeric(pars["fr"]),
    segment = Tmin
  )

  for (i in 1:length(segments)) {
    Tmax <- segments[i]
    times <- seq(Tmin, Tmax, by = 0.05)
    odes_times <- seq(0, Tmax - Tmin, by = 0.05)
    parms["segment"] <- segments[i]
    ## C compiled code_version ##
    out <- ode(
      y = rhs, parms = parms, times = odes_times, func = "derivs", initfunc = "initparms",
      nout = 9, outnames = c(
        "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
        "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
        "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
      ),
      dllname = paste("LogM", code_version, sep = "")
    )
    out <- as.data.frame(out)
    colnames(out) <- c(
      "time",
      "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
      "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
      "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected",
      "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
      "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
      "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
    )
    
    out$time <- c(times)
    out$segment <- paste0("Segment ", segments[i])
    output <- rbind(output, out)
    rhs <- tail(out, 1)
    # resetting of initial values
    ratio <- 1 / 2
    rhs <- rhs[, c(
      "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
      "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
      "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected"
    )] * ratio
    rhs <- setNames(as.numeric(rhs), names(rhs))
    Tmin <- Tmax
  }
  # 最後の行を取得
  last <- tail(out, 1)
  
  # ratio の設定
  ratio <- 1 / 2
  
  # 特定の列に対して ratio を掛け算
  cols_to_modify <- which(!names(last) %in% c("time", "segment"))
  last[, cols_to_modify] <- last[, cols_to_modify] * ratio  
  # 'segment' 列の値を 'Segment 7' に設定
  last$segment <- "Segment 7"
  
  # 修正された last データフレームを表示
  print(last)
  output <- rbind(output, last)
  return(output)
}

SSR <- function(pars) {
  val <- 0.0
  out <- ODEs_devided(pars)
  compare_times <- c(2, 3, 4, 5, 6)
  cell_take_out <- 1 / 2
  selected_columns <- c(
    "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
    "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
    "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected",
    "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
    "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
    "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
  )
  for (i in 1:3) {
    times <- all_data[[i]][, 1]
    datavalues <- log10(all_data[[i]][, 1:ncol(all_data[[i]])] + 10^-20)
    datavalues$time <- times
    # データが変なやつを邪魔にならないように修正
    datavalues[datavalues < -1] <- 0
    for (j in t(compare_times)) {
      # データのフィルタリング
      row <- out[out$time == j & out$segment == paste("Segment ", as.character(j), sep = ""), ]
      row <- row[, -which(names(row) == "segment")]
      
      datavalues_row <- datavalues[datavalues$time == j, ]
      row_selected <- row[, selected_columns]
      datavalues_row_selected <- datavalues_row[, selected_columns]
      
      difference <- (log10(row_selected * cell_take_out + 10^-20) - datavalues_row_selected)^2
      val <- val + sum(difference, na.rm = TRUE)
    }
  }
  return(val)
}

## log-SSR ##
SSRlog <- function(lpars) {
  SSR(10^(lpars))
}

set_constr <- function(low_vec, up_vec) {
  n_vec <- length(low_vec)
  stopifnot(length(up_vec) == n_vec)
  ## for constrOptim(): we can use the same rule to interexchange components so that we do not have to do it.##
  AL <- diag(n_vec)
  AU <- -diag(n_vec)
  Amat <- rbind(AL, AU) # 'ui'
  bvec <- c(low_vec, -up_vec) # 'ci'
  return(list(low.val = low_vec, up.val = up_vec, Amat = Amat, bvec = bvec))
}

## sensitivity function ##
sensrfun <- function(pars) {
  return(ODEs_devided(pars))
}

###################################################
## Parameter fitting
###################################################

# 制限条件を設定しているらしい
if (exec_fit == "yes") {
  lbtvec <- lbt * pars
  ubtvec <- ubt * pars
  loglbtvec <- log10(lbtvec)
  logubtvec <- log10(ubtvec)
  cbox <- set_constr(loglbtvec, logubtvec)
  
  # 大数のスケールでパラメータを探索する
  est <- constrOptim(
    theta = log10(pars), f = SSRlog, method = "Nelder-Mead", grad = NULL,
    ui = cbox$Amat, ci = cbox$bvec,
    control = list(trace = 0, maxit = 500, reltol = 1e-16)
  )
  fit_mixture <- 10^(est$par)
  save(fit_mixture, file = paste("fit_mixture", code_version, ".Rdata", sep = ""))
}

write.table(fit_mixture, "fit_mixture.txt", sep = "\t")

###################################################
## Plot fitted curve
###################################################

# Function for parameter validation

pltsetd <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "bold"),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 17, color = "black", face = "bold"),
    axis.text.y = element_text(size = 17, color = "black", face = "bold"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
  )


#########################
##テストプロット
#########################
# 共通の設定
ypmin <- 0
ypmax <- 6.0
fuchi <- "00000000"
calpha <- "AA"
xlabel <- "days"
ylabel <- "cells(Log10)"

# データフレームの準備
# pathデータ
fitted <- ODEs_devided(fit_mixture)
# f=ODEs_devided(fit)
# f_log10 <- f %>%
#   mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
# 培地分割前後
ex_fitting_points <- fitted %>%
  filter(time %in% c(0.25, 0.5, 1, 2, 3, 4, 5, 6, 7) & segment == paste0("Segment ", time))
fitting_points <- ex_fitting_points %>%
  mutate(across(
    .cols = -c(time, segment), # timeとsegment列を除外
    .fns = ~ if_else(time %in% c(0.25, 0.5, 1, 2, 3, 4, 5, 6, 7), .x / 2, .x)
  ))
fitted_log10 <- fitted %>%
  mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
ex_fitting_points_log10 <- ex_fitting_points %>%
  mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
fitting_points_log10 <- fitting_points %>%
  mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
titex <- paste("estimated pars (SSR = ", SSR(pars))
data_list <- vector(mode = "list", length = length(datatypes))
auxiliary_df <- fitted_log10 %>%
  group_by(time) %>%
  # セグメントが複数存在するtimeのみを抽出
  filter(n_distinct(segment) > 1) %>%
  # vertical列にtimeの値を入れる
  mutate(vertical = time) %>%
  ungroup()
for (i in seq_along(all_data)) {
  ## experimental data ##
  i_data <- all_data[[i]] # リストのi番目の要素にアクセス
  data_list[[i]] <- data.frame(
    time = i_data[, "time"],
    log10(i_data[, 2:ncol(i_data)] + 10^-20)
  )
}
ted1 <- data_list[[1]]
ted2 <- data_list[[2]]
ted3 <- data_list[[3]]

# 全てのカラムに対しての処理、プロット
plot_names <- c(
  "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
  "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
  "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected",
  "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
  "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
  "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
)
color_mapping <- list(
  xib = "#0000FF", # blue
  xip = "#800080", # purple
  xlr = "#FF0000", # red
  xln = "#808080", # grey
  Fb = "#0000FF", # blue
  Fp = "#800080", # 青紫
  Fr = "#FF0000" # red
)

# plot_namesの各エントリに対応するカラーを取得
plot_colors <- sapply(plot_names, function(name) {
  prefix <- strsplit(name, "_")[[1]][1] # "xib_JQ1_uninfected" の "xib" 部分を抽出
  color_mapping[[prefix]] # 対応するカラーを取得
})


for (i in 1:length(plot_names)) {
  plot_name <- plot_names[i]
  print(plot_name)
  plot_color <- plot_colors[i]
  
  sfn <- paste0("srplot", code_version, "", plot_name, ".png")
  png(sfn, width = 800, height = 800, bg = "transparent")
  
  maintitle <- paste(stringr::str_to_title(plot_name), "cells")
  
  # データフレームを結合し、指定した色の平均と標準誤差を計算
  ted_combined <- bind_rows(data_list) %>%
    group_by(time) %>%
    summarize(
      mean_value = mean(.data[[plot_name]], na.rm = TRUE),
      se = sd(.data[[plot_name]], na.rm = TRUE) / sqrt(n()),
      ymin = mean_value - se,
      ymax = mean_value + se
    )
  
  # プロット
  plot <- ggplot() +
    scale_size_manual(values = c(1, 1)) +
    scale_color_manual(values = c(fuchi, fuchi)) +
    geom_path(data = fitted_log10, aes(x = time, y = .data[[plot_name]]), color = plot_color, lwd = 2, alpha = 0.4) +
    # geom_path(data = f_log10, aes(x = time, y = .data[[plot_name]]), color = plot_color, lwd = 2, alpha = 0.2) +
    geom_point(data = fitting_points_log10, aes(x = time, y = .data[[plot_name]]), colour = plot_color, size = 5, shape = 3) +
    geom_point(data = ted_combined, aes(x = time, y = mean_value), colour = plot_color, size = 5) +
    geom_errorbar(data = ted_combined, aes(x = time, ymin = ymin, ymax = ymax), width = 0.2, colour = plot_color) +
    xlab(xlabel) +
    ylab(ylabel) +
    labs(title = maintitle) +
    pltsetd +
    scale_x_continuous(breaks = seq(0, 7, by = 1), labels = c("day 0", "day 1", "day 2", "day 3", "day 4", "day 5", "day 6", "day 7")) +
    scale_y_continuous(limits = c(0, 8))
  
  print(plot)
  dev.off()
}
print("前半終わり！後半も頑張ろう！")
print(fit_mixture)

###################################################
## ここがMCMCだよ！！！！！！！！！！！！！！！！！
###################################################
load(paste("fit_mixture", code_version, ".Rdata", sep = ""))
current_day <- Sys.time()
print(current_day)

library(parallel)


tail <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15")
chain <- 3
if (exec_MCMC == "yes") {
  mcmc_pars <- fit_mixture
  
  # mclapplyを使用してMCMCの実行
  results <- mclapply(1:chain, function(c) {
    sttime <- Sys.time()
    cat(sprintf("chain: %d, st: %s\n", c, sttime))
    
    MCMC_mixture <- modMCMC(
      f = SSRlog, p = log10(mcmc_pars),
      lower = log10(lbt * (mcmc_pars)),
      upper = log10(ubt * (mcmc_pars)),
      niter = mcmc_iter, burninlength = mcmc_burnin,
      jump = 0.1 * log10(fit_mixture), var0 = NULL, wvar0 = 0.1, updatecov = 50
    )
    MCMC_mixture$pars <- 10^(MCMC_mixture$pars)
    
    save(MCMC_mixture, file = paste0("MCMC_mixture", code_version, c, ".Rdata"))
    edtime <- Sys.time()
    cat(sprintf("chain: %d, ed: %s\n", c, edtime))
    
    return(MCMC_mixture)
  }, mc.cores = detectCores() - 1)
  
  # 結果の保存や処理
  names(results) <- paste("MCMC_mixture", 1:chain, sep = "")
}


# 平均値と分散値を計算
m_mcmc <- do.call(rbind, lapply(results, function(x) apply(x$pars, 2, mean)))
v_mcmc <- do.call(rbind, lapply(results, function(x) apply(x$pars, 2, var)))

w_mcmc <- 1 / mcmc_burnin / (chain - 1) * apply(v_mcmc, 2, sum)
b_mcmc <- mcmc_burnin / (chain - 1) * apply(m_mcmc, 2, var)
vhat_mcmc <- ((mcmc_burnin - 1) * w_mcmc + b_mcmc) / mcmc_burnin
rhat_mcmc <- sqrt(vhat_mcmc / w_mcmc)


mcmc_list <- list()
mcmclog_list <- list()
for (c in 1:chain) {
  if (plot_MCMC == "yes") {
    load(paste("MCMC_mixture", code_version, c, ".Rdata", sep = ""))
    png(paste("MCMC_mixture", code_version, c, ".png", sep = ""), width = 1200, height = 1200)
    plot(MCMC_mixture, Full = TRUE)
    pairs(MCMC_mixture, nsample = 1000)
    dev.off()
    png(paste("MCMC_mixture_log10", code_version, c, ".png", sep = ""), width = 1200, height = 1200)
    MCMC_mixture_log10 <- MCMC_mixture
    MCMC_mixture_log10$pars <- log10(MCMC_mixture_log10$pars)
    plot(MCMC_mixture_log10, Full = TRUE)
    pairs(MCMC_mixture_log10, nsample = 1000)
    dev.off()
  }
  if (exec_coda == "yes") {
    library(coda)
    load(paste("MCMC_mixture", code_version, c, ".Rdata", sep = ""))
    cobj <- as.mcmc(MCMC_mixture$pars)
    pdf(paste("MCMC_mixture_conv", code_version, c, ".pdf", sep = ""))
    plot(cobj)
    dev.off()
    mcmc_list[[c]] <- cobj
    cobjlog <- log10(cobj)
    pdf(paste("MCMC_mixture_conv_log10", code_version, c, ".pdf", sep = ""))
    plot(cobjlog)
    dev.off()
    mcmclog_list[[c]] <- cobjlog
  }
}
mcmc_list <- mcmc.list(mcmc_list)
mcmclog_list <- mcmc.list(mcmclog_list)
# ゲルマン-ルービン統計量の計算結果をデータフレームに変換
gelman_result <- as.data.frame(gelman.diag(mcmc_list)$psrf)
gelmanlog_result <- as.data.frame(gelman.diag(mcmclog_list)$psrf)
# テキストファイルとして保存
write.table(gelman_result, "gelman_rubin_results.txt", sep = "\t")
write.table(gelmanlog_result, "gelman_rubin_results_log10.txt", sep = "\t")

csv <- MCMC_mixture$pars[order(MCMC_mixture$SS),][1:MCMC_mixture$naccepted,]
csv <- apply(csv,2,function(x){quantile(x,c(0.01,0.025,0.5,0.975,0.99))})
csvp <-apply(csv,2,function(x){paste(x[2],"-",x[4],seP="")})
csvppp <- rbind(csv,csvp)
write.csv(csvppp,paste("95credible_intervals",".csv",sep=""))
###################################################
## Sensitivity ranges
###################################################
load(paste("fit_mixture", code_version, ".Rdata", sep = ""))

sR_pars <- fit_mixture
load(paste("MCMC_mixture", code_version, "1.Rdata", sep = ""))

# fatal #
sR <- sensRange(func = sensrfun, parms = sR_pars, parInput = MCMC_mixture$pars)

# 列名の変換関数
convert_col_name <- function(col_name) {
  # 最後の"."を数える
  num_dots <- nchar(gsub("[^.]", "", col_name))
  
  # 最後の"."を削除する
  if (num_dots > 1) {
    return(sub("\\.(?!.*\\.)", "00", col_name, perl = TRUE))
  } else {
    return(col_name)
  }
}


# 元のデータ mat の列名を変換
new_col_names <- sapply(colnames(sR), convert_col_name)
colnames(sR) <- new_col_names
# セグメントを削除
sR <- sR[, !grepl("segment", names(sR))]
# なぜか文字列になっているので変換
sR[] <- lapply(sR, function(x) as.numeric(as.character(x)))
# save(sR, file = paste("MCMC_globsens", code_version, ".Rdata", sep = ""))
# pdf(paste("MCMC_globsens", code_version, ".pdf", sep = ""))
# plot(summary(sR), xlab = "time")
# dev.off()


pltsetd <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", face = "bold"),
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 17, color = "black", face = "bold"),
    axis.text.y = element_text(size = 17, color = "black", face = "bold"),
    axis.title = element_text(size = 18, color = "black", face = "bold"),
  )


load(paste("fit_mixture", code_version, ".Rdata", sep = ""))
results_text <- SSR(fit_mixture)
write.table(results_text, "SSR.txt", sep = "\t")

#### Paper ####
#### JSMB ####
#### JSMB ####
# 共通の設定
ypmin <- 0
ypmax <- 6.0
fuchi <- "00000000"
calpha <- "54"
xlabel <- "days"
ylabel <- "cells(Log10)"

pltsetd <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black", family = "Arial"),  # 全体のフォントファミリをArialに
    axis.ticks = element_line(colour = "black"),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 25, color = "black", hjust = 0.5, family = "Arial"),  # タイトルにArialを適用
    axis.text.x = element_text(size = 17, color = "black", family = "Arial"),  # x軸テキストにArialを適用
    axis.text.y = element_text(size = 17, color = "black", family = "Arial"),  # y軸テキストにArialを適用
    axis.title = element_text(size = 18, color = "black", family = "Arial")  # 軸タイトルにArialを適用
    
  )


# データフレームの準備
# pathデータ
# MCMC の結果から受理されたサンプルを取得
accepted_samples <- MCMC_mixture$pars[order(MCMC_mixture$SS),][1:MCMC_mixture$naccepted,]

# 各パラメータの中央値を計算
medians <- apply(accepted_samples, 2, function(x) quantile(x, 0.5))
fitted <- ODEs_devided(medians)

# 培地分割前後
# ex_fitting_points <- fitted %>%
#   filter(time %in% c(0.25, 0.5, 1, 2, 3, 4, 5, 6, 7) & segment == paste0("Segment ", time))
# fitting_points <- ex_fitting_points %>%
#   mutate(across(
#     .cols = -c(time, segment), # timeとsegment列を除外
#     .fns = ~ if_else(time %in% c(0.25, 0.5, 1, 2, 3, 4, 5, 6, 7), .x / 2, .x)
#   ))
fitted_log10 <- fitted %>%
  mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
# ex_fitting_points_log10 <- ex_fitting_points %>%
#   mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
# fitting_points_log10 <- fitting_points %>%
#   mutate(across(.cols = -c(time, segment), .fns = ~ log10(. + 10^-20)))
titex <- paste("estimated pars (SSR = ", SSR(pars))
data_list <- vector(mode = "list", length = length(datatypes))
auxiliary_df <- fitted_log10 %>%
  group_by(time) %>%
  # セグメントが複数存在するtimeのみを抽出
  filter(n_distinct(segment) > 1) %>%
  # vertical列にtimeの値を入れる
  mutate(vertical = time) %>%
  ungroup()
for (i in seq_along(all_data)) {
  ## experimental data ##
  i_data <- all_data[[i]] # リストのi番目の要素にアクセス
  data_list[[i]] <- data.frame(
    time = i_data[, "time"],
    log10(i_data[, 2:ncol(i_data)] + 10^-20)
  )
}
ted1 <- data_list[[1]]
ted2 <- data_list[[2]]
ted3 <- data_list[[3]]

# 全てのカラムに対しての処理、プロット
plot_names <- c(
  "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
  "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
  "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected",
  "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
  "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
  "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
)
# color_mapping <- list(
#   xib = "#0000FF", # blue
#   xip = "#800080", # purple
#   xlr = "#FF0000", # red
#   xln = "#808080", # grey
#   Fb = "#0000FF", # blue
#   Fp = "#800080", # 青紫
#   Fr = "#FF0000" # red
# )
color_mapping <- list(
  xib = "#72C7F3", # blue
  xip = "#A753A6", # purple
  xlr = "#F566AF", # red
  xln = "#A3A3A3", # grey
  Fb = "#72C7F3", # blue
  Fp = "#A753A6", # 青紫
  Fr = "#F566AF" # red
)

# plot_namesの各エントリに対応するカラーを取得
plot_colors <- sapply(plot_names, function(name) {
  prefix <- strsplit(name, "_")[[1]][1] # "xib_JQ1_uninfected" の "xib" 部分を抽出
  color_mapping[[prefix]] # 対応するカラーを取得
})
plots=c()

for (i in 1:length(plot_names)) {
  plot_name <- plot_names[i]
  print(plot_name)
  plot_color <- plot_colors[i]
  name_lst <- colnames(sR)
  gids <- grep(paste0("^", plot_name), name_lst)
  mat <- sR[, gids]
  
  # "xib_JQ1_uninfected"を含む要素のインデックスを見つける
  blue_indices <- grep("xib_JQ1_infected", name_lst)
  
  # "xib_JQ1_uninfected"文字列を含む要素から数字だけを抽出する
  blue_values <- as.numeric(gsub("xib_JQ1_infected", "", name_lst[blue_indices]))
  
  # blue_valuesを出力
  times <- blue_values
  
  # matが数値でない要素をNAに変換
  mat <- apply(mat, 2, function(x) as.numeric(x))
  
  # 再度quantileを計算
  ymean <- log10(apply(mat, 2, mean, na.rm = TRUE))
  yCIlow <- log10(apply(mat, 2, function(x) quantile(x, 0.025, na.rm = TRUE)))
  yCIhigh <- log10(apply(mat, 2, function(x) quantile(x, 0.975, na.rm = TRUE)))
  ylow <- pmax(as.numeric(yCIlow), 0)
  yup <- as.numeric(yCIhigh)
  
  dc <- data.frame(x = times, ymin = ylow, ymax = yup)
  
  # sfn <- paste0("srplot", code_version, "", plot_name, ".png")
  # png(sfn, width = 800, height = 800, bg = "transparent")
  
  maintitle <- paste(stringr::str_to_title(plot_name), "cells")
  
  # データフレームを結合し、指定した色の平均と標準誤差を計算
  ted_combined <- bind_rows(data_list) %>%
    group_by(time) %>%
    summarize(
      mean_value = mean(.data[[plot_name]], na.rm = TRUE),
      se = sd(.data[[plot_name]], na.rm = TRUE) / sqrt(n()),
      ymin = mean_value - se,
      ymax = mean_value + se
    )
  
  
  # プロット
  plot <- ggplot() +
    geom_ribbon(data = dc, aes(x = x, ymin = ymin, ymax = ymax), fill = paste(plot_color, calpha, sep = "")) +
    geom_line(
      data = data.frame(x = rep(times, 2), y = c(ylow, yup), curves = gl(2, length(times), label = c("5%", "95%"))),
      aes(x = x, y = y, color = curves), lwd = 1
    ) +
    scale_size_manual(values = c(1, 1)) +
    scale_color_manual(values = c(fuchi, fuchi)) +
    geom_path(data = fitted_log10, aes(x = time, y = .data[[plot_name]]), color = plot_color, lwd = 2, alpha = 0.6) +
    # geom_point(data = fitting_points_log10, aes(x = time, y = .data[[plot_name]]), colour = plot_color, size = 5, shape = 3) +
    # WARNING:
    geom_point(data=ted1, aes(x = time, y = .data[[plot_name]]), colour = plot_color, size = 5) +
    geom_point(data=ted2, aes(x = time, y = .data[[plot_name]]), colour = plot_color, size = 5) +
    geom_point(data=ted3, aes(x = time, y = .data[[plot_name]]), colour = plot_color, size = 5) +
    # geom_point(data = fitting_points_log10, aes(x = time, y = .data[[plot_name]]), colour = plot_color, size = 5, shape = 3) +
    # geom_point(data = ted_combined, aes(x = time, y = mean_value), colour = plot_color, size = 5) +
    # geom_errorbar(data = ted_combined, aes(x = time, ymin = ymin, ymax = ymax), width = 0.2, colour = plot_color) +
    xlab("") +
    ylab("") +
    pltsetd +
    scale_x_continuous(limits = c(1, 6), breaks = seq(1, 6, by = 1), labels = c(1,2,3,4,5,6)) +
    scale_y_continuous(limits = c(0, 6)) +
    theme(axis.text.x = element_text(size = 45, family = "Arial"),
          axis.text.y = element_text(size = 45, family = "Arial"),
          axis.ticks = element_line(size = 2, color = "black"),          # ティックスの線の太さと色を指定
          axis.ticks.length = unit(0.4, "cm"),
          axis.line=element_line(size=2,color="black"))
  
  # プロットをリストに追加
  plots[[plot_name]] <- plot
  
  # プロットをファイルに保存
  png_filename <- paste0(plot_name, "data.png")
  png(filename = png_filename, width = 500, height = 500, res = 96)
  print(plot)
  dev.off()
}

# プロットの名前の順序を指定
plot_order <- c(
  "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
  "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
  "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected"

)

# plots リストからプロットを抽出する
# プロットの名前を `plot_order` に基づいて抽出
sorted_plots <- lapply(plot_order, function(name) {
  # plots リストから指定された名前のプロットを取得
  plots[[name]]
})


# グリッド形式でプロットを表示
combined_plot <- grid.arrange(grobs = sorted_plots, ncol = 4, nrow = 3)


# すべてのプロットを一つのファイルに保存
ggsave("combined_plot_3point.png", plot = combined_plot, width = 28, height = 21)

# 蛍光量のプロットの名前の順序を指定
plot_order <- c(
  "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
  "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
  "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
)

# plots リストからプロットを抽出する
# プロットの名前を `plot_order` に基づいて抽出
sorted_plots <- lapply(plot_order, function(name) {
  # plots リストから指定された名前のプロットを取得
  plots[[name]]
})

# グリッド形式でプロットを表示
combined_plot <- grid.arrange(grobs = sorted_plots, ncol = 3, nrow = 3)

# すべてのプロットを一つのファイルに保存
ggsave("combined_plot_f_3point.png", plot = combined_plot, width = 21, height = 21)

# 論文用にプロットをまとめる
plot_order <- c(
  "xib_notreat_infected", "xip_notreat_infected", "xlr_notreat_infected", "xln_notreat_infected",
  "Fb_notreat_infected", "Fp_notreat_infected", "Fr_notreat_infected",
  
  "xib_JQ1_infected", "xip_JQ1_infected", "xlr_JQ1_infected", "xln_JQ1_infected",
  "Fb_JQ1_infected", "Fp_JQ1_infected", "Fr_JQ1_infected",
  
  "xib_stim_infected", "xip_stim_infected", "xlr_stim_infected", "xln_stim_infected",
  "Fb_stim_infected", "Fp_stim_infected", "Fr_stim_infected"
)

# plots リストからプロットを抽出する
# プロットの名前を `plot_order` に基づいて抽出
sorted_plots <- lapply(plot_order, function(name) {
  # plots リストから指定された名前のプロットを取得
  plots[[name]]
})
# 'expression()' を使用して下付き文字を設定
x_labels <- c(expression(I[b]), expression(I[p]), expression(L[r]), expression(L[n]), 
              expression(F[b]), expression(F[p]), expression(F[r]))

y_labels <- c("No Treatment", "JQ1", "Stimulation")

# テーマの文字サイズ設定
custom_theme <- theme(
  plot.title = element_text(size = 45, hjust = 0.5),   # タイトルを大きくして中央揃え
  axis.title.x = element_text(size = 45),             # x軸ラベルのフォントサイズ
  axis.title.y = element_text(size = 45),             # y軸ラベルのフォントサイズ
  axis.text = element_text(size = 45)                 # 軸の目盛りのフォントサイズ
)

# 各プロットにラベルを設定
labeled_plots <- list()
for (i in seq_along(sorted_plots)) {
  # x_labelとy_labelは特定のプロットにのみ設定
  plot <- sorted_plots[[i]] + custom_theme
  
  # 1行目にのみx軸のラベルを設定
  if (i <= 7) {
    plot <- plot + labs(title = x_labels[(i - 1) %% 7 + 1])
  }
  
  # 1列目にのみy軸のラベルを設定
  if (i %% 7 == 1) {
    plot <- plot + labs(y = y_labels[ceiling(i / 7)])
  }
  
  # 最下段にのみx軸（"Days"）のラベルを設定
  if (i > 14) {
    plot <- plot + labs(x = "Days")
  }
  
  labeled_plots[[i]] <- plot
}

# グリッド形式でプロットを表示
combined_plot <- grid.arrange(grobs = labeled_plots, ncol = 7, nrow = 3)

# すべてのプロットを一つのファイルに保存
ggsave("combined_plot_paper.png", plot = combined_plot, width = 49, height = 21)

### Reactivation Pars Plot ###
### 今日も頑張ろう、歯医者に遅れないよう ###
load(paste("MCMC_mixture", code_version, "1.Rdata", sep = ""))

# MCMC の結果から受理されたサンプルを取得
accepted_samples <- MCMC_mixture$pars[order(MCMC_mixture$SS),][1:MCMC_mixture$naccepted,]

# 各パラメータの中央値と95%信用区間を計算
medians <- apply(accepted_samples, 2, function(x) quantile(x, 0.5))
ci_lower <- apply(accepted_samples, 2, function(x) quantile(x, 0.025))
ci_upper <- apply(accepted_samples, 2, function(x) quantile(x, 0.975))

# # 推定パラメータデータ
# parameters <- c(
#   JQ1_u2_infected = 0.0048242996, JQ1_u3_infected = 0.0691892784, JQ1_u4_infected = 0.1546772026, 
#   JQ1_u5_infected = 0.0250831756, JQ1_u6_infected = 0.1109752314, 
#   notreat_u2_infected = 0.0130890721, notreat_u3_infected = 0.0719461630, notreat_u4_infected = 0.1121819278, 
#   notreat_u5_infected = 0.0738113596, notreat_u6_infected = 0.1709851848, 
#   stim_u2_infected = 0.1020910958, stim_u3_infected = 1.4536138568, stim_u4_infected = 0.3495391886, 
#   stim_u5_infected = 0.0003574313, stim_u6_infected = 0.0064032638
# )

# # "u" が含まれているパラメータを抽出
# u_params <- parameters[grep("u[2-6]_infected", names(parameters))]

# tの値をn~n+1の形式に設定
t_values <- rep(c("Day1-2", "Day2-3", "Day3-4", "Day4-5", "Day5-6"), 3)

# ラベルを簡略化し、No Treatment, JQ1, Stimulationの順番で並べ替え
u_labels <- c(
  "notreat_u2_infected" , "notreat_u3_infected" , "notreat_u4_infected", 
  "notreat_u5_infected" , "notreat_u6_infected" , 
  "JQ1_u2_infected" , "JQ1_u3_infected" , "JQ1_u4_infected", 
  "JQ1_u5_infected" , "JQ1_u6_infected" ,
  "stim_u2_infected" , "stim_u3_infected" , "stim_u4_infected", 
  "stim_u5_infected" , "stim_u6_infected"
)

# カテゴリ分けのためにNo Treatment, JQ1, Stimulationで色を設定
group <- rep(c("No Treatment", "JQ1", "Stimulation"), each = 5)


# データフレームにまとめる
df <- data.frame(
  t = c("Day1-2","Day2-3","Day3-4","Day4-5","Day5-6"),
  label = u_labels,
  group = factor(group, levels = c("No Treatment", "JQ1", "Stimulation")),
  median = medians[(u_labels)],  # 対応する中央値
  ci_lower = ci_lower[(u_labels)],  # 95% CI 下限
  ci_upper = ci_upper[(u_labels)]   # 95% CI 上限
)
reactivation_plot <- ggplot(df, aes(x = t, y = median, fill = group)) + pltsetd +
  geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.7) +  # 棒グラフを描画
  # geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(0.9)) +  # エラーバー
  scale_fill_manual(values = c("No Treatment" = "white", "JQ1" = "#00C800", "Stimulation" = "#DDCA00"),
                    guide = guide_legend(title = NULL)) +
  labs(x = "", y = expression("Reactivation rate (" * day^{-1} * ")"))+
  scale_y_continuous(limits = c(0, NA)) +  # 通常スケールのy軸
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(size = 0.5, color = "black"),
    axis.text.x = element_text(size = 10.5, family = "Arial"),
    axis.text.y = element_text(size = 10.5, family = "Arial"),
    axis.title.x = element_text(size = 10.5, family = "Arial"),
    axis.title.y = element_text(size = 10.5, family = "Arial"),
    plot.title = element_text(size = 10.5, hjust = 0.5, family = "Arial"),
    axis.ticks.y = element_line(size = 0.5, color = "black"),
    axis.ticks.length.y = unit(0.2, "cm"),
    legend.position = c(0.85, 0.7),
    legend.text = element_text(size = 10.5, family = "Arial"),
    aspect.ratio = 0.25
  )

reactivation_plot

ggsave("reactivation_pars_2025.png", plot = reactivation_plot, width = 10, height = 2.5)
