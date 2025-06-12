library(tidyverse)
library(tidyr)
library(ggplot2)
pltsetd <- theme_classic() +
    theme( text = element_text(size = 14, color = "black", face = "bold"),
           axis.ticks = element_line(colour = "black"),
           axis.line = element_line(colour = "black"),
           legend.title = element_blank(),
           legend.text = element_blank(),
           legend.position = "none" , 
           plot.title = element_text(size = 25, color = "black", face = "bold", hjust = 0.5),
           axis.text.x = element_text(size = 17, color = "black", face = "bold"),
           axis.text.y = element_text(size = 17, color = "black", face = "bold"),
           axis.title = element_text(size = 18, color = "black", face = "bold"),
    )

rm(list = ls(all = TRUE)) # メモリ解放
med <- "JQ1"
parent_folder <- file.path("data", paste0(med, "_infected"))
setwd(parent_folder)
folders <- c("R1", "R2", "R3")

# 活性化領域を分けるための設定
# 青と無色の赤蛍光の最大、赤と無色の青色蛍光の最大の値を使用

blue_line <- 10**2.91812
red_line <- 10**2.93463

# ファイル名と計測時点
file_names <- c(
    "0",
    "0.25",
    "0.5",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6"
)
times <- c(0,0.25,0.5, 1, 2, 3, 4, 5, 6)
df_names <- data.frame(file_names, times)

# データのマージと計算


for (folder in folders){
    total_data <- data.frame()
    for (df_name in df_names$file_names) {
        file_path <- paste0(parent_folder, "/", folder,"/",df_name, ".csv") # ファイルのパスを設定
        data <- read.csv(file_path)
        data <- data %>%
            mutate(
            Blue = 10^Blue.log,
            Red = 10^Red.log,
            Intensity_nonLog = case_when(
                !is.na(Angle) ~ 10^Intensity,
                TRUE ~ 0
            )
        )
        
        # 色ごとにデータを分割
        data$Color <- case_when(
            data$Angle==0 ~ "Blue",
            data$Angle<90 & data$Angle>0 ~ "Purple",
            data$Angle==90 ~ "Red",
            TRUE ~ "NoColor"
        )
        
        # 計測時点とファイル名を追加
        data$Time <- df_names$times[df_names$file_names == df_name]
        data$FileName <- df_name
        
        # データをマージ
        total_data <- bind_rows(total_data, data)
    }
    # 各色ごとに細胞の蛍光量の合計と細胞の個数を計算
    summary_data <- total_data %>%
        group_by(Time, Color) %>%
        summarise(
            Intensity = sum(Intensity_nonLog, na.rm = TRUE),
            Cell_Count = n() # 細胞の個数
        )
    write.table(summary_data, file = paste("summary_",folder,".txt",sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
}
# 各色ごとに細胞の蛍光量の合計を計算

