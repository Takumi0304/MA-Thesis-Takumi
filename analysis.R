

library(tidyverse)
library(lme4)
library(lmerTest)
library(MASS)

#---- Load and pre-process data: ----

read_csv("results_final_RT.csv", col_names=FALSE) |>
  dplyr::select(
    time       = X1,
    subj       = X2,
    lable      = X6,
    chunk      = X10,
    text       = X11,
    list       = X14,
    item       = X15,
    cond       = X16,
    correct    = X17,
    mc_subj_nr = X18,
    type       = X19,
    rt         = X20,
    newline    = X21,
    sent       = X22) |>
  mutate(
    correct       = ifelse(correct=="F", FALSE, TRUE),
    newline       = ifelse(newline==0, FALSE, TRUE),
    mc_subject_nr = factor(mc_subj_nr, levels=c("SG", "PL")),
    cond          = as.integer(cond),
    mc_subj_nr    = as.factor(mc_subj_nr),
    item          = as.integer(item),
    chunk         = as.integer(chunk),
    list          = as.integer(list)) |>
  filter(type=="item") |>
  dplyr::select(-type, -newline, -lable) |>
  mutate(
    region = case_when(
      chunk==1     ~ "mc_subj",
      chunk%in%2:5 ~ "rc",
      chunk==6     ~ "mc_verb",
      chunk%in%7:8 ~ "spillover",
      TRUE         ~ NA),
    rc_type = case_when(
      cond %in% 1:2 ~ "SRC",
      cond %in% 3:4 ~ "ORC",
      TRUE          ~ NA),
    rc_type = factor(rc_type, levels=c("SRC", "ORC")),
    nr_match = case_when(
      cond %in% c(1,3) ~ "match",
      cond %in% c(2,4) ~ "mismatch",
      TRUE             ~ NA),
    nr_match = factor(nr_match, levels=c("match", "mismatch"))) -> d
head(d)

table(d$correct)



# 被験者数を数える（重複なしの被験者ID数をカウント）
num_participants_RT <- d %>%
  distinct(subj) %>%  # subj列のユニークな値を抽出
  nrow()

print(paste("RTデータの被験者数:", num_participants_RT))


#---- Sanity checks: ----

##---- Latin square lists balanced? ----
with(d, table(d$list))

#---- TODO Add more sanity checks to make sure that our assumptions about ----
## the data actually hold.




#---- Analysis: Reading Time ----
##---- Descriptive stats: ----

###---- Region 2: RC: ----
d |>
  filter(region=="rc") |>
  group_by(subj, item) |>
  summarize(
    rt         = sum(rt),
    rc_type    = first(rc_type),
    nr_match   = first(nr_match),
    mc_subj_nr = first(mc_subj_nr),
    text       = paste(text, collapse=" ")) -> d.rc

head(d.rc)

d.rc |>
  group_by(nr_match, rc_type, subj) |>
  summarize(
    rt = exp(mean(log(rt)))) |>
  summarize(
    se  = sd(log(rt))/sqrt(n()),
    rt  = exp(mean(log(rt))),
    cil = exp(log(rt)-2*se),
    ciu = exp(log(rt)+2*se),
    .groups="keep") |>
  ungroup() -> summary.rc

summary.rc

ggplot(summary.rc, aes(x=interaction(nr_match, rc_type), y=rt, ymin=cil, ymax=ciu)) +
  geom_point(size=2) +
  scale_y_continuous(limits=c(0, 6000)) +
  geom_errorbar() +
  theme_bw()


###---- Region 3: MC verb: ----
d |>
  filter(region=="mc_verb") |>
  group_by(subj, item) |>
  summarize(
    rt         = sum(rt),
    rc_type    = first(rc_type),
    nr_match   = first(nr_match),
    mc_subj_nr = first(mc_subj_nr),
    text       = paste(text, collapse=" ")) -> d.mc_verb

head(d.mc_verb)

d.mc_verb |>
  group_by(nr_match, rc_type, subj) |>
  summarize(
    rt = exp(mean(log(rt)))) |>
  summarize(
    se  = sd(log(rt))/sqrt(n()),
    rt  = exp(mean(log(rt))),
    cil = exp(log(rt)-2*se),
    ciu = exp(log(rt)+2*se),
    .groups="keep") |>
  ungroup() -> summary.mc_verb

summary.mc_verb

ggplot(summary.mc_verb, aes(x=interaction(nr_match, rc_type), y=rt, ymin=cil, ymax=ciu)) +
  geom_point(size=2) +
  scale_y_continuous(limits=c(0, 1500)) +
  geom_errorbar() +
  theme_bw()



###---- Region 4: Spillover region: ----
d |>
  filter(region=="spillover") |>
  group_by(subj, item) |>
  summarize(
    rt         = sum(rt),
    rc_type    = first(rc_type),
    nr_match   = first(nr_match),
    mc_subj_nr = first(mc_subj_nr),
    text       = paste(text, collapse=" ")) -> d.spillover

head(d.spillover)

d.spillover |>
  group_by(nr_match, rc_type, subj) |>
  summarize(
    rt = exp(mean(log(rt)))) |>
  summarize(
    se  = sd(log(rt))/sqrt(n()),
    rt  = exp(mean(log(rt))),
    cil = exp(log(rt)-2*se),
    ciu = exp(log(rt)+2*se),
    .groups="keep") |>
  ungroup() -> summary.spillover

summary.spillover

ggplot(summary.spillover, aes(x=interaction(nr_match, rc_type), y=rt, ymin=cil, ymax=ciu)) +
  geom_point(size=2) +
  scale_y_continuous(limits=c(0, 6000)) +
  geom_errorbar() +
  theme_bw()



##---- Inferential stats: ----

####---- Region 2: RC / LMM ----
contrasts(d.rc$nr_match) <- contr.sdif(2)
contrasts(d.rc$rc_type)  <- contr.sdif(2)

lmer(log(rt) ~ nr_match * rc_type + (1|subj) + (1|item), d.rc) -> m1
summary(m1)

####---- Region 2: RC / Table----
library(dplyr)
library(gt)

# 手入力データ
table_lmm_rc <- tribble(
  ~Term,                            ~Estimate, ~Std_Error, ~df,  ~t_value,  ~p_value,
  "Intercept (SRC, match)",            8.338,      0.100,  20.350,   83.500,  "<.001",
  "Number (mismatch – match)",        -0.044,      0.028,  714.157,  -1.574,  "0.116",
  "RC Type (ORC – SRC)",               0.110,      0.028,  699.156,   3.974,  "<.001",
  "Interaction: RC × Number",         -0.044,      0.055,  699.156,  -0.801,  "0.423"
)

table_gt_rc <- table_lmm_rc %>%
  gt() %>%
  tab_header(
    title = md("**Results of Linear Mixed Model (Region 2: Relative Clause)**")
  ) %>%
  cols_label(
    Term = md("**Predictor / Fixed effect**"),
    Estimate = md("**Estimate**"),
    Std_Error = md("**Standard error**"),
    df = md("**df**"),
    t_value = md("**t value**"),
    p_value = md("**p value**")
  ) %>%
  fmt_number(
    columns = c(Estimate, Std_Error, t_value),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = c(p_value),
      rows = p_value == "<.001"
    )
  ) %>%
  cols_align(
    align = "center",
    columns = c(Estimate, Std_Error, df, t_value, p_value)
  ) %>%
  opt_table_outline() %>%
  opt_align_table_header("center") %>%
  tab_options(
    table.font.names = "Arial",
    table.font.size = 12,
    data_row.padding = px(4)
  )

gtsave(table_gt_rc, filename = "table_lmm_rc.png")



####---- Region 2: Graph ----
# まずは d のカラムを確認
colnames(d)

# 二つあるEnglish_Proficiencyのどちらを使うかの見極め
table(d$English_Proficiency.x)
table(d$English_Proficiency.y)

# どちらも同じだったので統合
# 統合列を作成し、不要な列は削除
d <- d %>%
  mutate(English_Proficiency = coalesce(English_Proficiency.x, English_Proficiency.y)) %>%
  dplyr::select(-English_Proficiency.x, -English_Proficiency.y)

d.rc <- d %>%
  filter(region == "rc") %>%
  group_by(subj, item) %>%
  summarize(
    rt = sum(rt),
    rc_type = first(rc_type),
    nr_match = first(nr_match),
    mc_subj_nr = first(mc_subj_nr),
    text = paste(text, collapse = " "),
    English_Proficiency = first(English_Proficiency),
    .groups = "drop"
  )

# Make a graph
library(ggplot2)

# nr_match と rc_type を結合して条件ラベルを作る
d.rc <- d.rc %>%
  mutate(Condition = interaction(nr_match, rc_type, sep = "-"))

# English_Proficiencyを因子にして表示順を整える（必要なら）
d.rc$English_Proficiency <- factor(d.rc$English_Proficiency, levels = c("B1", "B2", "C1", "C2"))

# 条件と英語力ごとの平均RTと95%CIを計算
summary.rc <- d.rc %>%
  group_by(Condition, English_Proficiency) %>%
  summarize(
    mean_rt = exp(mean(log(rt))),
    se_log = sd(log(rt)) / sqrt(n()),
    cil = exp(log(mean_rt) - 2 * se_log),
    ciu = exp(log(mean_rt) + 2 * se_log),
    .groups = "drop"
  )

# バープロット作成
ggplot(summary.rc, aes(x = Condition, y = mean_rt, fill = English_Proficiency)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = cil, ymax = ciu),
                position = position_dodge(width = 0.8),
                width = 0.2) +
  labs(x = "Condition", y = "Reading Time (ms)", fill = "English Proficiency") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###---- Region 3: MC verb: ----
contrasts(d.mc_verb$nr_match) <- contr.sdif(2)
contrasts(d.mc_verb$rc_type)  <- contr.sdif(2)

lmer(log(rt) ~ nr_match * rc_type + (1|subj) + (1|item), d.mc_verb) -> m2
summary(m2)

####---- Table-Graph Region 3: MC Verb: ----
library(dplyr)
library(gt)

## Manually inputting the values
table_lmm_mc_verb <- tribble(
  ~Term,                            ~Estimate, ~Std_Error, ~df,  ~t_value,  ~p_value,
  "Intercept (SRC, match)",            6.768,      0.076,   27.645, 88.475,  "<.001",
  "Number (mismatch – match)",         0.038,      0.036,  703.699,  1.058,  "0.290",
  "RC Type (ORC – SRC)",               0.006,      0.035,  697.731,  0.165,  "0.869",
  "Interaction: RC × Number",         -0.013,      0.070,  697.731, -0.178,  "0.859"
)


# Apply gt and make the table neat
table_lmm_mc_verb %>%
  gt() %>%
  tab_header(
    title = "Results of LMM (Region 3: Main Clause Verb)"
  ) %>%
  fmt_number(
    columns = c(Estimate, Std_Error, t_value),
    decimals = 3
  ) %>%
  cols_label(
    Term = "Predictor / Fixed effect",
    Estimate = "Estimate",
    Std_Error = "Standard error",
    df = "df",
    t_value = "t value",
    p_value = "p value"
  )

# 1. making gt table based on table before and save it in object
table_gt_mc_verb <- table_lmm_mc_verb %>%
  gt() %>%
  tab_header(
    title = "Results of LMM (Region 3: Main Clause Verb)"
  ) %>%
  fmt_number(
    columns = c(Estimate, Std_Error, t_value),
    decimals = 3
  ) %>%
  cols_label(
    Term = "Predictor / Fixed effect",
    Estimate = "Estimate",
    Std_Error = "Standard error",
    df = "df",
    t_value = "t value",
    p_value = "p value"
  )

# 2. Save as PNG file
gtsave(table_gt_mc_verb, filename = "table_lmm_mc_verb.png")


###---- Region 4: Spillover region: ----
contrasts(d.spillover$nr_match) <- contr.sdif(2)
contrasts(d.spillover$rc_type)  <- contr.sdif(2)

lmer(log(rt) ~ nr_match * rc_type + (1|subj) + (1|item), d.spillover) -> m3
summary(m3)

####---- Table-Graph Region 4: Spillover: ----
library(dplyr)
library(gt)

# Manually inputting the values
table_lmm_spillover <- tribble(
  ~Term,                            ~Estimate, ~Std_Error, ~df,  ~t_value,  ~p_value,
  "Intercept (SRC, match)",            7.710,      0.098,  24.103,  78.287,  "<.001",
  "Number (mismatch – match)",        -0.060,      0.030,  702.455, -2.003,  "0.046",
  "RC Type (ORC – SRC)",               0.021,      0.030,  697.621,  0.700,  "0.484",
  "Interaction: RC × Number",          0.064,      0.059,  697.621,  1.085,  "0.278"
)


# Apply gt and make the table neat
table_lmm_spillover %>%
  gt() %>%
  tab_header(
    title = "Results of LMM (Region 4: Spillover)"
  ) %>%
  fmt_number(
    columns = c(Estimate, Std_Error, t_value),
    decimals = 3
  ) %>%
  cols_label(
    Term = "Predictor / Fixed effect",
    Estimate = "Estimate",
    Std_Error = "Standard error",
    df = "df",
    t_value = "t value",
    p_value = "p value"
  )

# 1. making gt table based on table before and save it in object
table_gt_spillover <- table_lmm_spillover %>%
  gt() %>%
  tab_header(
    title = "Results of LMM (Region 4: Spillover)"
  ) %>%
  fmt_number(
    columns = c(Estimate, Std_Error, t_value),
    decimals = 3
  ) %>%
  cols_label(
    Term = "Predictor / Fixed effect",
    Estimate = "Estimate",
    Std_Error = "Standard error",
    df = "df",
    t_value = "t value",
    p_value = "p value"
  )

# 2. Save as PNG file
gtsave(table_gt_spillover, filename = "table_lmm_spillover.png")



#---- Analysis: Comprehension Accuracy ----


##---- Install necessary packages ----
#install.packages("tidyverse")  # データ処理用
#install.packages("readr")      # CSV読み込み用
#install.packages("gt")         # 表作成用

##---- Load necessary packages ----
library(tidyverse)
library(readr)
library(ggplot2)
library(dplyr)
library(stringr)

##---- Load the CSV file ----
###---- Set the Names of Columns----
column_names <- c("Result_Reception_Timestamp", 
                  "Participant_IP_Hash", 
                  "Controller_Name", 
                  "Item_Order?", 
                  "Inner_Element_0", 
                  "Label", 
                  "Latin_Square_NULL", 
                  "Penn_Element_Type", 
                  "Penn_Element_Name", 
                  "Process_Chunk_Order", 
                  "Value_Chunk_Content", 
                  "EventTime_RT_Timestamp",
                  "Code", 
                  "List",
                  "Item",
                  "Condition",
                  "Correct_Answer",
                  "MC_Subj_Number",
                  "Type",
                  "Reading_Time",
                  "Comment_0",
                  "Full_Sentence",
                  "Comment"
)

###---- Transform the CSV File ----
# Skipping First 20 lines
# Naming the Columns Based on "column_names"
data <- read_csv("results_final_Comp.csv", skip = 20, col_names = column_names, na = c("NULL", "Wait success", ""))

str(data)
problems(data)

###---- Show the Beginning Part of the Data ----
head(data)

##---- Clean the Mother Data ----
###---- See the Overview of the Data ----
glimpse(data)

###---- Change the Data Type ----
data <- data %>%
  mutate(
    Item_Order = suppressWarnings(as.numeric(`Item_Order?`)),
    Inner_Element_0 = suppressWarnings(as.numeric(Inner_Element_0)),
    EventTime_RT_Timestamp = suppressWarnings(as.numeric(EventTime_RT_Timestamp)),
    Item = suppressWarnings(as.numeric(Item)),
    Condition = suppressWarnings(as.numeric(Condition)),
    Penn_Element_Name = as.character(Penn_Element_Name),
    Correct_Answer = as.character(Correct_Answer),
    MC_Subj_Number = as.character(MC_Subj_Number),
    Type = as.character(Type),
    
    # List の値が 1, 2, 3, 4 の場合のみ numeric に変換
    List = case_when(
      List %in% c("1", "2", "3", "4") ~ suppressWarnings(as.numeric(List)),
      TRUE ~ NA_real_
    ),
    
    # Reading_Time をまず character 型に変換し、NULL や空白を NA に置き換える
    Reading_Time = as.character(Reading_Time), 
    Reading_Time = na_if(Reading_Time, "NULL"),  
    Reading_Time = na_if(Reading_Time, "Wait success"),  
    Reading_Time = na_if(Reading_Time, ""),  
    
    # 数値型に変換（警告を抑制）
    Reading_Time_numeric = suppressWarnings(as.numeric(Reading_Time))
  )

###---- Delete Unnecessary Columns ----
selected_data <- data %>%
  select(-Controller_Name, -Inner_Element_0, -Latin_Square_NULL, -Penn_Element_Type, -Comment_0, -Comment)


glimpse(selected_data)






# cleaned_data_comp を作成
cleaned_data_comp <- selected_data %>%
  filter(Type == "item") %>%  # ターゲット文のみを抽出
  mutate(
    Correct = case_when(
      Process_Chunk_Order == "PressedKey" & Value_Chunk_Content == Correct_Answer ~ 1,
      Process_Chunk_Order == "PressedKey" & Value_Chunk_Content != Correct_Answer ~ 0,
      TRUE ~ NA_real_  # それ以外はNA
    )
  ) %>%
  filter(!is.na(Correct))  # `NA` の行を削除




###---- Extract the Necessary Columns ----
# Only when Process_Chunk_Order == "PressedKey", take the value of Value_Chunk_Content and compare with Correct_Answer
# Type (only extract "item" = target sentences)

# Typeが"item"であるデータに絞り、正答率を計算
accuracy_per_participant <- cleaned_data_comp %>%
  filter(Type == "item") %>%  # Typeが"item" = 実際に分析したいターゲット文に絞る
  group_by(Participant_IP_Hash) %>%  # 被験者ごとにグループ化
  summarise(accuracy = mean(Correct, na.rm = TRUE))  # 被験者ごとの正答率

#print(accuracy_per_participant)

# 被験者の数を確認
num_participants <- nrow(accuracy_per_participant)  

#print(paste("被験者数:", num_participants))

# 被験者が21人であることを確認（パイロットでは19人かどうかを確認したい）
if (num_participants != 21) {
  warning("被験者数が21人ではありません！データを確認してください。")
}

# 全体の記述統計を計算（被験者ごとのデータを使う）
descriptive_accuracy_all <- accuracy_per_participant %>%
  summarise(
    mean_accuracy = mean(accuracy, na.rm = TRUE),
    median_accuracy = median(accuracy, na.rm = TRUE),
    min_accuracy = min(accuracy, na.rm = TRUE),
    max_accuracy = max(accuracy, na.rm = TRUE),
    sd_accuracy = ifelse(num_participants > 1, sd(accuracy, na.rm = TRUE), NA)  # 標準偏差
  )

# 結果を表示
print(descriptive_accuracy_all)


##---- Descriptive Statistics of Comprehension Accuracy by Conditions ----

###---- Descriptive Stats about Comp Analysis: CONDITION ----
# Descriptive Stats of data_clean (Comp Anallysis: CONDITION)
descriptive_accuracy_condition <- cleaned_data_comp %>%
  group_by(Condition) %>%  # Conditionごとにグループ化
  summarise(
    mean_accuracy = mean(Correct, na.rm = TRUE),   # 平均正答率
    median_accuracy = median(Correct, na.rm = TRUE), # 中央値
    min_accuracy = min(Correct, na.rm = TRUE),     # 最小値
    max_accuracy = max(Correct, na.rm = TRUE),     # 最大値
    sd_accuracy = sd(Correct, na.rm = TRUE)        # 標準偏差
  )

# Show Descriptive Stats by Condition
print(descriptive_accuracy_condition)






# Memo: cleaned_data_compにはselected_dataにはない"Correct"という列があったりするData Frame

###---- Descriptive Stats about Comp Analysis: CONDITION ----
# Descriptive Stats of data_clean (Comp Anallysis: CONDITION)
descriptive_accuracy_condition <- cleaned_data_comp %>%
  group_by(Condition) %>%  # Conditionごとにグループ化
  summarise(
    mean_accuracy = mean(Correct, na.rm = TRUE),   # 平均正答率
    median_accuracy = median(Correct, na.rm = TRUE), # 中央値
    min_accuracy = min(Correct, na.rm = TRUE),     # 最小値
    max_accuracy = max(Correct, na.rm = TRUE),     # 最大値
    sd_accuracy = sd(Correct, na.rm = TRUE)        # 標準偏差
  )

# Show Descriptive Stats by Condition
print(descriptive_accuracy_condition)

num_participants <- cleaned_data_comp %>%
  distinct(Participant_IP_Hash) %>%
  nrow()

print(paste("被験者数:", num_participants))


# Visualize the dipersion in participants / Box plot
library(ggplot2)

accuracy_per_participant <- cleaned_data_comp %>%
  group_by(Participant_IP_Hash) %>%
  summarise(accuracy = mean(Correct, na.rm = TRUE))

ggplot(accuracy_per_participant, aes(x = "", y = accuracy)) +
  geom_boxplot(outlier.shape = NA, fill = "lightgray") +
  geom_jitter(width = 0.1, height = 0, size = 2, color = "steelblue") +
  labs(
    title = "Distribution of Comprehension Accuracy by Participant",
    y = "Accuracy (Proportion Correct)",
    x = ""
  ) +
  theme_minimal()

# bar plot
library(ggplot2)
library(scales)

accuracy_per_participant <- cleaned_data_comp %>%
  group_by(Participant_IP_Hash) %>%
  summarise(accuracy = mean(Correct, na.rm = TRUE)) %>%
  arrange(desc(accuracy)) %>%
  mutate(
    Subject_Num = row_number(),
    Subject_ID = paste0("Subject ", Subject_Num)
  )

p <- ggplot(accuracy_per_participant, aes(x = factor(Subject_ID, levels = Subject_ID), y = accuracy)) +
  geom_col(fill = "gray50") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Participant",
    y = "Accuracy"
  ) +
  theme_minimal(base_family = "Times New Roman", base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray80", size = 0.3),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10))
  )

print(p)

ggsave("comprehension_accuracy_by_participant.png", plot = p, width = 6, height = 4, dpi = 300)

# プロット表示
print(p)

# 画像保存（保存場所は作業ディレクトリ）
ggsave("comprehension_accuracy_by_participant.png", plot = p, width = 6, height = 4, dpi = 300)



##---- GLMM Comprehension Accuracy ----
#install.packages("lme4")
#install.packages("lmerTest")

library(lme4)
library(lmerTest)
library(ggplot2)

#English_Proficiencyという列がないので以下この列を新たに作る名前は"English_Proficiency"
names(cleaned_data_comp)

library(dplyr)
library(stringr)

library(dplyr)
library(stringr)
library(lme4)
library(lmerTest)
library(ggplot2)

# 1. proficiency_dataの作成時はdplyr::selectを明示的に使う
proficiency_data <- selected_data %>%
  filter(
    str_detect(Process_Chunk_Order, "English_Proficiency"),
    Value_Chunk_Content == "checked"
  ) %>%
  mutate(
    English_Proficiency = str_replace(Process_Chunk_Order, "English_Proficiency_", "")
  ) %>%
  dplyr::select(Participant_IP_Hash, English_Proficiency) %>%
  distinct()


# 2. cleaned_data_comp に英語力列をマージ
cleaned_data_comp <- cleaned_data_comp %>%
  left_join(proficiency_data, by = "Participant_IP_Hash")

cleaned_data_comp <- cleaned_data_comp %>%
  left_join(proficiency_data, by = "Participant_IP_Hash")

head(cleaned_data_comp$English_Proficiency)

# 3. モデル作成
# GLMM model without interaction
glmm_model_comp_02 <- glmer(
  Correct ~ Condition + English_Proficiency + MC_Subj_Number + 
    (1 | Item),
  data = cleaned_data_comp,
  family = binomial(link = "logit")
)

summary(glmm_model_comp_02)

# GLMM model with interaction
glmm_model_comp_03 <- glmer(
  Correct ~ Condition * English_Proficiency + MC_Subj_Number + 
    (1 | Item),
  data = cleaned_data_comp,
  family = binomial(link = "logit")
)

summary(glmm_model_comp_03)


# proficiency_data が空でないか
nrow(proficiency_data)

# cleaned_data_comp の ID と一致しているか
intersect(cleaned_data_comp$Participant_IP_Hash, proficiency_data$Participant_IP_Hash)




##----- Further Analysis Item ----
library(dplyr)

item_condition_accuracy <- cleaned_data_comp %>%
  group_by(Item, Condition) %>%
  summarise(
    N = n(),
    Accuracy = mean(Correct),
    SD = sd(Correct),
    .groups = "drop"
  ) %>%
  arrange(Item, Condition)

head(item_condition_accuracy, 10) 

# 例：Item 23 のみ表示（視覚的に確認しやすく）
item_condition_accuracy %>% filter(Item == 23)


# Screen the top 5 hardest Items by conditions
item_condition_accuracy %>%
  group_by(Condition) %>%
  slice_min(Accuracy, n = 5, with_ties = FALSE) %>%
  ungroup()

