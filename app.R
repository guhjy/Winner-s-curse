# 載入必要套件
library(shiny)
library(stats) # 用於 power.t.test
library(ggplot2) # 用於繪製圖表
library(dplyr)   # 用於數據整理

# --- 輔助函數：執行單次研究模擬 ---
# 此函數模擬一次研究 (t 檢定) 並返回 p 值和觀察到的 Cohen's d
run_study <- function(true_effect, n, sd = 1) {
  # 確保 n 是有效的整數且至少為 2
  if (is.na(n) || n < 2) {
    return(list(p.value = NA, d_observed = NA))
  }
  n <- round(n) # 確保 n 是整數

  # 模擬兩組數據
  group1 <- rnorm(n, mean = 0, sd = sd)
  # The 'true_effect' here determines the mean difference for data generation
  group2 <- rnorm(n, mean = true_effect, sd = sd)

  # 執行獨立樣本 t 檢定
  # 使用 tryCatch 處理 t.test 可能因常數數據等原因失敗的情況
  t_test_result <- tryCatch({
    t.test(group2, group1, var.equal = TRUE)
  }, error = function(e) {
    # 如果 t.test 失敗，返回 NA p 值
    list(p.value = NA)
  })

  # 如果 t.test 失敗，p.value 會是 NA
  p_value <- t_test_result$p.value
  d_observed <- NA # 預設 d 值為 NA

  # 僅在 t.test 成功且 p 值非 NA 時計算觀察到的 Cohen's d
  if (!is.na(p_value)) {
    mean1 <- mean(group1)
    mean2 <- mean(group2)
    var1 <- var(group1)
    var2 <- var(group2)
    # 計算合併標準差
    # 檢查變異數是否有效
    if (is.na(var1) || is.na(var2) || var1 < 0 || var2 < 0) {
      sd_pooled <- NA
    } else {
      # 避免 n=1 的情況 (雖然前面已檢查 n>=2)
      if (n > 1) {
        # Ensure denominators are not zero before division
        if ((n + n - 2) > 0) {
           sd_pooled <- sqrt(((n - 1) * var1 + (n - 1) * var2) / (n + n - 2))
        } else {
           sd_pooled <- NA # Avoid division by zero if n=1 for both groups somehow
        }
      } else {
        sd_pooled <- NA # 無法計算合併標準差
      }
    }

    # 處理 sd_pooled 為 0 或 NA 或 無限 的罕見情況
    if (is.na(sd_pooled) || sd_pooled == 0 || !is.finite(sd_pooled)) {
      d_observed <- NA
    } else {
      d_observed <- (mean2 - mean1) / sd_pooled
    }
  }

  # 返回結果列表
  return(list(p.value = p_value, d_observed = d_observed))
}


# --- Shiny UI 定義 ---
ui <- fluidPage(
  # 應用程式標題
  titlePanel("獨立樣本 t 檢定：p 值可重複性、統計功效與贏家詛咒 (輸入樣本數，計算功效)"),

  sidebarLayout(
    sidebarPanel(
      h4("模擬參數設定"),
      # *** 修改：輸入每組樣本數 (n)，預設值改為 30 ***
      numericInput("sample_size", "每組樣本數 (n):",
                   value = 30, min = 5, max = 1000, step = 5), # 允許設定 n, 預設值為 30
      # 顯示計算出的功效
      h5("計算出的統計功效 (Power):"),
      strong(textOutput("calculated_power_display", inline = TRUE)),
      hr(), # 分隔線
      # *** 修改：真實效應大小，預設值改為 0.2 ***
      sliderInput("true_effect_HA", "真實效應大小 (Cohen's d, 對立假設情境):",
                  min = 0.05, max = 1.5, value = 0.2, step = 0.05), # 預設值為 0.2
      sliderInput("alpha", "顯著水準 (alpha):",
                  min = 0.001, max = 0.1, value = 0.05, step = 0.001),
      numericInput("num_simulations", "模擬次數:",
                   value = 20000, min = 100, max = 50000, step = 100),
      actionButton("run_sim", "執行模擬", icon = icon("play"), class = "btn-success"),
      hr(),
      p(em("注意：此版本修改了模擬邏輯。重複研究現在基於原始研究觀察到的效應大小進行模擬。"))
    ),

    mainPanel(
      tabsetPanel(
        # 結果與解釋的分頁 (內容大致不變)
        tabPanel("模擬結果與圖表",
                 h4("情境一：對立假設 (H_A) 為真 (存在真實效應)"),
                 verbatimTextOutput("results_HA"),
                 hr(),
                 h4("情境二：虛無假設 (H_0) 為真 (無真實效應)"),
                 verbatimTextOutput("results_H0"),
                 hr(),
                 h4("原始研究 p 值分佈圖 (直方圖)"), # Changed title slightly
                 plotOutput("p_value_distribution_plot")
        ),
        tabPanel("分層重複成功率比較",
                 h4("情境一：對立假設 (H_A) 為真 (存在真實效應)"),
                 p("比較不同初始 p 值區間的重複成功率 (重複研究 p < alpha 的比例)："),
                 verbatimTextOutput("stratified_results_HA"),
                 hr(),
                 h4("情境二：虛無假設 (H_0) 為真 (無真實效應)"),
                 p("比較不同初始偽陽性 p 值區間的重複成功率 (重複研究 p < alpha 的比例，即偽陽性重複率)："),
                 verbatimTextOutput("stratified_results_H0")
        ),
        tabPanel("概念解釋",
                 h4("關於此模擬"),
                 p("這個應用程式模擬了成對的「原始研究」與「重複研究」，並探討兩種主要情境："),
                 tags$ul(
                   tags$li(strong("情境一 (對立假設 H_A 為真)："), "假設母體中存在一個真實的效應，其大小由側邊欄的「真實效應大小」參數決定。"),
                   tags$li(strong("情境二 (虛無假設 H_0 為真)："), "假設母體中不存在真實的效應（真實效應大小為零）。")
                 ),
                 p("您可以調整側邊欄的參數（現在包含「每組樣本數」），然後點擊「執行模擬」按鈕，觀察參數變化對模擬結果（包括計算出的功效）的影響。"),
                 p(strong("重要修改："), "在此版本中，每個「重複研究」的數據生成是基於其配對的「原始研究」所 *觀察到* 的效應大小 (Cohen's d)，而不是基於預設的真實效應大小。這模擬了一種情況，即研究者根據初步發現來設計或預期重複研究的結果。如果原始研究的效應量無法計算，則退回使用預設的真實效應大小。"),

                 h4("關鍵概念與結果解讀"),
                 strong("統計功效 (Statistical Power)："),
                 # *** 修改：描述反映計算出的功效 ***
                 p("指當對立假設 (H_A) 為真（即真實效應存在）時，研究能夠成功偵測到此效應（也就是得到 p < alpha）的機率。功效越高越好。它主要取決於真實效應大小、樣本數、alpha 水準以及數據變異性。在此應用程式中，您設定「每組樣本數」，我們會根據其他參數計算出預期的「統計功效」。"),
                 p(em("在情境一 (對立假設 H_A 為真) 中："), "模擬結果會估計實際的功效（顯示為「原始研究顯著的比例」），這個值應接近側邊欄顯示的「計算出的統計功效」。"),

                 # 其他解釋大致相同，但對 "目標功效" 或 "計算出的 n" 的引用已隱含更改
                 strong("重複成功率 (Replication Rate)："),
                 p("指一個報告了顯著結果 (p < alpha) 的原始研究，在相同條件下被重複進行時，能夠再次得到顯著結果的機率。"),
                 p(em("在此模擬的特定設置下："), "重複研究的成功與否現在受到原始研究觀察到的效應大小（可能因隨機誤差而偏離真實值）的影響。"),
                 p(em("在情境一 (對立假設 H_A 為真) 中："), "重複成功率會受到原始研究的功效以及觀察到的效應量（贏家詛咒可能使其膨脹）的雙重影響。"),
                 p(em("在情境二 (虛無假設 H_0 為真) 中："), "原始研究的「顯著」結果是偽陽性。重複研究現在基於這個（因隨機性產生的）非零觀察效應量進行模擬。重複研究再次得到「顯著」結果的機率可能會高於 alpha 水準，因為它是基於一個非零的（儘管是隨機的）效應預期。"),
                 p(em("在「分層重複成功率比較」分頁中："), "觀察初始 p 值與重複成功率的關係。在 H_A 下，關係可能更複雜。在 H_0 下，由於重複研究基於原始的隨機效應，重複率可能不再恆等於 alpha，且可能與原始 p 值（或其對應的觀察效應量）相關。"),

                 strong("贏家詛咒 (Winner's Curse)："),
                 p("指那些因為達到統計顯著而被挑選出來的研究（「贏家」），其觀察到的效應量往往會高估真實效應量的現象。這是因為，除了真實效應外，隨機誤差也「幫助」了這些研究達到顯著，而這種「有幫助」的隨機誤差通常會使觀察到的效應看起來更大。"),
                 p(em("在情境一 (對立假設 H_A 為真) 中："), "比較「所有 H_A 研究」和「*僅限顯著的* H_A 研究」的平均觀察效應量 (d)。後者（「贏家」）的平均 d 通常會更高，這就體現了贏家詛咒。"),
                 p(em("在情境二 (虛無假設 H_0 為真) 中："), "類似現象也會發生。偽陽性的研究（「贏家」）必然因為隨機誤差而觀察到一個偏離真實值零的效應量。我們關注的是「平均*絕對*觀察效應量 |d|」，對於偽陽性研究，這個值會顯著大於零。"),

                 strong("向平均值迴歸 (Regression Toward the Mean)："),
                 p("當您重複一個初次觀察到極端結果（例如與顯著性相關的非常大的效應量）的研究時，重複研究所得到的結果，其效應量大小更可能接近母體的真實平均值（通常不如第一次那麼極端）。"),
                 p(em("在情境一 (對立假設 H_A 為真) 中："), "比較原始「贏家」研究的平均 d，與它們對應的「重複研究」的平均 d。由於重複研究現在基於原始觀察到的 d，這個現象的表現可能會與原始版本不同，但總體趨勢可能仍然存在，因為單次的隨機高估不太可能在重複中持續。"),
                 p(em("在情境二 (虛無假設 H_0 為真) 中："), "比較原始偽陽性研究的平均*絕對* |d|，與它們對應的重複研究的平均*絕對* |d|。重複研究的平均 |d| 通常較小，更接近真實值零。"),

                 strong("P 值分佈圖："),
                 p("此圖顯示了在對立假設為真和虛無假設為真兩種情境下，模擬出的原始研究 p 值的分布情況。"),
                 p(em("當虛無假設 (H_0) 為真時："), "理論上 p 值應呈現均勻分佈（圖上的直方圖或密度曲線應大致是平坦的）。低於 alpha 的 p 值比例應接近 alpha。使用直方圖時，每個長條的高度應大致相等。"), # Added note for histogram
                 p(em("當對立假設 (H_A) 為真時："), "p 值會傾向於集中在較低的值。低於 alpha 的 p 值比例反映了研究的統計功效。分佈的形狀取決於功效的大小。")
        )
      )
    )
  )
)

# --- Shiny Server 邏輯 ---
server <- function(input, output, session) {

  # --- 響應式表達式：根據輸入的樣本數計算功效 ---
  calculated_power <- reactive({
    # 確保效應大小為正且 n 有效
    req(input$true_effect_HA > 0, input$sample_size >= 2)
    # 使用 tryCatch 處理可能的錯誤
    power_result <- tryCatch({
      power.t.test(
        n = input$sample_size, # 使用輸入的樣本數
        delta = input$true_effect_HA,
        sd = 1,
        sig.level = input$alpha,
        type = "two.sample",
        alternative = "two.sided",
        power = NULL # 讓函數解出 power
      )$power # 提取計算出的 power
    }, error = function(e) {
      NA # 如果出錯，返回 NA
    })

    return(power_result)
  })

  # --- 輸出：顯示計算出的功效 ---
  output$calculated_power_display <- renderText({
    power_val <- calculated_power()
    if (is.na(power_val)) {
      "無法計算 (請檢查參數)"
    } else {
      sprintf("%.1f%%", power_val * 100) # 格式化為百分比
    }
  })

  # --- 響應式表達式：執行模擬 (當按鈕被點擊時) ---
  simulation_results <- eventReactive(input$run_sim, {
    # *** 修改：直接從輸入獲取樣本數 ***
    n_to_use <- input$sample_size
    # 執行前確保 n 有效
    req(!is.na(n_to_use), n_to_use >= 2)

    # 顯示通知
    showNotification("正在執行模擬...", type = "message", duration = 3)

    # 獲取其他參數
    true_effect_HA <- input$true_effect_HA
    alpha <- input$alpha
    num_sim <- input$num_simulations
    true_effect_H0 <- 0
    sd_val <- 1
    calculated_pow <- calculated_power() # 獲取計算出的功效以供顯示

    # --- 執行成對模擬 (H_A 和 H_0) ---
    # 初始化數據框
    results_HA_pairs <- data.frame(p_original = numeric(num_sim), d_original = numeric(num_sim),
                                   p_replication = numeric(num_sim), d_replication = numeric(num_sim))
    results_H0_pairs <- data.frame(p_original = numeric(num_sim), d_original = numeric(num_sim),
                                   p_replication = numeric(num_sim), d_replication = numeric(num_sim))

    # 使用進度條
    withProgress(message = '模擬進行中', value = 0, {
      for (i in 1:num_sim) {
        # --- H_A 模擬對 ---
        # 1. 執行原始研究
        res_orig_HA <- run_study(true_effect = true_effect_HA, n = n_to_use, sd = sd_val)

        # 2. *** MODIFICATION: Determine effect size for replication ***
        #    Use the observed effect size (d_observed) from the original study.
        #    Handle NA case: if d_observed is NA, fall back to the input true_effect_HA.
        effect_for_rep_HA <- if (!is.na(res_orig_HA$d_observed) && is.finite(res_orig_HA$d_observed)) {
          res_orig_HA$d_observed
        } else {
          true_effect_HA # Fallback if d_observed is NA or infinite
        }

        # 3. 執行重複研究，使用從原始研究確定的效應大小
        res_rep_HA <- run_study(true_effect = effect_for_rep_HA, n = n_to_use, sd = sd_val)

        # 4. 儲存結果
        results_HA_pairs$p_original[i] <- res_orig_HA$p.value
        results_HA_pairs$d_original[i] <- res_orig_HA$d_observed
        results_HA_pairs$p_replication[i] <- res_rep_HA$p.value
        results_HA_pairs$d_replication[i] <- res_rep_HA$d_observed

        # --- H_0 模擬對 ---
        # 1. 執行原始研究 (真實效應為 0)
        res_orig_H0 <- run_study(true_effect = true_effect_H0, n = n_to_use, sd = sd_val)

        # 2. *** MODIFICATION: Determine effect size for replication ***
        #    Use the observed effect size (d_observed) from the original study.
        #    Handle NA case: if d_observed is NA, fall back to the true_effect_H0 (which is 0).
        effect_for_rep_H0 <- if (!is.na(res_orig_H0$d_observed) && is.finite(res_orig_H0$d_observed)) {
          res_orig_H0$d_observed
        } else {
          true_effect_H0 # Fallback if d_observed is NA or infinite (effectively 0)
        }

        # 3. 執行重複研究，使用從原始研究確定的效應大小
        res_rep_H0 <- run_study(true_effect = effect_for_rep_H0, n = n_to_use, sd = sd_val)

        # 4. 儲存結果
        results_H0_pairs$p_original[i] <- res_orig_H0$p.value
        results_H0_pairs$d_original[i] <- res_orig_H0$d_observed
        results_H0_pairs$p_replication[i] <- res_rep_H0$p.value
        results_H0_pairs$d_replication[i] <- res_rep_H0$d_observed

        # 更新進度
        incProgress(1/num_sim, detail = paste("完成", i, "/", num_sim))
      }
    }) # 結束 withProgress

    # 移除模擬中可能產生的 NA 值 (主要來自 p 值)
    results_HA_pairs <- results_HA_pairs[complete.cases(results_HA_pairs$p_original, results_HA_pairs$p_replication), ]
    results_H0_pairs <- results_H0_pairs[complete.cases(results_H0_pairs$p_original, results_H0_pairs$p_replication), ]
    # d_original 或 d_replication 可能仍有 NA，如果計算失敗，但不影響 p 值分析
    num_valid_sim_HA <- nrow(results_HA_pairs)
    num_valid_sim_H0 <- nrow(results_H0_pairs)


    # --- 分析 H_A 結果 ---
    significant_originals_HA <- results_HA_pairs[results_HA_pairs$p_original < alpha, ]
    num_significant_HA <- nrow(significant_originals_HA)
    output_HA_text <- ""
    output_stratified_HA <- "無顯著原始研究可供分層分析。\n"

    if (num_valid_sim_HA > 0) {
      prop_sig_HA <- mean(results_HA_pairs$p_original < alpha) * 100
      output_HA_text <- paste0(output_HA_text, sprintf("原始研究顯著的比例：%.2f%% (經驗功效，基於 %d 次有效模擬)\n", prop_sig_HA, num_valid_sim_HA))
      # *** 修改：顯示計算出的功效 ***
      if (!is.na(calculated_pow)) {
         output_HA_text <- paste0(output_HA_text, sprintf("  (理論計算功效 = %.1f%%)\n", calculated_pow * 100))
      } else {
         output_HA_text <- paste0(output_HA_text, "  (無法計算理論功效)\n")
      }
    } else {
      output_HA_text <- "無有效的對立假設模擬結果。\n"
    }

    if (num_significant_HA > 0) {
      # 整體重複率
      rate_rep_HA <- mean(significant_originals_HA$p_replication < alpha) * 100
      output_HA_text <- paste0(output_HA_text, sprintf("整體重複成功率 (給定原始研究顯著)：%.2f%%\n", rate_rep_HA))
      output_HA_text <- paste0(output_HA_text, "  (注意：重複研究基於原始觀察效應量模擬)\n")


      # 贏家詛咒 (分析原始研究的 d)
      # Filter out potential NA/Inf values in d_original before calculating mean
      avg_d_all_HA <- mean(results_HA_pairs$d_original, na.rm = TRUE)
      avg_d_sig_HA <- mean(significant_originals_HA$d_original, na.rm = TRUE)
      output_HA_text <- paste0(output_HA_text, "\n--- 贏家詛咒分析 (對立假設 H_A - 原始研究) ---\n")
      output_HA_text <- paste0(output_HA_text, sprintf("平均觀察效應量 (d) - 所有研究：%.3f\n", avg_d_all_HA))
      output_HA_text <- paste0(output_HA_text, sprintf("平均觀察效應量 (d) - *僅限顯著* 研究 ('贏家')：%.3f\n", avg_d_sig_HA))
      if(!is.na(avg_d_sig_HA) && !is.na(avg_d_all_HA) && avg_d_sig_HA > avg_d_all_HA) output_HA_text <- paste0(output_HA_text, "  -> 顯著研究的平均效應量被高估。\n")

      # 向平均值迴歸 (比較原始顯著研究的 d 和它們對應的重複研究的 d)
      # Filter out potential NA/Inf values before calculating mean
      avg_d_rep_for_winners <- mean(significant_originals_HA$d_replication, na.rm = TRUE)
      output_HA_text <- paste0(output_HA_text, "\n--- 向平均值迴歸分析 (對立假設 H_A) ---\n")
      output_HA_text <- paste0(output_HA_text, sprintf("原始'贏家'研究的平均觀察效應量 (d)：%.3f\n", avg_d_sig_HA))
      output_HA_text <- paste0(output_HA_text, sprintf("對應重複研究的平均觀察效應量 (d)：%.3f\n", avg_d_rep_for_winners))
      if(!is.na(avg_d_rep_for_winners) && !is.na(avg_d_sig_HA) && avg_d_rep_for_winners < avg_d_sig_HA) output_HA_text <- paste0(output_HA_text, "  -> 重複研究效應量向平均值迴歸。\n")

      # 分層重複成功率計算 (H_A)
      bin1_HA <- significant_originals_HA[significant_originals_HA$p_original >= 0.01 & significant_originals_HA$p_original < 0.05, ]
      bin2_HA <- significant_originals_HA[significant_originals_HA$p_original >= 0.001 & significant_originals_HA$p_original < 0.01, ]
      bin3_HA <- significant_originals_HA[significant_originals_HA$p_original < 0.001, ]

      output_stratified_HA <- ""
      if (nrow(bin1_HA) > 0) {
        rate1 <- mean(bin1_HA$p_replication < alpha) * 100
        output_stratified_HA <- paste0(output_stratified_HA, sprintf(" - 原始 p 值介於 [0.01, 0.05) (共 %d 次): 重複成功率 = %.2f%%\n", nrow(bin1_HA), rate1))
      } else { output_stratified_HA <- paste0(output_stratified_HA, " - 原始 p 值介於 [0.01, 0.05): 無此情況。\n") }
      if (nrow(bin2_HA) > 0) {
        rate2 <- mean(bin2_HA$p_replication < alpha) * 100
        output_stratified_HA <- paste0(output_stratified_HA, sprintf(" - 原始 p 值介於 [0.001, 0.01) (共 %d 次): 重複成功率 = %.2f%%\n", nrow(bin2_HA), rate2))
      } else { output_stratified_HA <- paste0(output_stratified_HA, " - 原始 p 值介於 [0.001, 0.01): 無此情況。\n") }
      if (nrow(bin3_HA) > 0) {
        rate3 <- mean(bin3_HA$p_replication < alpha) * 100
        output_stratified_HA <- paste0(output_stratified_HA, sprintf(" - 原始 p 值 < 0.001 (共 %d 次): 重複成功率 = %.2f%%\n", nrow(bin3_HA), rate3))
      } else { output_stratified_HA <- paste0(output_stratified_HA, " - 原始 p 值 < 0.001: 無此情況。\n") }

    } else {
      output_HA_text <- paste0(output_HA_text, "\n在對立假設情境下，無顯著的原始研究可供分析重複性/贏家詛咒。\n")
    }

    # --- 分析 H_0 結果 ---
    significant_originals_H0 <- results_H0_pairs[results_H0_pairs$p_original < alpha, ]
    num_significant_H0 <- nrow(significant_originals_H0)
    output_H0_text <- ""
    output_stratified_H0 <- "無偽陽性原始研究可供分層分析。\n"

    if (num_valid_sim_H0 > 0) {
      prop_sig_H0 <- mean(results_H0_pairs$p_original < alpha) * 100
      output_H0_text <- paste0(output_H0_text, sprintf("原始研究顯著的比例 (偽陽性)：%.2f%% (經驗第一類型錯誤率，基於 %d 次有效模擬)\n", prop_sig_H0, num_valid_sim_H0))
      output_H0_text <- paste0(output_H0_text, sprintf("  (應接近 alpha = %.1f%%)\n", alpha * 100))
    } else {
      output_H0_text <- "無有效的虛無假設模擬結果。\n"
    }

    if (num_significant_H0 > 0) {
      # 偽陽性重複率
      rate_rep_H0 <- mean(significant_originals_H0$p_replication < alpha) * 100
      output_H0_text <- paste0(output_H0_text, sprintf("整體偽陽性重複率 (給定原始研究為偽陽性)：%.2f%%\n", rate_rep_H0))
      output_H0_text <- paste0(output_H0_text, "  (注意：重複研究基於原始偽陽性觀察效應量模擬)\n")
      # output_H0_text <- paste0(output_H0_text, sprintf("  (可能不再等於 alpha = %.1f%%)\n", alpha * 100)) # Commented out as it might be confusing

      # H0 下的贏家詛咒現象 (分析原始研究的 d)
      # Filter out potential NA/Inf values before calculating mean
      avg_d_all_H0 <- mean(results_H0_pairs$d_original, na.rm = TRUE)
      avg_abs_d_sig_H0 <- mean(abs(significant_originals_H0$d_original), na.rm = TRUE)
      output_H0_text <- paste0(output_H0_text, "\n--- 虛無假設 (H_0) 情境下的贏家詛咒類似現象分析 ---\n")
      output_H0_text <- paste0(output_H0_text, sprintf("平均觀察效應量 (d) - 所有研究：%.3f (應接近 0)\n", avg_d_all_H0))
      output_H0_text <- paste0(output_H0_text, sprintf("平均 *絕對* 效應量 |d| - *僅限偽陽性* 研究 ('贏家')：%.3f\n", avg_abs_d_sig_H0))
      if(!is.na(avg_abs_d_sig_H0) && !is.na(avg_d_all_H0) && avg_abs_d_sig_H0 > abs(avg_d_all_H0) + 0.01) output_H0_text <- paste0(output_H0_text, "  -> 偽陽性結果因隨機誤差顯示非零效應量。\n")

      # H0 下的向平均值迴歸 (比較原始偽陽性研究的 |d| 和它們對應的重複研究的 |d|)
      # Filter out potential NA/Inf values before calculating mean
      avg_abs_d_rep_for_winners <- mean(abs(significant_originals_H0$d_replication), na.rm = TRUE)
      output_H0_text <- paste0(output_H0_text, "\n--- 虛無假設 (H_0) 情境下的向平均值迴歸分析 ---\n")
      output_H0_text <- paste0(output_H0_text, sprintf("原始偽陽性研究的平均 *絕對* 效應量 |d|：%.3f\n", avg_abs_d_sig_H0))
      output_H0_text <- paste0(output_H0_text, sprintf("對應重複研究的平均 *絕對* 效應量 |d|：%.3f\n", avg_abs_d_rep_for_winners))
      if(!is.na(avg_abs_d_rep_for_winners) && !is.na(avg_abs_d_sig_H0) && avg_abs_d_rep_for_winners < avg_abs_d_sig_H0) output_H0_text <- paste0(output_H0_text, "  -> 重複研究的絕對效應量向零迴歸。\n")

      # 分層偽陽性重複率計算 (H_0)
      bin1_H0 <- significant_originals_H0[significant_originals_H0$p_original >= 0.01 & significant_originals_H0$p_original < 0.05, ]
      bin2_H0 <- significant_originals_H0[significant_originals_H0$p_original >= 0.001 & significant_originals_H0$p_original < 0.01, ]
      bin3_H0 <- significant_originals_H0[significant_originals_H0$p_original < 0.001, ]

      output_stratified_H0 <- ""
      if (nrow(bin1_H0) > 0) {
        rate1_H0 <- mean(bin1_H0$p_replication < alpha) * 100
        output_stratified_H0 <- paste0(output_stratified_H0, sprintf(" - 原始 p 值介於 [0.01, 0.05) (共 %d 次): 偽陽性重複率 = %.2f%%\n", nrow(bin1_H0), rate1_H0))
      } else { output_stratified_H0 <- paste0(output_stratified_H0, " - 原始 p 值介於 [0.01, 0.05): 無此偽陽性情況。\n") }
      if (nrow(bin2_H0) > 0) {
        rate2_H0 <- mean(bin2_H0$p_replication < alpha) * 100
        output_stratified_H0 <- paste0(output_stratified_H0, sprintf(" - 原始 p 值介於 [0.001, 0.01) (共 %d 次): 偽陽性重複率 = %.2f%%\n", nrow(bin2_H0), rate2_H0))
      } else { output_stratified_H0 <- paste0(output_stratified_H0, " - 原始 p 值介於 [0.001, 0.01): 無此偽陽性情況。\n") }
      if (nrow(bin3_H0) > 0) {
        rate3_H0 <- mean(bin3_H0$p_replication < alpha) * 100
        output_stratified_H0 <- paste0(output_stratified_H0, sprintf(" - 原始 p 值 < 0.001 (共 %d 次): 偽陽性重複率 = %.2f%%\n", nrow(bin3_H0), rate3_H0))
      } else { output_stratified_H0 <- paste0(output_stratified_H0, " - 原始 p 值 < 0.001: 無此偽陽性情況。\n") }
      output_stratified_H0 <- paste0(output_stratified_H0, sprintf("   (注意：在 H0 下，重複率基於原始觀察效應量，可能不再等於 alpha = %.1f%%)\n", alpha * 100))

    } else {
      output_H0_text <- paste0(output_H0_text, "\n在虛無假設情境下，無顯著的（偽陽性）原始研究。\n")
    }

    # 返回包含所有結果的列表
    list(
      text_HA = output_HA_text,
      text_H0 = output_H0_text,
      stratified_HA = output_stratified_HA,
      stratified_H0 = output_stratified_H0,
      p_values_HA = results_HA_pairs$p_original, # Keep original p-values for distribution plot
      p_values_H0 = results_H0_pairs$p_original  # Keep original p-values for distribution plot
    )
  })

  # --- 輸出：渲染文字結果 ---
  output$results_HA <- renderText({ simulation_results()$text_HA })
  output$results_H0 <- renderText({ simulation_results()$text_H0 })
  output$stratified_results_HA <- renderText({ simulation_results()$stratified_HA })
  output$stratified_results_H0 <- renderText({ simulation_results()$stratified_H0 })

  # --- 輸出：渲染 p 值分佈圖 ---
  output$p_value_distribution_plot <- renderPlot({
    res <- simulation_results() # 獲取模擬結果
    # 確保 p 值數據可用 (使用 req 函數)
    req(length(res$p_values_HA) > 0 || length(res$p_values_H0) > 0)
    alpha_val <- input$alpha # 獲取當前的 alpha 值

    # 準備繪圖數據
    plot_data_list <- list()
    if(length(res$p_values_HA) > 0) {
        plot_data_list$HA <- data.frame(p_value = res$p_values_HA, Hypothesis = "對立假設 (H_A) 為真")
    }
    if(length(res$p_values_H0) > 0) {
        plot_data_list$H0 <- data.frame(p_value = res$p_values_H0, Hypothesis = "虛無假設 (H_0) 為真")
    }
    # Check if plot_data_list is empty before binding
    if (length(plot_data_list) == 0) {
        return(NULL) # Return NULL if no data to plot
    }
    plot_data <- bind_rows(plot_data_list)

    # Check if plot_data is empty or has no valid p_values after binding
    if (nrow(plot_data) == 0 || all(is.na(plot_data$p_value))) {
        return(NULL) # Return NULL if no valid data
    }

    # 確保 Hypothesis 是因子類型，並設定正確的層級順序 (如果需要)
    present_levels <- intersect(c("對立假設 (H_A) 為真", "虛無假設 (H_0) 為真"), unique(plot_data$Hypothesis))
    if (length(present_levels) > 0) {
      plot_data$Hypothesis <- factor(plot_data$Hypothesis, levels = present_levels)
    } else {
      warning("Hypothesis factor levels could not be set properly.")
    }

    # Define colors making sure names match the factor levels present in the data
    color_map <- c("對立假設 (H_A) 為真" = "skyblue", "虛無假設 (H_0) 為真" = "salmon")
    active_colors <- color_map[present_levels] # Use only colors for levels present


    # *** 使用 ggplot2 繪製直方圖 (顯示密度) ***
    ggplot(plot_data, aes(x = p_value, fill = Hypothesis)) +
      # Use geom_histogram, map y to density, set binwidth and boundary
      # Use position="identity" to overlay histograms with transparency
      geom_histogram(aes(y = ..density..), binwidth = 0.05, boundary = 0, alpha = 0.7, position = "identity", na.rm = TRUE) +
      # Add vline for alpha level
      geom_vline(xintercept = alpha_val, linetype = "dashed", color = "red", linewidth = 1) +
      # Add annotation for alpha
      annotate("text", x = alpha_val, y = Inf, label = paste("alpha =", alpha_val), hjust = -0.1, vjust = 1.5, color = "red", size = 4) +
      # Use the dynamically created color map and limits
      scale_fill_manual(values = active_colors,
                        limits = present_levels,
                        drop = FALSE) +
      labs(
        title = "原始研究 P 值分佈 (直方圖)", # Updated title
        x = "P 值",
        y = "密度", # Y-axis is density
        fill = "假設情境"
      ) +
      theme_minimal(base_family = "sans") +
      theme(
        plot.title = element_text(hjust = 0.5, size=16),
        legend.position = "bottom"
      ) +
      coord_cartesian(xlim = c(0, 1)) # 限制 x 軸範圍為 0 到 1

  })

}

# --- 執行 Shiny 應用程式 ---
shinyApp(ui = ui, server = server)
