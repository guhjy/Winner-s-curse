# 載入必要套件
library(shiny)
library(stats) # 用於 power.t.test
library(ggplot2) # 用於繪製更美觀的圖表
library(dplyr)   # 用於數據整理

# --- 輔助函數：執行單次研究模擬 ---
# 此函數模擬一次研究 (t 檢定) 並返回 p 值和觀察到的 Cohen's d
run_study <- function(true_effect, n, sd = 1) {
  # 模擬兩組數據
  group1 <- rnorm(n, mean = 0, sd = sd)
  group2 <- rnorm(n, mean = true_effect, sd = sd)

  # 執行獨立樣本 t 檢定
  t_test_result <- t.test(group2, group1, var.equal = TRUE)

  # 計算觀察到的 Cohen's d
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  var1 <- var(group1)
  var2 <- var(group2)
  # 計算合併標準差
  sd_pooled <- sqrt(((n - 1) * var1 + (n - 1) * var2) / (n + n - 2))
  # 處理 sd_pooled 為 0 或 NA 的罕見情況
  if (sd_pooled == 0 || is.na(sd_pooled) || !is.finite(sd_pooled)) {
    d_observed <- NA
  } else {
    d_observed <- (mean2 - mean1) / sd_pooled
  }

  # 返回結果列表
  return(list(p.value = t_test_result$p.value, d_observed = d_observed))
}


# --- Shiny UI 定義 ---
ui <- fluidPage(
  titlePanel("模擬實驗：研究重複性、統計功效與贏家詛咒"),

  sidebarLayout(
    sidebarPanel(
      h4("模擬參數設定"),
      sliderInput("n_per_group", "每組樣本數 (n):",
                  min = 10, max = 200, value = 30, step = 5),
      sliderInput("true_effect_HA", "真實效應大小 (Cohen's d, H_A 情境):",
                  min = 0, max = 1.5, value = 0.5, step = 0.05),
      sliderInput("alpha", "顯著水準 (alpha):",
                  min = 0.001, max = 0.1, value = 0.05, step = 0.001),
      numericInput("num_simulations", "模擬次數:",
                   value = 2000, min = 100, max = 20000, step = 100),
      actionButton("run_sim", "執行模擬", icon = icon("play"), class = "btn-success"),
      hr(),
      h4("理論統計功效 (H_A 情境)"),
      # 顯示基於輸入參數計算的理論功效
      textOutput("theoretical_power_display")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("模擬結果與圖表",
                 h4("情境一：H_A 為真 (存在真實效應)"),
                 # 顯示 H_A 情境的文字結果
                 verbatimTextOutput("results_HA"),
                 hr(),
                 h4("情境二：H_0 為真 (無真實效應)"),
                 # 顯示 H_0 情境的文字結果
                 verbatimTextOutput("results_H0"),
                 hr(),
                 h4("原始研究 p 值分佈圖"),
                 # 顯示 p 值分佈圖
                 plotOutput("p_value_distribution_plot")
        ),
        tabPanel("概念解釋",
                 h4("關於此模擬"),
                 p("這個應用程式模擬了成對的「原始研究」與「重複研究」，並探討兩種主要情境："),
                 tags$ul(
                   tags$li(strong("情境一 (H_A 為真)："), "假設母體中存在一個真實的效應，其大小由側邊欄的「真實效應大小」參數決定。"),
                   tags$li(strong("情境二 (H_0 為真)："), "假設母體中不存在真實的效應（真實效應大小為零）。")
                 ),
                 p("您可以調整側邊欄的參數，然後點擊「執行模擬」按鈕，觀察參數變化對模擬結果的影響。"),

                 h4("關鍵概念與結果解讀"),
                 strong("統計功效 (Statistical Power)："),
                 p("指當備擇假設 (H_A) 為真（即真實效應存在）時，研究能夠成功偵測到此效應（也就是得到 p < alpha）的機率。功效越高越好。它主要取決於真實效應大小、樣本數、alpha 水準以及數據變異性。側邊欄會根據您設定的 H_A 參數，顯示計算出的*理論*統計功效。"),
                 p(em("在情境一 (H_A 為真) 中："), "模擬結果會估計實際的功效（顯示為「原始研究顯著的比例」）。對於那些原始研究顯著的案例，「整體重複成功率」也應該接近研究本身的功效（假設重複研究設計相同）。"),

                 strong("重複成功率 (Replication Rate)："),
                 p("指一個報告了顯著結果 (p < alpha) 的原始研究，在相同條件下被重複進行時，能夠再次得到顯著結果的機率。"),
                 p(em("在情境一 (H_A 為真) 中："), "這個成功率主要由研究的統計功效決定。"),
                 p(em("在情境二 (H_0 為真) 中："), "此時，原始研究的「顯著」結果其實是偽陽性（第一類型錯誤）。重複研究再次得到「顯著」結果（也就是另一個偽陽性）的機率，就等於您設定的 alpha 水準。"),

                 strong("贏家詛咒 (Winner's Curse)："),
                 p("指那些因為達到統計顯著而被挑選出來的研究（「贏家」），其觀察到的效應量往往會高估真實效應量的現象。這是因為，除了真實效應外，隨機誤差也「幫助」了這些研究達到顯著，而這種「有幫助」的隨機誤差通常會使觀察到的效應看起來更大。"),
                 p(em("在情境一 (H_A 為真) 中："), "比較「所有 H_A 研究」和「*僅限顯著的* H_A 研究」的平均觀察效應量 (d)。後者（「贏家」）的平均 d 通常會更高，這就體現了贏家詛咒。"),
                 p(em("在情境二 (H_0 為真) 中："), "類似現象也會發生。偽陽性的研究（「贏家」）必然因為隨機誤差而觀察到一個偏離真實值零的效應量。我們關注的是「平均*絕對*觀察效應量 |d|」，對於偽陽性研究，這個值會顯著大於零。"),

                 strong("趨中迴歸 (Regression Toward the Mean)："),
                 p("當您重複一個初次觀察到極端結果（例如與顯著性相關的非常大的效應量）的研究時，重複研究所得到的結果，其效應量大小更可能接近母體的真實平均值（通常不如第一次那麼極端）。"),
                 p(em("在情境一 (H_A 為真) 中："), "比較原始「贏家」研究的平均 d，與它們對應的「重複研究」的平均 d。重複研究的平均 d 通常較小，更接近真實效應大小。"),
                 p(em("在情境二 (H_0 為真) 中："), "比較原始偽陽性研究的平均*絕對* |d|，與它們對應的重複研究的平均*絕對* |d|。重複研究的平均 |d| 通常較小，更接近真實值零。"),

                 strong("P 值分佈圖："),
                 p("此圖顯示了在 H_A 為真和 H_0 為真兩種情境下，模擬出的原始研究 p 值的分布情況。"),
                 p(em("當 H_0 為真時："), "理論上 p 值應呈現均勻分佈（圖上的直方圖或密度曲線應大致是平坦的）。低於 alpha 的 p 值比例應接近 alpha。"),
                 p(em("當 H_A 為真時："), "p 值會傾向於集中在較低的值。低於 alpha 的 p 值比例反映了研究的統計功效。分佈的形狀取決於功效的大小。")
        )
      )
    )
  )
)

# --- Shiny Server 邏輯 ---
server <- function(input, output, session) {

  # --- 響應式表達式：計算理論功效 ---
  theoretical_power <- reactive({
    # 使用 tryCatch 處理可能的計算錯誤 (例如 n 過小)
    tryCatch({
      power_result <- power.t.test(
        n = input$n_per_group,
        delta = input$true_effect_HA,
        sd = 1, # 假設 sd=1，與模擬一致
        sig.level = input$alpha,
        type = "two.sample",
        alternative = "two.sided"
      )
      # 返回格式化的功效百分比
      sprintf("%.1f%%", power_result$power * 100)
    }, error = function(e) {
      "計算錯誤" # 若計算失敗則返回錯誤訊息
    })
  })

  # --- 輸出：顯示理論功效 ---
  output$theoretical_power_display <- renderText({
    theoretical_power()
  })

  # --- 響應式表達式：執行模擬 (當按鈕被點擊時) ---
  simulation_results <- eventReactive(input$run_sim, {
    # 顯示模擬進行中的通知
    showNotification("正在執行模擬...", type = "message", duration = 3)

    # 獲取使用者輸入的參數
    n <- input$n_per_group
    true_effect_HA <- input$true_effect_HA
    alpha <- input$alpha
    num_sim <- input$num_simulations
    true_effect_H0 <- 0
    sd_val <- 1 # 與功效計算一致

    # --- 執行 H_A 模擬 ---
    results_HA_list <- replicate(num_sim, run_study(true_effect = true_effect_HA, n = n, sd = sd_val), simplify = FALSE)
    # --- 執行 H_0 模擬 ---
    results_H0_list <- replicate(num_sim, run_study(true_effect = true_effect_H0, n = n, sd = sd_val), simplify = FALSE)

    # --- 整理 H_A 結果 ---
    results_HA <- data.frame(
      p_original = sapply(results_HA_list, `[[`, "p.value"),
      d_original = sapply(results_HA_list, `[[`, "d_observed")
      # 我們只需要原始研究的 p 值和 d 值來進行後續分析和繪圖
      # 如果需要重複研究的數據，可以類似地提取
    )
    # 為了計算重複率和迴歸平均值，我們需要模擬成對的研究
    results_HA_pairs <- data.frame(
        p_original = numeric(num_sim), d_original = numeric(num_sim),
        p_replication = numeric(num_sim), d_replication = numeric(num_sim)
    )
     for (i in 1:num_sim) {
      res_orig <- run_study(true_effect = true_effect_HA, n = n, sd = sd_val)
      res_rep <- run_study(true_effect = true_effect_HA, n = n, sd = sd_val)
      results_HA_pairs$p_original[i] <- res_orig$p.value
      results_HA_pairs$d_original[i] <- res_orig$d_observed
      results_HA_pairs$p_replication[i] <- res_rep$p.value
      results_HA_pairs$d_replication[i] <- res_rep$d_observed
    }
    significant_originals_HA <- results_HA_pairs[results_HA_pairs$p_original < alpha & !is.na(results_HA_pairs$p_original), ]
    num_significant_HA <- nrow(significant_originals_HA)


    # --- 整理 H_0 結果 ---
     results_H0_pairs <- data.frame(
        p_original = numeric(num_sim), d_original = numeric(num_sim),
        p_replication = numeric(num_sim), d_replication = numeric(num_sim)
    )
     for (i in 1:num_sim) {
      res_orig <- run_study(true_effect = true_effect_H0, n = n, sd = sd_val)
      res_rep <- run_study(true_effect = true_effect_H0, n = n, sd = sd_val)
      results_H0_pairs$p_original[i] <- res_orig$p.value
      results_H0_pairs$d_original[i] <- res_orig$d_observed
      results_H0_pairs$p_replication[i] <- res_rep$p.value
      results_H0_pairs$d_replication[i] <- res_rep$d_observed
    }
    significant_originals_H0 <- results_H0_pairs[results_H0_pairs$p_original < alpha & !is.na(results_H0_pairs$p_original), ]
    num_significant_H0 <- nrow(significant_originals_H0)

    # --- 格式化 H_A 文字結果 ---
    output_HA_text <- ""
    # 使用 results_HA_pairs 來計算顯著比例，因為它代表了獨立的原始研究
    prop_sig_HA <- mean(results_HA_pairs$p_original < alpha, na.rm = TRUE) * 100
    output_HA_text <- paste0(output_HA_text, sprintf("原始研究顯著的比例：%.2f%% (經驗功效)\n", prop_sig_HA))

    if (num_significant_HA > 0) {
      # 重複成功率
      rate_rep_HA <- mean(significant_originals_HA$p_replication < alpha, na.rm = TRUE) * 100
      output_HA_text <- paste0(output_HA_text, sprintf("整體重複成功率 (給定原始研究顯著)：%.2f%%\n", rate_rep_HA))

      # 贏家詛咒
      avg_d_all_HA <- mean(results_HA_pairs$d_original, na.rm = TRUE) # 使用 pairs 的原始 d
      avg_d_sig_HA <- mean(significant_originals_HA$d_original, na.rm = TRUE)
      output_HA_text <- paste0(output_HA_text, "\n--- 贏家詛咒分析 (H_A) ---\n")
      output_HA_text <- paste0(output_HA_text, sprintf("平均觀察效應量 (d) - 所有研究：%.3f\n", avg_d_all_HA))
      output_HA_text <- paste0(output_HA_text, sprintf("平均觀察效應量 (d) - *僅限顯著* 研究 ('贏家')：%.3f\n", avg_d_sig_HA))
      if(!is.na(avg_d_sig_HA) && !is.na(avg_d_all_HA) && avg_d_sig_HA > avg_d_all_HA) {
        output_HA_text <- paste0(output_HA_text, "  -> 顯著研究的平均效應量被高估。\n")
      }

      # 趨中迴歸
      avg_d_rep_for_winners <- mean(significant_originals_HA$d_replication, na.rm = TRUE)
      output_HA_text <- paste0(output_HA_text, "\n--- 趨中迴歸分析 (H_A) ---\n")
      output_HA_text <- paste0(output_HA_text, sprintf("重複研究的平均觀察效應量 (d) (針對原始'贏家')：%.3f\n", avg_d_rep_for_winners))
       if(!is.na(avg_d_rep_for_winners) && !is.na(avg_d_sig_HA) && avg_d_rep_for_winners < avg_d_sig_HA) {
        output_HA_text <- paste0(output_HA_text, "  -> 重複研究效應量趨向平均值。\n")
      }
    } else {
      output_HA_text <- paste0(output_HA_text, "\n在 H_A 情境下，無顯著的原始研究可供分析重複性/贏家詛咒。\n")
    }

    # --- 格式化 H_0 文字結果 ---
    output_H0_text <- ""
    prop_sig_H0 <- mean(results_H0_pairs$p_original < alpha, na.rm = TRUE) * 100 # 使用 pairs 的原始 p
    output_H0_text <- paste0(output_H0_text, sprintf("原始研究顯著的比例 (偽陽性)：%.2f%% (經驗第一類型錯誤率)\n", prop_sig_H0))
    output_H0_text <- paste0(output_H0_text, sprintf("  (應接近 alpha = %.1f%%)\n", alpha * 100))

    if (num_significant_H0 > 0) {
       # 偽陽性重複率
      rate_rep_H0 <- mean(significant_originals_H0$p_replication < alpha, na.rm = TRUE) * 100
      output_H0_text <- paste0(output_H0_text, sprintf("整體偽陽性重複率 (給定原始研究為偽陽性)：%.2f%%\n", rate_rep_H0))
      output_H0_text <- paste0(output_H0_text, sprintf("  (應接近 alpha = %.1f%%)\n", alpha * 100))

      # H0 下的贏家詛咒現象
      avg_d_all_H0 <- mean(results_H0_pairs$d_original, na.rm = TRUE) # 使用 pairs 的原始 d
      avg_abs_d_sig_H0 <- mean(abs(significant_originals_H0$d_original), na.rm = TRUE)
      output_H0_text <- paste0(output_H0_text, "\n--- H0 情境下的贏家詛咒類似現象分析 ---\n")
      output_H0_text <- paste0(output_H0_text, sprintf("平均觀察效應量 (d) - 所有研究：%.3f (應接近 0)\n", avg_d_all_H0))
      output_H0_text <- paste0(output_H0_text, sprintf("平均 *絕對* 效應量 |d| - *僅限偽陽性* 研究 ('贏家')：%.3f\n", avg_abs_d_sig_H0))
      if(!is.na(avg_abs_d_sig_H0) && avg_abs_d_sig_H0 > abs(avg_d_all_H0) + 0.01) { # 檢查是否顯著大於 0
         output_H0_text <- paste0(output_H0_text, "  -> 偽陽性結果因隨機誤差顯示非零效應量。\n")
      }

       # H0 下的趨中迴歸
      avg_abs_d_rep_for_winners <- mean(abs(significant_originals_H0$d_replication), na.rm = TRUE)
      output_H0_text <- paste0(output_H0_text, "\n--- H0 情境下的趨中迴歸分析 ---\n")
      output_H0_text <- paste0(output_H0_text, sprintf("重複研究的平均 *絕對* 效應量 |d| (針對原始偽陽性)：%.3f\n", avg_abs_d_rep_for_winners))
       if(!is.na(avg_abs_d_rep_for_winners) && !is.na(avg_abs_d_sig_H0) && avg_abs_d_rep_for_winners < avg_abs_d_sig_H0) {
         output_H0_text <- paste0(output_H0_text, "  -> 重複研究的絕對效應量趨向零。\n")
      }
    } else {
       output_H0_text <- paste0(output_H0_text, "\n在 H_0 情境下，無顯著的（偽陽性）原始研究。\n")
    }

    # 返回包含文字結果和 p 值數據的列表，用於繪圖
    list(
      text_HA = output_HA_text,
      text_H0 = output_H0_text,
      p_values_HA = results_HA_pairs$p_original, # 使用 pairs 的原始 p 值
      p_values_H0 = results_H0_pairs$p_original  # 使用 pairs 的原始 p 值
    )
  })

  # --- 輸出：渲染文字結果 ---
  output$results_HA <- renderText({
    res <- simulation_results() # 觸發模擬
    req(res) # 確保結果可用
    res$text_HA
  })

  output$results_H0 <- renderText({
    res <- simulation_results() # 觸發模擬
    req(res) # 確保結果可用
    res$text_H0
  })

  # --- 輸出：渲染 p 值分佈圖 ---
  output$p_value_distribution_plot <- renderPlot({
      res <- simulation_results() # 獲取模擬結果
      req(res) # 確保結果可用
      alpha_val <- input$alpha # 獲取當前的 alpha 值

      # 準備繪圖數據
      plot_data <- data.frame(
          p_value = c(res$p_values_HA, res$p_values_H0),
          Hypothesis = factor(rep(c("H_A 為真", "H_0 為真"), each = length(res$p_values_HA)))
      )

      # 使用 ggplot2 繪製密度圖
      ggplot(plot_data, aes(x = p_value, fill = Hypothesis)) +
          geom_density(alpha = 0.6) + # 繪製密度曲線，設定透明度
          geom_vline(xintercept = alpha_val, linetype = "dashed", color = "red", size = 1) + # 加入 alpha 垂直線
          annotate("text", x = alpha_val, y = Inf, label = paste("alpha =", alpha_val), hjust = -0.1, vjust = 1.5, color = "red", size = 4) + # 標註 alpha 線
          scale_fill_manual(values = c("H_A 為真" = "skyblue", "H_0 為真" = "salmon")) + # 設定顏色
          labs(
              title = "原始研究 P 值分佈",
              x = "P 值",
              y = "密度",
              fill = "假設情境" # 圖例標題
          ) +
          theme_minimal(base_family = "sans") + # 使用簡潔主題，確保中文字體可用 (sans 通常可以)
          theme(
              plot.title = element_text(hjust = 0.5, size=16), # 標題置中加大
              legend.position = "bottom" # 圖例放底部
          ) +
          coord_cartesian(xlim = c(0, 1)) # 確保 x 軸範圍是 0 到 1

  })

}

# --- 執行 Shiny 應用程式 ---
shinyApp(ui = ui, server = server)

