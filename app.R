# Load necessary library
library(shiny)
library(stats) # For power.t.test

# --- Helper Function: Run a single study simulation ---
# This function simulates one study (t-test) and returns p-value and observed Cohen's d
run_study <- function(true_effect, n, sd = 1) {
  # Simulate two groups
  group1 <- rnorm(n, mean = 0, sd = sd)
  group2 <- rnorm(n, mean = true_effect, sd = sd)

  # Perform t-test
  t_test_result <- t.test(group2, group1, var.equal = TRUE)

  # Calculate observed Cohen's d
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  var1 <- var(group1)
  var2 <- var(group2)
  # Pooled standard deviation
  sd_pooled <- sqrt(((n - 1) * var1 + (n - 1) * var2) / (n + n - 2))
  # Handle rare case of zero sd_pooled
  if (sd_pooled == 0 || is.na(sd_pooled)) {
    d_observed <- NA
  } else {
    d_observed <- (mean2 - mean1) / sd_pooled
  }

  # Return results
  return(list(p.value = t_test_result$p.value, d_observed = d_observed))
}


# --- Shiny UI Definition ---
ui <- fluidPage(
  titlePanel("Simulation: Replication, Power & Winner's Curse"),

  sidebarLayout(
    sidebarPanel(
      h4("Simulation Parameters"),
      sliderInput("n_per_group", "Sample Size per Group (n):",
                  min = 10, max = 200, value = 30, step = 5),
      sliderInput("true_effect_HA", "True Effect Size (Cohen's d) for H_A:",
                  min = 0, max = 1.5, value = 0.5, step = 0.05),
      sliderInput("alpha", "Significance Level (alpha):",
                  min = 0.001, max = 0.1, value = 0.05, step = 0.001),
      numericInput("num_simulations", "Number of Simulations:",
                   value = 2000, min = 100, max = 20000, step = 100), # Reduced default for faster interaction
      actionButton("run_sim", "Run Simulation", icon = icon("play"), class = "btn-success"),
      hr(),
      h4("Theoretical Power"),
      textOutput("theoretical_power_display"),
      hr(),
      h4("Location Info"),
      # Displaying location info as requested in context
      p(strong("Current Location:"), "屏東縣屏東市 (Pingtung City, Pingtung County, Taiwan)"),
      p(strong("Current Time:"), textOutput("current_time_display", inline = TRUE))

    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Simulation Results",
                 h4("Scenario 1: H_A is True (Real Effect Exists)"),
                 verbatimTextOutput("results_HA"),
                 hr(),
                 h4("Scenario 2: H_0 is True (No Real Effect)"),
                 verbatimTextOutput("results_H0")
        ),
        tabPanel("Explanations",
                 h4("About This Simulation"),
                 p("This application simulates pairs of 'original' and 'replication' studies under two scenarios:"),
                 tags$ul(
                   tags$li(strong("Scenario 1 (H_A True):"), "Assumes a real effect exists in the population with the specified 'True Effect Size'."),
                   tags$li(strong("Scenario 2 (H_0 True):"), "Assumes there is no real effect in the population (true effect size is zero).")
                 ),
                 p("Adjust the parameters in the sidebar and click 'Run Simulation' to see how they influence the results."),

                 h4("Key Concepts & Interpretation"),
                 strong("Statistical Power:"),
                 p("The probability that a study will detect a real effect (i.e., find p < alpha) when that effect actually exists (H_A is true). Higher power is better. It depends on effect size, sample size, alpha, and data variability. The sidebar shows the *theoretical* power calculated for the H_A scenario based on your inputs."),
                 p(em("In Scenario 1 (H_A True):"), "The simulation estimates power empirically ('Proportion of significant original studies'). The 'Overall Replication Rate' for significant studies should also be close to the power, assuming the replication study has the same design."),

                 strong("Replication Rate:"),
                 p("The probability that a study reporting a significant result (p < alpha) will again produce a significant result when repeated under identical conditions."),
                 p(em("In Scenario 1 (H_A True):"), "This rate is primarily driven by the study's statistical power."),
                 p(em("In Scenario 2 (H_0 True):"), "Here, an original 'significant' result is a false positive (Type I error). The chance of getting *another* false positive in the replication study is simply the alpha level."),

                 strong("Winner's Curse:"),
                 p("The phenomenon where studies selected for significance ('winners') tend to overestimate the true effect size. This happens because random chance contributing to the significance also tends to inflate the observed effect."),
                 p(em("In Scenario 1 (H_A True):"), "Compare the 'Average observed effect size (d)' for 'All H_A studies' vs. 'Only significant H_A studies'. The latter ('winners') typically shows a higher average d, illustrating the curse."),
                 p(em("In Scenario 2 (H_0 True):"), "A similar effect occurs. False positive studies ('winners') must, by chance, show an effect size different from the true value of zero. We look at the 'Average *absolute* observed effect size |d|', which will be greater than zero for these false positives."),

                 strong("Regression Toward the Mean:"),
                 p("When you repeat a study that initially yielded an extreme result (like a very large effect size associated with significance), the result in the replication study is likely to be closer to the true average effect size."),
                 p(em("In Scenario 1 (H_A True):"), "Compare the average d for the original 'winners' with the average d in their corresponding 'replication studies'. The replication average d is usually smaller, closer to the true effect size."),
                 p(em("In Scenario 2 (H_0 True):"), "Compare the average *absolute* d for the original false positives with the average *absolute* d in their replications. The replication average |d| is usually smaller, closer to the true value of zero.")
        )
      )
    )
  )
)

# --- Shiny Server Logic ---
server <- function(input, output, session) {

  # Reactive expression to calculate theoretical power
  theoretical_power <- reactive({
    # Use tryCatch for cases where power calculation might fail (e.g., n too small)
    tryCatch({
      power_result <- power.t.test(
        n = input$n_per_group,
        delta = input$true_effect_HA,
        sd = 1, # Assuming sd=1 as in simulation
        sig.level = input$alpha,
        type = "two.sample",
        alternative = "two.sided"
      )
      # Return formatted power
      sprintf("%.1f%%", power_result$power * 100)
    }, error = function(e) {
      "Calculation Error" # Return error message if calculation fails
    })
  })

  # Display theoretical power
  output$theoretical_power_display <- renderText({
    theoretical_power()
  })

  # Reactive expression to run simulation when button is clicked
  simulation_results <- eventReactive(input$run_sim, {
    # Show a notification that simulation is running
    showNotification("Running simulations...", type = "message", duration = 2)

    # Get inputs
    n <- input$n_per_group
    true_effect_HA <- input$true_effect_HA
    alpha <- input$alpha
    num_sim <- input$num_simulations
    true_effect_H0 <- 0
    sd_val <- 1 # Consistent with power calculation

    # --- Run HA Simulation ---
    results_HA <- data.frame(p_original = numeric(num_sim), d_original = numeric(num_sim),
                             p_replication = numeric(num_sim), d_replication = numeric(num_sim))
    for (i in 1:num_sim) {
      res_orig <- run_study(true_effect = true_effect_HA, n = n, sd = sd_val)
      res_rep <- run_study(true_effect = true_effect_HA, n = n, sd = sd_val)
      results_HA$p_original[i] <- res_orig$p.value
      results_HA$d_original[i] <- res_orig$d_observed
      results_HA$p_replication[i] <- res_rep$p.value
      results_HA$d_replication[i] <- res_rep$d_observed
    }
    significant_originals_HA <- results_HA[results_HA$p_original < alpha & !is.na(results_HA$p_original), ]
    num_significant_HA <- nrow(significant_originals_HA)

    # --- Run H0 Simulation ---
    results_H0 <- data.frame(p_original = numeric(num_sim), d_original = numeric(num_sim),
                             p_replication = numeric(num_sim), d_replication = numeric(num_sim))
    for (i in 1:num_sim) {
      res_orig <- run_study(true_effect = true_effect_H0, n = n, sd = sd_val)
      res_rep <- run_study(true_effect = true_effect_H0, n = n, sd = sd_val)
      results_H0$p_original[i] <- res_orig$p.value
      results_H0$d_original[i] <- res_orig$d_observed
      results_H0$p_replication[i] <- res_rep$p.value
      results_H0$d_replication[i] <- res_rep$d_observed
    }
    significant_originals_H0 <- results_H0[results_H0$p_original < alpha & !is.na(results_H0$p_original), ]
    num_significant_H0 <- nrow(significant_originals_H0)

    # --- Format HA Results ---
    output_HA <- ""
    prop_sig_HA <- (num_significant_HA / num_sim) * 100
    output_HA <- paste0(output_HA, sprintf("Proportion of significant original studies: %.2f%% (Empirical Power)\n", prop_sig_HA))

    if (num_significant_HA > 0) {
      # Replication Rate
      rate_rep_HA <- mean(significant_originals_HA$p_replication < alpha, na.rm = TRUE) * 100
      output_HA <- paste0(output_HA, sprintf("Overall Replication Rate (given original sig.): %.2f%%\n", rate_rep_HA))

      # Winner's Curse
      avg_d_all_HA <- mean(results_HA$d_original, na.rm = TRUE)
      avg_d_sig_HA <- mean(significant_originals_HA$d_original, na.rm = TRUE)
      output_HA <- paste0(output_HA, "\n--- Winner's Curse Analysis (H_A) ---\n")
      output_HA <- paste0(output_HA, sprintf("Avg. observed effect size (d) - All studies: %.3f\n", avg_d_all_HA))
      output_HA <- paste0(output_HA, sprintf("Avg. observed effect size (d) - *Significant* studies ('Winners'): %.3f\n", avg_d_sig_HA))
      if(!is.na(avg_d_sig_HA) && !is.na(avg_d_all_HA) && avg_d_sig_HA > avg_d_all_HA) {
        output_HA <- paste0(output_HA, "  -> Effect size overestimated in significant studies.\n")
      }

      # Regression to Mean
      avg_d_rep_for_winners <- mean(significant_originals_HA$d_replication, na.rm = TRUE)
      output_HA <- paste0(output_HA, "\n--- Regression to Mean Analysis (H_A) ---\n")
      output_HA <- paste0(output_HA, sprintf("Avg. observed effect size (d) in Replications (for original 'Winners'): %.3f\n", avg_d_rep_for_winners))
       if(!is.na(avg_d_rep_for_winners) && !is.na(avg_d_sig_HA) && avg_d_rep_for_winners < avg_d_sig_HA) {
        output_HA <- paste0(output_HA, "  -> Replication effect size regressed towards the mean.\n")
      }
    } else {
      output_HA <- paste0(output_HA, "\nNo significant original studies found in H_A scenario to analyze replication/winner's curse.\n")
    }

    # --- Format H0 Results ---
    output_H0 <- ""
    prop_sig_H0 <- (num_significant_H0 / num_sim) * 100
    output_H0 <- paste0(output_H0, sprintf("Proportion of significant original studies (False Positives): %.2f%% (Empirical Type I Error Rate)\n", prop_sig_H0))
    output_H0 <- paste0(output_H0, sprintf("  (Should be close to alpha = %.1f%%)\n", alpha * 100))


    if (num_significant_H0 > 0) {
       # Replication Rate (False Positives)
      rate_rep_H0 <- mean(significant_originals_H0$p_replication < alpha, na.rm = TRUE) * 100
      output_H0 <- paste0(output_H0, sprintf("Overall False Positive Replication Rate (given original FP): %.2f%%\n", rate_rep_H0))
      output_H0 <- paste0(output_H0, sprintf("  (Should be close to alpha = %.1f%%)\n", alpha * 100))

      # Winner's Curse (H0 version)
      avg_d_all_H0 <- mean(results_H0$d_original, na.rm = TRUE)
      avg_abs_d_sig_H0 <- mean(abs(significant_originals_H0$d_original), na.rm = TRUE)
      output_H0 <- paste0(output_H0, "\n--- Winner's Curse Analysis (H_0) ---\n")
      output_H0 <- paste0(output_H0, sprintf("Avg. observed effect size (d) - All studies: %.3f (Should be ~0)\n", avg_d_all_H0))
      output_H0 <- paste0(output_H0, sprintf("Avg. *absolute* effect size |d| - *False Positive* studies ('Winners'): %.3f\n", avg_abs_d_sig_H0))
      if(!is.na(avg_abs_d_sig_H0) && avg_abs_d_sig_H0 > abs(avg_d_all_H0) + 0.01) { # Check if clearly > 0
         output_H0 <- paste0(output_H0, "  -> False positives show non-zero effect size due to chance.\n")
      }

       # Regression to Mean (H0 version)
      avg_abs_d_rep_for_winners <- mean(abs(significant_originals_H0$d_replication), na.rm = TRUE)
      output_H0 <- paste0(output_H0, "\n--- Regression to Mean Analysis (H_0) ---\n")
      output_H0 <- paste0(output_H0, sprintf("Avg. *absolute* effect size |d| in Replications (for original FPs): %.3f\n", avg_abs_d_rep_for_winners))
       if(!is.na(avg_abs_d_rep_for_winners) && !is.na(avg_abs_d_sig_H0) && avg_abs_d_rep_for_winners < avg_abs_d_sig_H0) {
         output_H0 <- paste0(output_H0, "  -> Replication absolute effect size regressed towards zero.\n")
      }

    } else {
       output_H0 <- paste0(output_H0, "\nNo false positive original studies found in H_0 scenario.\n")
    }

    # Return list of results
    list(HA = output_HA, H0 = output_H0)
  })

  # Render results
  output$results_HA <- renderText({
    res <- simulation_results() # Trigger simulation
    req(res) # Ensure results are available
    res$HA
  })

  output$results_H0 <- renderText({
    res <- simulation_results() # Trigger simulation
    req(res) # Ensure results are available
    res$H0
  })

  # Display current time reactively
   output$current_time_display <- renderText({
    invalidateLater(1000, session) # Update every second
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  })

}

# Run the Shiny app
shinyApp(ui = ui, server = server)
