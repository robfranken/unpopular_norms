library(ggplot2)
library(dplyr)
#install.packages("ggpattern")
library(ggpattern)

# Function to plot conformist utility with filled regions, red reference line, and patterned overlap
plot_conformist_utility <- function(s = 10, e = 0, w = 30, z = 60, lambda1 = 0, lambda2 = 0) {
  proportions <- seq(0, 1, length.out = 100)  # proportion of trend-following neighbors p
  
  # decreasing marginal utility returns function
  f_p <- function(p, z, lambda) {
    if (lambda == 0) return(z * p) # linear case
    z * (1 - exp(-lambda * p)) / (1 - exp(-lambda))
  }
  
  # compute utilities
  utility_follow_trend <- e + f_p(proportions, z, lambda1)  
  utility_no_trend <- s + f_p(1 - proportions, w, lambda2)
  
  # new dataframe
  df <- data.frame(
    Proportion_Trend_Followers = proportions,
    Utility_Trend = utility_follow_trend,
    Utility_No_Trend = utility_no_trend
  )
  
  # define fill colors based on which utility is higher
  df <- df %>%
    mutate(Fill_Color = ifelse(Utility_Trend > Utility_No_Trend, "B=1", "B=0"))
  
  # find maximum utility for resisting the trend
  max_resist_utility <- max(df$Utility_No_Trend)
  max_resist_p <- df$Proportion_Trend_Followers[which.max(df$Utility_No_Trend)]
  
  # find region where the red line crosses the black curve
  df_pattern <- df %>%
    filter(Utility_Trend >= max_resist_utility)  # only where black line is above red line
  
  # Plot using ggplot2
  ggplot(df, aes(x = Proportion_Trend_Followers)) +
    geom_ribbon(aes(ymin = pmin(Utility_Trend, Utility_No_Trend), 
                    ymax = pmax(Utility_Trend, Utility_No_Trend), 
                    fill = Fill_Color), alpha = 0.3) +  # Fill between curves
    geom_ribbon(data = df_pattern, 
                aes(ymin = max_resist_utility, ymax = Utility_Trend), 
                fill = "black", pattern = "stripe", pattern_color = "red", pattern_density = 0.1, alpha = 0.5) +  # Striped overlap area
    geom_line(aes(y = Utility_Trend, color = "B=1"), size = 1.2) +
    geom_line(aes(y = Utility_No_Trend, color = "B=0"), size = 1.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") +  # Threshold φ = 0.5
    geom_segment(aes(x = max_resist_p, xend = 1, y = max_resist_utility, yend = max_resist_utility),
                 linetype = "dotted", color = "red", size = 1) +  # red reference line
    scale_x_continuous(breaks = seq(0, 1, 0.1), 
                       sec.axis = sec_axis(~ 1 - ., breaks = seq(0, 1, 0.1), name = "1-p")) +  # secondary x-axis
    scale_fill_manual(values = c("B=1" = "black", "B=0" = "gold"), guide = "none") +  # remove fill legend
    scale_color_manual(values = c("B=1" = "black", "B=0" = "gold")) +  # Line colors
    coord_cartesian(ylim = c(0, 25)) +  # Fix y-axis limits
    labs(
      x = "p",
      y = "Utility",
      col = NULL 
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(color = "black"),  # color for the 'p' axis label
      axis.title.x.top = element_text(color = "gold"),  # color for the '1-p' axis label
      legend.position = "bottom"
    )
}

# Call the function to generate the plot
#a <- plot_conformist_utility(s=10,w=30,z=60,lambda1=0,lambda2=0) + theme(legend.position = "none")
b <- plot_conformist_utility(s=10, e=5, w=12, z=14,lambda1=5.5,lambda2=0.5)
c <- plot_conformist_utility(s=15,w=40, z=50,lambda1=3,lambda2=1.8)

cowplot::plot_grid( b, c, labels = c("(A)", "(B)"), ncol = 2, label_size = 14)
