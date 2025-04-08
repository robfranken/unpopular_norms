#custom functions
#by Rob Franken
#last edited: 26-02-2025

# function to install/load packages
fpackage.check <- function(packages) {
  lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}

# function to save data with time stamp in correct directory
fsave <- function(x, file, location = "./data/processed/", ...) {
  if (!dir.exists(location))
    dir.create(location)
  datename <- substr(gsub("[:-]", "", Sys.time()), 1, 8)
  totalname <- paste(location, datename, file, sep = "")
  print(paste("SAVED: ", totalname, sep = ""))
  save(x, file = totalname)
}

# function to load R-objects under new names
fload <- function(fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# print objects (tibble / data.frame) nicely on screen in .Rmd.
fshowdf <- function(x, ...) {
  knitr::kable(x, digits = 2, "html", ...) %>%
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover")) %>%
    kableExtra::scroll_box(width = "100%", height = "300px")
}

# visualize graph in the "multi-population"
fplot_graph <- function(graph, main=NULL, layout_algo=NULL, 
                        col1 = "#800080", col2 =  "#FFD700", legend=TRUE) {
  plot(graph,
       main = main,
       layout = layout_algo,
       vertex.label = NA,
       vertex.size = degree(graph)*.5 + 3, # node size based on degree
       vertex.color = ifelse(V(graph)$role == "trendsetter", col1, col2),
       edge.width = 0.5,
       
       edge.color = "darkgrey")
  #add legend
  if(legend==TRUE) {
    legend("bottomleft",
           legend = c("Trendsetter", "Conformist"),
           pch = 21,
           col = c(col1,col2),
           pt.bg = c(col1,col2),
           pt.cex = 3,
           bty = "n")
  }
}

# utility function 
futility <- function(agent_id, choice, agents, network, params) {
  # get ego and his local neighborhood
  ego <- agents[agent_id, ]
  neighbors <- neighbors(network, ego$id)
  alters <- agents[as.numeric(neighbors), ]
  n <- nrow(alters)
  
  # proportion of neighbors who adopted (p) or resisted (1-p) the trend
  p <- sum(alters$choice == 1) / n
  q <- 1 - p  # proportion of trend-resistant neighbors
  
  # calculate diminishing returns function
  diminishing_returns <- function(x, v, lambda) {
    v * (1 - exp(-lambda * x)) / (1 - exp(-lambda))
  }
  
  # calculate expected utility (depending on alters' choices in prior round)
  if(ego$role == "conformist") {
    if (choice == 0) { # resist the trend
      choice_payoff <- params$s # fixed utility for resisting
      coordination_payoff <- diminishing_returns(q, params$w, params$lambda2)
    } else { #follow the trend
      choice_payoff <- 0
      coordination_payoff <- diminishing_returns(p, params$z, params$lambda1)
    }
  } else { # trendsetters only care about following
    if (choice == 1) {
      choice_payoff <- params$e
    } else {
      choice_payoff <- 0
    }
    coordination_payoff <- 0
  }
  
  # return total utility and counts of neighbors
  return(list(utility = choice_payoff + coordination_payoff,
              n1 = sum(alters$choice == 1),  # count of trend-adopting neighbors
              n0 = sum(alters$choice == 0))) # count of trend-resistant neighbors
}

# create a degree sequence based on a distribution
fdegseq <- function(n, dist = "power-law", alpha, k_min = 1, k_max = n - 1, seed = NULL) {
  
  if (!is.null(seed)) { 
    set.seed(seed)
  }
  
  # generate a degree sequence based on a power-law distribution of the form p(k) ∝ k^{-alpha}
  
  # create probability distribution
  probs <- (1 / (k_min:k_max))^alpha
  
  # normalize probabilities
  probs <- probs / sum(probs) 
  
  # sample a degree sequence
  degseq <- sample(k_min:k_max, size = n, replace = TRUE, prob = (1/(k_min:k_max))^alpha)
  
  # correct the degree sequence if its sum is odd (necessary for the configuration model)
  if (sum(degseq) %% 2 != 0) {
    degseq[1] <- degseq[1] + 1
  }
  
  if (dist == "power-law") {
    # if the specified distribution type is power-law, return the degree sequence
    return(degseq)
    
  } else if (dist == "log-normal") {
    # if the specified distribution type is log-normal, generate a degree sequence following this distribution;
    # but with a same mean degree <k> as its 'scale-free' counterpart. the spread alpha = 1. 
    
    # calculate mean degree in sequence
    mean_deg <- mean(degseq)
    
    # generate raw log-normal degree sequence
    raw_degseq <- rlnorm(n = n, meanlog = log(mean_deg), sdlog = 1)
    
    # re-scale the degree sequence to match the target mean degree
    scaling_fctr <- mean_deg / mean(raw_degseq)
    scaled_degseq <- raw_degseq * scaling_fctr
    
    # apply bounds [k_min, k_max] and round to integer values
    bounded_degseq <- pmin(pmax(round(scaled_degseq), k_min), k_max)
    
    # since the mean may shift after rounding and bounding, re-scale again:
    scaling_fctr2 <- mean_deg / mean(bounded_degseq)
    degseq <- pmin(pmax(round(bounded_degseq * scaling_fctr2), k_min), k_max)
    
    # correct the degree sequence if its sum is odd (necessary for the configuration model)
    if (sum(degseq) %% 2 != 0) {
      degseq[1] <- degseq[1] + 1
    }
    return(degseq)
  } 
  else {
    stop("Invalid distribution type. Please choose 'power-law' or 'log-normal'.")
  }
}

# create a "random" network with a minimum degree restriction 
frandkmin <- function(n, p, k_min = 3, max_iter = 1e3, verbose = TRUE, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)  
  }
  
  # ER network
  network <- erdos.renyi.game(n, p)
  
  iter <- 0
  while (iter < max_iter) {
    deg <- degree(network)
    
    # look for nodes whose degree is less than k_min
    nodes_below_kmin <- which(deg < k_min)
    
    if (length(nodes_below_kmin) == 0) {
      if (verbose) print(paste("Iteration", iter + 1, ": All nodes have degree >= k_min. Stopping..."))
      
      break
    }
    
    # process nodes that dont meet k_min
    for (i in nodes_below_kmin) {
      if (deg[i] >= k_min) next 
      if (verbose) print(paste("Node", i, "has degree", deg[i], "< k_min. Stealing degree..."))
      
      # find a node j with degree > k_min
      neighbors_i <- neighbors(network, i)  # get neighbors of node i
      candidates_j <- setdiff(which(deg > k_min), nodes_below_kmin)  # nodes with degree > k_min and not already below k_min
      
      if (length(candidates_j) == 0) {
        if (verbose) print(paste("No candidate j found for node", i, ". Skipping."))
        next
      }
      
      j <- sample(candidates_j, 1)  # sample one node j from candidates
      
      if (verbose) print(paste("Node", j, "is a candidate for stealing degree."))
      
      # find a neighbor m of node j to "steal" the edge from
      neighbors_j <- neighbors(network, j)
      m <- sample(neighbors_j, 1)  # Randomly select one of j's neighbors
      
      if (deg[m] <= k_min) {
        if (verbose) print(paste("Node", m, "has degree <= k_min. Skipping this neighbor."))
        next
      }
      
      if (verbose) print(paste("Node", j, "is connected to node", m, ". Stealing edge."))
      
      # check if the edge exists between j and m before attempting to remove it
      edge_exists <- any(ends(network, E(network))[,1] == j & ends(network, E(network))[,2] == m |
                           ends(network, E(network))[,1] == m & ends(network, E(network))[,2] == j)
      
      if (edge_exists) {
        # remove the edge between j and m, and add an edge between i and m
        network <- delete_edges(network, get.edge.ids(network, c(j, m)))  # use edge ID instead
        network <- add_edges(network, c(i, m))    # ddd the edge between i and m
        
        # update the degree vector
        deg <- degree(network)
        
        if (verbose) print(paste("Node", i, "now has degree", deg[i]))
      } else {
        if (verbose) print(paste("No edge exists between node", j, "and node", m, "Skipping..."))
      }
    }
    
    # increment iteration counter
    iter <- iter + 1
    if (verbose) print(paste("Iteration", iter, "completed."))
  }
  
  return(network)
}

# rewiring algorithm to manipulate a network's degree assortativity (without changing the degree distribution)
frewire_r <- function(network, target_r, max_iter = 1e5, tol = 0.01, verbose = TRUE) {
  
  current_r <- assortativity_degree(network)
  iteration <- 1
  
  if (verbose) {
    cat("Target assortativity coefficient:", target_r, "\n")
    cat("Starting assortativity coefficient:", current_r, "\n")
    cat("Tolerance:", tol, "\n")
  }
  
  while (abs(current_r - target_r) > tol && iteration < max_iter) {
    
    # get network edges
    edges <- E(network)
    # to edgelist
    edge_list <- ends(network, edges)
    
    # randomly select two pairs of connected nodes
    idx1 <- sample(1:nrow(edge_list), 1)
    idx2 <- sample(1:nrow(edge_list), 1)
    
    # extract node indices
    u1 <- edge_list[idx1, 1] # node 1 of first edge
    v1 <- edge_list[idx1, 2] # node 2 of first edge
    u2 <- edge_list[idx2, 1] # etc
    v2 <- edge_list[idx2, 2] 
    
    # check if the two pairs of connected nodes (u1, v1; u2, v2) are disjoint
    if (length(unique(c(u1, v1, u2, v2))) == 4) {
      # check if there is already an edge across the node-pairs
      # ensure no loops and no duplicate edges
      if (!are_adjacent(network, u1, u2) && !are_adjacent(network, v1, v2) && u1 != v2 && u2 != v1) {
        
        # perform the edge swap (u1,v1) <-> (u2,v2) becomes (u1,v2) <-> (u2,v1)
        new_network <- network # copy network
        
        # check if the new edges already exist to avoid duplicates
        if (!are_adjacent(new_network, u1, v2) && !are_adjacent(new_network, u2, v1)) {
          # add edges
          new_network <- add_edges(new_network, c(u1, v2, u2, v1))
          # remove edges
          new_network <- delete_edges(new_network, get.edge.ids(new_network, c(u1, v1, u2, v2)))
          
          # new assortativity
          new_r <- assortativity_degree(new_network)
          
          # accept tie swap if it brings us closer to the target assortativity
          if (abs(new_r - target_r) < abs(current_r - target_r)) {
            current_r <- new_r
            network <- new_network
            if (verbose) { 
              cat("Rewiring at iteration", iteration, "brought assortativity closer to target! Current assortativity coefficient:", new_r, "\n")
            }
          }
        }
      }
    }
    iteration <- iteration + 1
  }
  
  if (verbose) {
    cat("Final assortativity coefficient:", current_r, "\n")
    if (abs(current_r - target_r) <= tol) {
      cat("Target reached within tolerance.\n")
    } else {
      cat("Reached maximum iterations without meeting target.\n")
    }
  }
  
  return(network)
}

# algorithm to manipulate the degree-trait correlation
fdegtraitcor <- function(network) {
  roles <- ifelse(V(network)$role == "trendsetter", 1, 0)
  degrees <- degree(network)
  return(list(cor = cor(roles, degrees), roles = roles, degrees = degrees))
}

#swapping function to adjust degree-trait correlation
fswap_rho <- function(network, target_rho, max_iter = 1e3, tol = 0.05, verbose = TRUE) {
  
  current <- fdegtraitcor(network)
  iteration <- 1
  best_network <- network
  best_rho <- current$cor
  
  if (verbose) {
    cat("Target degree-trait correlation:", target_rho, "\n")
    cat("Starting degree-trait correlation:", current$cor, "\n")
    cat("Tolerance:", tol, "\n\n")
  }
  
  while (iteration <= max_iter) {
    # check if we are already within tolerance
    if (abs(current$cor - target_rho) <= tol) {
      if (verbose) cat("Target reached within tolerance at iteration", iteration, ".\n")
      break
    }
    
    # randomly select nodes for swapping
    v1 <- sample(which(current$roles == 1), 1)
    v0 <- sample(which(current$roles == 0), 1)
    
    # get degrees of selected nodes
    k1 <- current$degrees[v1]
    k0 <- current$degrees[v0]
    
    # swap roles if condition k_v0 > k_v1 is met
    if (k0 > k1) {
      current$roles[v1] <- 0
      current$roles[v0] <- 1
      
      # update graph roles
      V(network)$role <- ifelse(current$roles == 1, "trendsetter", "conformist")
      
      # recalculate degree-trait correlation
      current <- fdegtraitcor(network)
      
      # check if this is the closest correlation to the target so far
      if (abs(current$cor - target_rho) < abs(best_rho - target_rho)) {
        best_network <- network
        best_rho <- current$cor
        if (verbose) {
          cat("Trait-swapping at iteration", iteration, "brought correlation closer to target! Current correlation:", current$cor, "\n")
        }
      }
    }
    iteration <- iteration + 1
  }
  
  # check if the final correlation is worse than the best correlation
  final_correlation <- current$cor
  if (abs(final_correlation - target_rho) > abs(best_rho - target_rho)) {
    if (verbose) {
      cat("\nWarning: Final iteration made the correlation worse. Reverting to best observed correlation.\n")
    }
  }
  
  if (verbose) {
    cat("\nFinal degree-trait correlation:", best_rho, "\n")
    if (abs(best_rho - target_rho) <= tol) {
      cat("Target reached within tolerance.\n")
    } else if (iteration > max_iter) {
      cat("Reached maximum iterations without meeting target.\n")
    }
  }
  return(best_network)
}

fcalculate_majority_illusion <- function(network, threshold = 0.49) {
  roles <- V(network)$role
  actions <- V(network)$action
  
  #initialize counter for majority illusion
  mi_count <- 0
  
  #loop over conformists
  for (v in V(network)) {
    if (roles[v] == "conformist") {
      neighbors <- neighbors(network, v)
      trend_neighbors <- sum(actions[neighbors] == 1)
      prop_trend <- trend_neighbors / length(neighbors)
      
      if (prop_trend > threshold) { #under "weak influence", more than half of the neighborhood suffices; so set threshold at .5
        mi_count <- mi_count + 1
      }
    }
  }
  # return fraction of conformists who have majority illusion
  return(mi_count / sum(roles == "conformist"))
}

# function for our evolutionary model
fabm <- function(network = network, # the generated network
                 rounds = 10, # number of timesteps/rounds
                 choice_rule = "deterministic", # choice update rule
                 utility_fn = futility, # the utility function
                 params = list(s=15, e=10, w=40, z=50, lambda1=4.3, lamda2=1.8), # utility parameters
                 mi_threshold = ifelse(params$z > 50, .49, .50), # under the "strong influence" condition (i.e., 50<z<90) half of neighbors adopting the trend is sufficient; under the "weak influence" condition (i.e., 30<z<50) *more* than half of neighbors adopting the trend is needed to be influenced.)
                 histories = FALSE, # return decision history
                 outcome = TRUE, # return outcomes
                 plot = FALSE ) { # return plot
  
  # make an agents dataframe
  agents <- tibble(
    id = 1:length(network),
    role = V(network)$role,
    preference = ifelse(role == "trendsetter", 1, 0), # 1 = follow trend, 0 = not follow
    choice = NA,
    beta = rgamma(length(network), shape = 2, scale = 1.5)) # Gamma-distributed noise level (for conformists' decisions under probabilistic choice rule)
  
  # initialize decision history
  decision_history <- tibble()
  
  # also initialize an equilibrium flag
  equilibrium_reached <- FALSE
  equilibrium_t <- NA
  
  # simulation loop
  for (t in 1:rounds) {
    if (t == 1) {
      # round 1: agents make decisions based on private preferences (no social information available yet)
      agents <- agents %>%
        mutate(choice = preference)
    } else {
      # round t > 1: agents make decisions based on social information from neighbors (decisions at t-1)
      agents <- agents %>%
        rowwise() %>%
        mutate(
          util_1 = utility_fn(id, 1, agents, network, list(s = params$s, e = params$e, w = params$w, z = params$z, lambda1 = params$lambda1, lambda2 = params$lambda2))$utility,
          util_0 = utility_fn(id, 0, agents, network, list(s = params$s, e = params$e, w = params$w, z = params$z, lambda1 = params$lambda1, lambda2 = params$lambda2))$utility,
          
          #compute probability of following the trend for conformists
          prob_follow = 1 / ( 1 + exp(-beta * (util_1 - util_0))),
          
          # make decision based on expected utility
          choice = ifelse(
            role == "trendsetter",
            ifelse(util_1 > util_0, 1, 0), #always deterministic
            ifelse(
              choice_rule == "deterministic",
              ifelse(util_1 > util_0, 1, 0), # deterministic if specified
              rbinom(1, 1, prob_follow) # probabilistic otherwise
            )
          )
        )
    }
    # store the decisions for this round
    decision_history <- bind_rows(decision_history, agents %>% 
                                    mutate(round = t))
  }
  
  # check for equilibrium
  # @RF: when introducing some stochasticity, this should be handled differently...
  for(t in 2:rounds ) { 
    prev_trend_followers <- mean(decision_history %>% filter(round == t - 1) %>% pull(choice) == 1)
    curr_trend_followers <- mean(decision_history %>% filter(round == t) %>% pull(choice) == 1)
    
    if (!equilibrium_reached && abs(curr_trend_followers - prev_trend_followers) < .Machine$double.eps) {
      equilibrium_reached <- TRUE
      equilibrium_t <- t - 1  
    }
  }
  
  # based on decision_history:
  
  # 1. calculate global majority illusion over rounds
  globMI <- numeric(rounds)
  for (t in 1:rounds) {
    if (t == 1) {
      # in round 1: no social information, so no MI
      globMI[t] <- NA
    } else {
      # rounds t > 1: calculate magnitude of majority illusion
      # first, make a copy of the network object
      exposure_network <- network
      # update the actions of actors, based on their choices in the previous round
      # after all, actors don't observe others' roles, but only their choices,
      # based on which they can infer their role
      V(exposure_network)$action <- decision_history[decision_history$round == t - 1,]$choice
      globMI[t] <- fcalculate_majority_illusion(exposure_network, threshold = mi_threshold)  
    }
  }
  
  # to long format
  MI <- tibble(
    round = 1:rounds,
    outcome = "majority_illusion",
    statistic = globMI
  )
  
  # 2. calculate the evolution of the unpopular norm
  UN <- decision_history %>%
    group_by(round) %>%
    summarise(
      follow_trend = mean(choice == 1, na.rm = TRUE)
    ) %>%
    pivot_longer(cols = c("follow_trend"),
                 names_to = "outcome",
                 values_to = "statistic")
  
  # bind
  plotdata <- bind_rows(MI, UN)
  
  if (plot) {
    fig <- ggplot(plotdata, aes(x=round, y=statistic, color=factor(outcome))) +
      geom_line() +
      geom_point() +
      scale_y_continuous(labels = scales::percent_format(scale = 100), limits = c(0,1)) +
      scale_x_continuous(breaks = seq(1, max(plotdata$round), by = 1)) +
      labs(
        title = "Evolution of an unpopular norm",
        subtitle = "`follow_trend` denotes the percentage of all agents that follow the trend.\n`majority_illusion` reflects the percentage of conformists whose neighbors meet or exceed the adoption threshold φ (i.e., 0.50),\nwith 'strong influence' referring to both meeting and exceeding the threshold, and 'weak influence' to exceeding it only.\nThe grey dashed line reflects the percentage of agents whose role is trendsetter. The purple circle depicts the equilibrium point.",
        x = "round",
        y = "% agents",
        color = "outcome") +
      theme(panel.grid.minor.x = element_blank(),
            #legend.position = "bottom",
            plot.subtitle = element_text(size = 8)) +
      geom_hline(yintercept = prop.table(table(agents$role))[2], linetype = "dashed", color = "darkgrey", size = 1)
    
    # add a circle around the equilibrium state 'follow_trend' statistic
    if (!is.na(equilibrium_t)) {
      fig <- fig + geom_point(
        data = plotdata %>% filter(round == equilibrium_t & outcome == "follow_trend"),
        aes(x = round, y = statistic), 
        shape = 1, size = 4, color = "purple", stroke = 2
      )
    }
  }
  
  # return outputs
  output <- list()
  if (histories) { 
    output$decision_history <- decision_history}
  
  if (outcome) { 
    output$outcomes <- plotdata
    output$equilibrium <- list(
      reached = equilibrium_reached,
      round = equilibrium_t,
      prop_follow_trend = if (!is.na(equilibrium_t)) {
        mean(decision_history %>% filter(round == equilibrium_t) %>% pull(choice) == 1)
      } else { NA }
    )
  }
  if (plot) { output$plot <- fig}
  return(output)
}

# function to create a plot displaying the outputs (% norm followers) across simulated networks
fcreate__outcome_plot <- function(p_t_level) {
  data_filtered <- data %>% filter(p_t == p_t_level)
  
  ggplot(data_filtered, aes(x = actual_rho, y = actual_r, color = prop_trend)) +
    geom_point(size = 1.5, alpha = 0.7) +
    facet_grid2(
      rows = vars(dist),
      cols = vars(alpha),
      labeller = labeller(alpha = function(x) paste0("\u03B1=", x))
    ) +
    scale_color_gradientn(
      colors = c("lightblue", "yellow", "red", "black"),
      values = scales::rescale(c(0, 0.2, 0.5, 0.75, 0.9, 1)),
      limits = c(0, 1),
      name = "% agents\nfollowing"~italic(B)~""
    ) +
    labs(
      x = expression(rho[kx]),
      y = expression(r[kk]),
      title = paste("proportion minority =", p_t_level)
    ) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "grey90", color = "grey50")
    )
}

# function to display aggregated outcomes (probability of a negative cascade) across the parameter space
fcreate_heatmap <- function(p_t_level) {
  # filter data for the given level of p_t
  data_filtered <- data %>% filter(p_t == p_t_level)
  
  # summarize the data
  data_summary <- data_filtered %>%
    group_by(target_rho, target_r, alpha, dist) %>%
    summarize(
      `P(neg. cascade)` = mean(prop_trend == 1, na.rm = TRUE),  # mean of prop_trend == 1
      n = n(),  # count of observations in each group
      se_prob = sqrt((`P(neg. cascade)` * (1 - `P(neg. cascade)` )) / n),  # standard error for the proportion
      .groups = "drop"
    )
  
  # find the highest value for each facet
  max_values <- data_summary %>%
    group_by(alpha, dist) %>%
    filter(`P(neg. cascade)` == max(`P(neg. cascade)`)) %>%
    distinct(alpha, dist, .keep_all = TRUE) %>%
    ungroup()
  
  # create the heatmap
  ggplot(data_summary, aes(x = target_rho, y = target_r, fill = `P(neg. cascade)`)) +
    geom_tile(color = "white") +  #  with gridlines
    scale_fill_gradientn(
      colors = c("lightblue", "yellow", "red", "black"), 
      name = "P(neg. cascade)",
      limits = c(0, 1)  # ensure the fill scale ranges between 0 and 1
    ) +
    facet_grid2(
      rows = vars(dist),
      cols = vars(alpha),
      scales = "free",
      independent = TRUE,
      labeller = labeller(alpha = function(x) paste0("\u03B1=", x)) 
    ) +
    # add a line connecting the text box to the corresponding tile (drawn first)
    #geom_segment(data = max_values, 
    #             aes(x = target_rho, y = target_r, 
    #                 xend = target_rho - 0.3, yend = target_r),
    #             color = "black", size = 0.5) +  # thin line connecting text and tile
    # add a point at the center of the cell
    #geom_point(data = max_values, 
    #           aes(x = target_rho, y = target_r),
    #           color = "black", size = 1) +  # point (dot) at the center of the cell
    # add a label box with the highest value per panel (light fill color)
    #geom_label(data = max_values, 
  #           aes(label = paste0("P = ", round(`P(neg. cascade)`, 2))),
  #           fill = "gray90", color = "black", size = 4, 
  #           label.padding = unit(0.2, "lines"), label.size = 0.2, 
  #           nudge_x = -0.3, nudge_y = 0) +  # light gray background
  labs(
    x = expression(rho[kx]), 
    y = expression(r[kk]), 
    title = paste("Proportion trendsetters =", p_t_level)  # add a title for each level of p_t
  ) +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      legend.position = "right"  # move legend to the bottom for consistent layout
    )
}

# function to generate a gif displaying the propagation of the norm through a (static) network
fnetworkgif <- function(network, decision_history, rounds, fps = 2, output_dir = "./figures") {
  
  # colors for choices
  choice_color <- c("0" = brewer.pal(3, "Set3")[1],  
                    "1" = brewer.pal(3, "Set3")[2]) 
  
  # function to get the choice data for a specific round
  fround_data <- function(round_number, decision_history) {
    round_data <- decision_history %>%
      filter(round == round_number) %>%
      select(id, role, choice, round)
    return(round_data)
  }
  
  # function to update the network and set border color based on initial trendsetters
  fupdate_choices <- function(network, decision_history, round_number) {
    # get the decisions for the specific round
    round_data <- fround_data(round_number, decision_history)
    
    # set the node attribute based on actors' decisions
    V(network)$choice <- as.character(round_data$choice)
    
    # border properties to distinguish trendsetters
    V(network)$border_color <- ifelse(V(network)$role == "trendsetter", "black", "grey70")
    V(network)$border_size <- ifelse(V(network)$role == "trendsetter", 1, 0.5)
    return(network)
  }
  
  # function to generate a list of updated networks
  fupdate_nets <- function(network, decision_history, rounds) {
    snapshots <- list()
    
    for (round_number in 1:rounds) {
      updated_network <- fupdate_choices(network, decision_history, round_number)  
      snapshots[[round_number]] <- updated_network
    }
    return(snapshots)
  }
  
  # generate updated networks for all rounds
  evo <- fupdate_nets(network = network, decision_history = decision_history, rounds = rounds)
  
  # create a layout for each network
  xy <- layout_nicely(evo[[1]])
  
  # create a (static) plot list
  pList <- vector("list", length(evo))
  
  for (i in 1:length(evo)) {
    # calculate the percentage of b=1
    choice1_percentage <- sum(V(evo[[i]])$choice == "1") / length(V(evo[[i]])) * 100
    
    # update title to include percentage
    round_title <- paste("Round", i, " - % agents following <i>B</i> :", round(choice1_percentage, 2), "%")
    
    pList[[i]] <- ggraph(evo[[i]], layout = "manual", x = xy[, 1], y = xy[, 2]) +
      geom_edge_link0(edge_width = 0.3, edge_colour = "grey66") +
      geom_node_point(shape = 21, aes(fill = as.factor(choice), color = border_color), size = 4, stroke = V(evo[[i]])$border_size) +  
      scale_fill_manual(values = choice_color) +
      scale_color_identity() +
      theme_graph() +
      theme(legend.position = "none") +
      labs(title = round_title) +
      theme(plot.title = element_markdown())
  }
  
  # ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # save each plot to individual files in the specified directory
  for (i in 1:length(pList)) {
    ggsave(paste0(output_dir, "/plot_", i, ".png"), plot = pList[[i]], width = 6, height = 6)
  }
  
  # read saved images into a magick image object
  image_files <- list.files(path = output_dir, pattern = "plot_.*\\.png", full.names = TRUE)
  
  # sort the images by their numerical order (1 to X)
  image_files <- image_files[order(as.numeric(gsub(".*_(\\d+)\\.png", "\\1", image_files)))]
  
  img_list <- image_read(image_files)
  
  # create an animated gif
  gif <- image_animate(img_list, fps = fps)
  
  # save the gif to a file
  gif_output_path <- paste0(output_dir, "/animation.gif")
  image_write(gif, gif_output_path)
  
  # remove the temporary plot images
  file.remove(image_files)
  
  return(gif_output_path)
}