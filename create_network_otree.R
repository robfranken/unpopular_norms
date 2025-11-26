library(igraph)
library(jsonlite)  # for saving as JSON

# number of players
#n <- 6
n = 4

# generate skewed degree distribution (undirected, no self-loops)
set.seed(42)
g <- sample_degseq(c(2,2,2,2,4,4), method = "vl")

degseq <- sample(2:12, size = n, replace = TRUE, prob = (1/(2:12))^3)

if (sum(degseq)%%2 != 0) {
  degseq[1] <- degseq[1] + 1
}
g <- sample_degseq(degseq, method = "vl")

# convert to adjacency matrix
adj_matrix <- as.matrix(as_adjacency_matrix(g))

# proportion of minority players
p_t <- 0.15
n_m <- max(floor(n * p_t),1)  # number of minority players

# create role vector (0 = majority, 1 = minority)
role_vector <- rep(0, n)
#minority_indices <- sample(1:n, n_m)
#role_vector[minority_indices] <- 1

# make central agents the minority
role_vector[1] <- 1
role_vector[4] <- 1

# create a list to store the network data
net <- list(adj_matrix = adj_matrix, role_vector = role_vector)

# save the list as a JSON file
write_json(net, "network_test_n20.json")
