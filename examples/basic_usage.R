library(spatconnect)

x <- c(2.4, 1.7, 0.6, -0.2, -1.1, -1.8)
edges <- data.frame(i = 1:5, j = 2:6)

# Global superlevel and sublevel analyses.
g <- global_connectivity(
  x,
  edges,
  reference = "mean",
  scale = "sd",
  side = "both",
  n_perm = 999,
  seed = 123
)

print(g)
g$superlevel$curve
g$sublevel$curve

# Local activation--merging indicators and conditional inference.
l <- local_connectivity(
  x,
  edges,
  reference = "mean",
  scale = "sd",
  side = "both",
  n_perm = 999,
  seed = 123,
  progress = FALSE
)

l$superlevel$results
l$sublevel$results
