# see https://www.mjdenny.com/Preparing_Network_Data_In_R.html

library(visNetwork)
library(dplyr)
library(CVN)

load_kikme <- T

if (load_kikme) { 
  raw_input <- readRDS("../KiKme-analysis/data/KiKme-dataset.rds")
  data <- readRDS("../KiKme-analysis/results/KiKme-grid.rds") 
}

A <- as.matrix(data$adj_matrices[[6]][[1]])

nodes <- create_nodes_visnetwork(nrow(A), labels = raw_input$gene_labels)
edges <- create_edges_visnetwork(A)

subset_edges <- list(from = c(1, 114), to = c(42, 128))
edges <- set_attributes_to_edges_visnetwork(edges, subset_edges, width = c(8, 2), color = c("red", "blue"))

visnetwork(nodes, edges)

#nodes <- data.frame(id = 1:data$p)
#nodes$title <- raw_input$gene_labels # Text on click


create_list_nodes <- function(n_nodes, labels = 1:n_nodes) {
  nodes <- data.frame(id = 1:data$p)
  nodes$title <- labels
  return(nodes)
}


create_list_edges <- function(A) { 
  A[lower.tri(A)] <- 0
  edges <- as.data.frame(which(A != 0, arr.ind = T)) 
  edges %>% 
    rename(
      from = row,
      to = col
    ) %>% 
    arrange(from)
}

assign_attributes_to_edges <- function(edges, subset_edges, width = c(NA, NA), color = c(NULL, NULL)) { 
  edges$id <- 1:nrow(edges)
  
  ids <- edges %>% filter(
    from %in% subset_edges$from, 
    to %in% subset_edges$to
  )
  
  in_group <- ids$id
  out_group <- edges$id[-in_group]
  
  if (!is.na(width[1])) { 
    edges$width[1:nrow(edges)] <- NA
    edges$width[in_group] <- width[1]
    edges$width[out_group] <- width[2]
  }
  
  if (!is.null(color[1])) { 
    edges$color[1:nrow(edges)] <- NA
    edges$color[in_group] <- color[1]
    edges$color[out_group] <- color[2]
  }
  edges %>% select(-id)
}

visnetwork <- function(nodes, 
                       edges, 
                       node_titles = 1:nrow(adj_matrix), 
                       title = "", 
                       igraph_layout = "layout_in_circle") { 
  visNetwork(nodes, edges, width = "100%", main = list(text = title)) %>% 
    visIgraphLayout(layout = igraph_layout) %>%
    visOptions(highlightNearest = list(enabled = T, hover = T))
}


nodes <- create_list_nodes(nrow(A), labels = raw_input$gene_labels)
edges <- create_list_edges(A)

subset_edges <- list(from = c(1, 114), to = c(42, 128))
edges <- assign_attributes_to_edges(edges, subset_edges, width = c(8, 2), color = c("red", "blue"))

visnetwork(nodes, edges)
# 
# 
# 
# #edges$width <- c(rep(5, 5), rep(1, nrow(edges) - 5)) # line width
# #edges$color <- c(rep("red", 5), rep("blue", nrow(edges) - 5))    # line color  
# #links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
# #links$smooth <- FALSE    # should the edges be curved?
# #links$shadow <- FALSE    # edge shadow
# 
# subset_edges <- list(from = c(1, 114), to = c(42, 128))
# 
# assign_attributes_to_edges <- function(edges, subset_edges, width = c(NA, NA), color = c(NULL, NULL)) { 
#   edges$id <- 1:nrow(edges)
#   
#   ids <- edges %>% filter(
#     from %in% subset_edges$from, 
#     to %in% subset_edges$to
#   )
#   
#   in_group <- ids$id
#   out_group <- edges$id[-in_group]
#   
#   if (!is.na(width[1])) { 
#     edges$width[1:nrow(edges)] <- NA
#     edges$width[in_group] <- width[1]
#     edges$width[out_group] <- width[2]
#   }
#   
#   if (!is.null(color[1])) { 
#     edges$width[1:nrow(edges)] <- NA
#     edges$width[in_group] <- color[1]
#     edges$width[out_group] <- color[2]
#   }
#   edges %>% select(-id)
# }
# subset_edges <- list(from = c(1, 114), to = c(42, 128))
# edges <- assign_attributes_to_edges(edges, subset_edges, width = c(5, 1), color = c("red", "blue"))
# 
# edges$id <- 1:nrow(edges)
# 
# ids <- edges %>% filter(
#     from %in% subset_edges$from, 
#     to %in% subset_edges$to
#   )
# 
# in_group <- ids$from
# out_group <- edges$id[-in_group]
# 
# range <- data.frame(region = c("A", "A", "B"), start = c(1, 3, 20), end = c(5, 10, 100))
# site <- data.frame(region = c("A", "A", "B"), site = c(4, 8, 25))
# 
# edges <- create_list_edges(A)
# visNetwork(nodes, edges, width = "100%", main = list(text = "title")) %>% 
# visIgraphLayout(layout = "layout_in_circle") %>%
#   visOptions(highlightNearest = list(enabled = T, hover = T))
# 
# 
# 
# 
# nodes <- create_list_nodes(nrow(A), labels = raw_input$gene_labels)
# 
# visnetwork <- function(nodes, 
#                        edges, 
#                        node_titles = 1:nrow(adj_matrix), 
#                        title = "", 
#                        igraph_layout = "layout_in_circle") { 
#   visNetwork(nodes, edges, width = "100%", main = list(text = title)) %>% 
#     visIgraphLayout(layout = igraph_layout) %>%
#     visOptions(highlightNearest = list(enabled = T, hover = T))
# }
# 
# visnetwork(nodes, edges)
# 
# 
#  
# 
# 
# # library(ggraph)
# # library(igraph)
# # 
# # library(network)
# # library(sna)
# # library(ggplot2)
# # 
# # 
# # library(GGally)
# # 
# # net = rgraph(10, mode = "graph", tprob = 0.5)
# # net = network(net, directed = FALSE)
# # 
# # create_network_object <- function(adj_matrix) { 
# #    network::as.network(x = as.matrix(adj_matrix), # the network object
# #                        directed = FALSE, # specify whether the network is directed
# #                        loops = FALSE, # do we allow self ties (should not allow them)
# #                        matrix.type = "adjacency") # the type of input
# # }
# # 
# # M = as.matrix(cvn$adj_matrices[[2]][[1]])
# # 
# # net = create_network_object(cvn$adj_matrices[[2]][[1]])
# # 
# # 
# # # vertex names
# # network.vertex.names(net) = letters[1:10]
# # 
# # set.edge.attribute(net, "core", c(rep(2, 5), rep(.5,5)))
# # set.edge.attribute(net, "color", c(rep("tomato", 5), rep("black",43-5)))
# # 
# # ggnet2(net, node.size = 12, 
# #        node.color = "black", 
# #        edge.size = 1, 
# #        label = TRUE, 
# #        label.color = "white", 
# #        label.size = 5, 
# #        mode = "circle", 
# #        edge.color = "color")
# # 
# # 
# # bip = data.frame(event1 = c(1, 2, 1, 0),
# #                  event2 = c(0, 0, 3, 0),
# #                  event3 = c(1, 1, 0, 4),
# #                  row.names = letters[1:4])
# # 
# # # weighted bipartite network
# # bip = network(bip,
# #               matrix.type = "bipartite",
# #               ignore.eval = FALSE,
# #               names.eval = "weights")
# # 
# # net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
# # 
# # 
# # ggraph(net) +
# #   geom_edge_link() +   # add edges to the plot
# #   geom_node_point()    # add nodes to the plot


nodes <- CVN::create_nodes_visnetwork(n_nodes = 5, labels = LETTERS[1:5])

adj_matrix <- matrix(c(0, 1, 0, 1, 0, 
                       1, 0, 1, 0, 0, 
                       0, 1, 0, 0, 0, 
                       1, 0, 0, 0, 1,
                       0, 0, 0, 1, 0), ncol = 5)

edges <- CVN::create_edges_visnetwork(adj_matrix)

edges <- set_attributes_to_edges_visnetwork(edges, 
                                   subset_edges = list(from = c(1, 2), to = c(4, 3)),
                                   width = c(3, .5), 
                                   color = c("red", "blue"))

CVN::visnetwork(nodes, edges)

cvn <- readRDS("../KiKme-analysis/results/KiKme-grid.rds")
cvn <- CVN::strip_cvn(cvn)

p <- CVN::visnetwork_cvn(cvn)
