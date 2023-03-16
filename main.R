library(tidyverse)
library(reshape2)
library(ggraph)
library(tidygraph)

get_adjacent_weights <- function(mat, mat_row, mat_col, self_select_weight = 0, diagonals = FALSE) {

  n_rows <- dim(mat)[1] 
  n_cols <- dim(mat)[2]
  
  if(diagonals) {
    df <- expand.grid(row = c(mat_row, mat_row - 1, mat_row + 1),
                      col = c(mat_col, mat_col - 1, mat_col + 1))
  } else {
    df <- data.frame(row = c(mat_row, mat_row - 1, mat_row + 1, mat_row, mat_row),
                     col = c(mat_col, mat_col, mat_col, mat_col - 1, mat_col + 1))
  }
  
  df <- df %>%
    rowwise() %>%
    transmute(from_row = rep(mat_row),
           from_col = rep(mat_col),
           to_row = case_when(row == 0L ~ n_rows,
                              row == n_rows + 1 ~ 1L,
                           TRUE ~ as.integer(row)),
           to_col = case_when(col == 0L ~ n_cols,
                              col == n_cols + 1 ~ 1L,
                           TRUE ~ as.integer(col)),
           weight = mat[to_row, to_col])
  
  df$weight[1] <- self_select_weight
  
  return(df)
}


get_adjacency_table <- function(mat, self_select_weight = 0) {
  n_rows <- dim(mat)[1] 
  n_cols <- dim(mat)[2]
  
  if(length(self_select_weight) == 1) {
    self_select_weight <- matrix(self_select_weight, nrow = n_rows, ncol = n_cols)
  }
  
  adjacency_df_list <- c()
  i <- 1
  row <- 1
  col <- 1
  
  for(row in 1:n_rows) {
    for(col in 1:n_cols) {
      adjacency_df_list[[i]] <- get_adjacent_weights(mat, row, col, self_select_weight = self_select_weight[row, col])
      i <- i + 1
    }
    col <- 1
  }
  
  return(do.call("rbind", adjacency_df_list))
}

normalize_weights <- function(x) {
  x/sum(x)
}


random_walk <- function(weight_mat, steps = 1, self_select_weight = 0) {
  n_rows <- dim(weight_mat)[1] 
  n_cols <- dim(weight_mat)[2]
  
  if(length(self_select_weight) == 1) {
    self_select_weight <- matrix(self_select_weight, nrow = n_rows, ncol = n_cols)
  }
  
  start_row <- sample(1:n_rows, 1)
  start_col <- sample(1:n_cols, 1)
  start_prob <- 0 #weight_mat[start_row, start_col]
  
  step <- 2
  
  path_row <- c(start_row)
  path_col <- c(start_col)
  path_prob <- c(1/(n_rows*n_cols))
  
  max_steps <- steps + 1
  
  for(step in 2:max_steps) {
    adjacent_weights <- get_adjacent_weights(weight_mat, path_row[step - 1], path_col[step - 1], self_select_weight = self_select_weight[path_row[step - 1], path_col[step - 1]])
    
    next_step_index <- sample(nrow(adjacent_weights), 1, prob = normalize_weights(adjacent_weights$weight))
    
    path_row[step] <- adjacent_weights$to_row[next_step_index]
    path_col[step] <- adjacent_weights$to_col[next_step_index]
    path_prob[step] <- normalize_weights(adjacent_weights$weight)[next_step_index]
  }
  
  return(data.frame(row = path_row, col = path_col, prob = path_prob))
}

random_walk_graph <- function(weight_mat, steps = 1, self_select_weight = 0) {
  n_rows <- dim(weight_mat)[1] 
  n_cols <- dim(weight_mat)[2]
  
  if(length(self_select_weight) == 1) {
    self_select_weight <- matrix(self_select_weight, nrow = n_rows, ncol = n_cols)
  }
  
  start_row <- sample(1:n_rows, 1)
  start_col <- sample(1:n_cols, 1)
  start_prob <- 0 #weight_mat[start_row, start_col]
  
  step <- 1
  
  path_from_row <- c(start_row)
  path_from_col <- c(start_col)  
  path_to_row <- c()
  path_to_col <- c()
  path_prob <- c(1/(n_rows*n_cols))
  
  for(step in 1:steps) {
    if(step == 1) {
      
      adjacent_weights <- get_adjacent_weights(weight_mat, path_from_row[step], path_from_col[step], self_select_weight = self_select_weight[path_from_row[step], path_from_col[step]])
      next_step_index <- sample(nrow(adjacent_weights), 1, prob = normalize_weights(adjacent_weights$weight))
      
    } else {

      adjacent_weights <- get_adjacent_weights(weight_mat, path_to_row[step - 1], path_to_col[step - 1], self_select_weight = self_select_weight[path_to_row[step - 1], path_to_col[step - 1]])
      next_step_index <- sample(nrow(adjacent_weights), 1, prob = normalize_weights(adjacent_weights$weight))
      
      path_from_row[step] <- adjacent_weights$from_row[next_step_index]
      path_from_col[step] <- adjacent_weights$from_col[next_step_index]

    }
    path_to_row[step] <- adjacent_weights$to_row[next_step_index]
    path_to_col[step] <- adjacent_weights$to_col[next_step_index]
    path_prob[step] <- normalize_weights(adjacent_weights$weight)[next_step_index]
  }
  
  return(data.frame(from_row = path_from_row, from_col = path_from_col, to_row = path_to_row, to_col = path_to_col, prob = path_prob))
}


monte_carlo <- function(weight_mat, steps = 20, sims = 20, self_select_weight = 0) {
  walks <- c()
  for(sim in 1:sims) {
    walks[[sim]] <- random_walk(weight_mat, steps = steps, self_select_weight) %>%
      mutate(sim_index = as.factor(sim))
  }
  return(do.call("rbind", walks))
}

monte_carlo_graph <- function(weight_mat, steps = 20, sims = 20, self_select_weight = 0) {
  walks <- c()
  for(sim in 1:sims) {
    walks[[sim]] <- random_walk_graph(weight_mat, steps = steps, self_select_weight) %>%
      mutate(sim_index = as.factor(sim))
  }
  return(do.call("rbind", walks))
}

n_row <- 49
n_col <- 49

mat <- matrix(runif(n_row * n_col, min = 0, max = 1), nrow = n_row, ncol = n_col)

mat <- matrix(rep(1:7, 343), nrow = n_row, ncol = n_col)

random_walk_test <- random_walk(mat, steps = 100)

#adjacency_df <- get_adjacency_table(mat)

mc_test <- monte_carlo(mat, steps = 100, sims = 1)

mat_df <- melt(mat) %>%
  rename(row = Var1, col = Var2, weight = value)

ggplot(mat_df, aes(x = col, y = row)) + 
  geom_raster(aes(fill = weight)) + 
  #geom_line(data = mc_test, aes(x = col, y = row, color = sim_index), linewidth = 2, alpha = 0.5) +
  geom_segment(data = mc_test, aes(x = col, y = row, color = sim_index), linewidth = 2, alpha = 0.5) +
  scale_fill_gradient(low = "grey90", high = "blue") +
  labs(x = "Column", y = "Row", title = "Weight Matrix") +
  theme_bw() + theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.3),
                     axis.text.y = element_text(size = 9),
                     plot.title = element_text(size = 11)) +
  coord_fixed()

mc_test_graph_df <- monte_carlo_graph(mat, steps = 100, sims = 1)

mc_test_graph <- as_tbl_graph(mc_test_graph_df)

ggraph(mc_test_graph, layout = "grid") + 
  geom_edge_link(
    arrow = arrow(), 
    start_cap = circle(5, "mm"),
    end_cap = circle(5, "mm")
  ) + 
  geom_node_point(aes(colour = class), size = 8)
