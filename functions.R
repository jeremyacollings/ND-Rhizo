
##### FUNCTIONS ######

# Contrast Distributions --------------------------------------------------

get_contrasts <- function(exp, cont, sp.dims){
  contrast <- exp - cont
  medians <- c(median(contrast, na.rm = TRUE), 
               apply(contrast, sp.dims, median, na.rm = TRUE))
  lowers <- c(quantile(contrast, na.rm = TRUE, 0.025),
              apply(contrast, sp.dims, quantile, 0.025, na.rm = TRUE))
  uppers <- c(quantile(contrast, na.rm = TRUE, 0.975),
              apply(contrast, sp.dims, quantile, 0.975, na.rm = TRUE))
  pd <- c(sum(contrast > 0, na.rm = TRUE, 0.025)/length(contrast),
          apply(contrast, sp.dims, function(x) sum(x > 0, na.rm = TRUE)/length(x)))
  cbind.data.frame(medians, lowers, uppers, pd)
}

# Translating Arrays and Dataframes ---------------------------------------

array2df <- function(metric_array, metric){
  temp <- melt(metric_array, value.name = metric) %>%
    rename("draw" = "Var1", "focal" = "Var2", "comp" = "Var3", 
           "rhizo" = "Var4", "nitro" = "Var5")
  
  temp$pair <- apply(temp[,c("focal", "comp")], 1, function(x) paste(sort(x), collapse = "."))
  temp <- temp[!duplicated(paste(temp$pair, temp$rhizo, temp$nitro, temp$draw)),]
  temp <- temp[which(temp$focal != temp$comp),]
  temp
}

# Network Diagrams --------------------------------------------------------

pseudo_log <- function(x){
  log((x/2) + sqrt((x/2)^2 + 1))/log(10)
}

make_matrix <- function(mod, rh, eu){
  temp <- mod %>%
    spread_draws(alpha[s,c,r,e]) %>%
    group_by(s,c,r,e) %>%
    summarise(med = median(alpha), 
              low = quantile(alpha, 0.025), 
              up = quantile(alpha, 0.975)) %>%
    filter(r == rh, e == eu) %>%
    ungroup() %>%
    select(med) %>%
    unlist() %>%
    unname()
  
  matrix(temp, nrow = 3, ncol = 3, 
         byrow = TRUE)
}

plot_network <- function(mat, spp){
  edges <- data.frame(
    from = rep(1:3, each = 3),
    to = rep(1:3, times = 3),
    weight = as.vector(mat)
  )
  
  graph <- as_tbl_graph(edges, directed = TRUE)
  
  ggraph(graph, layout = "circle") +
    # interspecific edges
    geom_edge_fan(
      aes(width = weight, color = weight),
      arrow = arrow(length = unit(2, 'mm'), type = "closed"),
      start_cap = circle(10, 'mm'),
      end_cap = circle(10, 'mm'),
      angle_calc = "along"
    ) +
    # intraspecific edges
    geom_edge_loop(
      aes(width = weight, color = weight),
      arrow = arrow(length = unit(2, 'mm'), type = "closed"),
      position = "identity",
      start_cap = circle(10, 'mm'),
      end_cap = circle(10, 'mm'),
      angle = 90
    ) +
    geom_node_point(size = 20, color = "grey") +
    geom_node_text(aes(label = spp), color = "black", size = 6) +
    scale_edge_width(range = c(0.8, 2.5)) +
    scale_edge_alpha(range = c(0.4, 1)) +
    theme_graph() +
    theme(
      plot.margin = margin(30, 30, 30, 30),
      legend.position = "none"
    ) +
    coord_fixed(clip = "off") + 
    scale_edge_color_gradient2(low = "#FA9F42", 
                               mid = "white", 
                               high = "#2B4162")
}

