#Multiple GGridgeplots then combine them on one page using Cowplot
zero_ridges = ggplot(zero, aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 1.0, rel_min_height = 0.01) +
  xlab("Normalized fluorescence") +
  ylab("Generation") +
  ggtitle("Zero copy control") +
  theme_classic() +
  #scale_x_continuous("Normalized Fluorescence", limits=c(0.05,1.0), expand = c(0.01, 0), breaks = c(0.1, 0.55, 1.0)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.text = element_text(family="Arial", size = 12),#edit legend text font and size
    legend.title = element_blank(), #remove legend title
    legend.position = 'none', #remove the legend
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
one_ridges = ggplot(one, aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 1.0, rel_min_height = 0.01) +
  #xlab("normalized fluorescence") +
  ylab("Generation") +
  ggtitle("One copy control") +
  theme_classic() +
  scale_x_continuous(limits=c(0.0,2.5), breaks = c(0, 1, 2, 2.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
one_ridges

two_ridges = ggplot(two, aes(x = B2A_FSC, y = generation, fill = ..x.., height=..density..)) +
  geom_density_ridges_gradient(scale = 1.0, rel_min_height = 0.01) +
  #xlab("normalized fluorescence") +
  ylab("generation") +
  ggtitle("two copy control") +
  theme_classic() +
  scale_x_continuous(limits=c(0.0,2.5), breaks = c(0, 1, 2, 2.5)) +
  scale_y_discrete(expand = expansion(add = c(0.2, 1.0))) + #expands the graph space or else the top is cut off
  scale_fill_distiller(type = "seq", palette = 5, direction = 1, guide = "colourbar") + #makes it green
  theme(
    legend.position = 'none', #remove the legend
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(family="Arial", size = 10, color = "black"), #edit x-tick labels
    axis.text.y = element_text(family="Arial", size = 10, color = "black")
  )
two_ridges

cowplot::plot_grid(zero_ridges, one_ridges, two_ridges,labels = c("A"), nrow=1, align = "h")