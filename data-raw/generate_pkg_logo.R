# GENERATE PACKAGE LOGO =======================================================

# Libraries:
library(dplyr)        # For data manipulation
library(stringr)      # For string manipulation
library(tidyr)        # for data manipulation
library(purrr)        # For vectorised functions in dataframes
library(ggplot2)      # For plotting
library(RColorBrewer) # For colours in plots
library(hexSticker)   # For sticker creation
library(showtext)     # For custom font

# Loading Google fonts:
font_add_google("Inter", "inter", bold.wt = 800)
## Automatically use showtext to render text for future devices
showtext_auto()

# Set sig level:
alpha <- 0.05

## Load wrapper functions =====================================================

source("data-raw/compute_size_curves.R")
source("data-raw/compute_power_curves.R")

# Plot m = 27, n = 21, OR %in% 1 like in paper:

plot_m <- 27
plot_n <- 21
plor_or <- 1

# Testing set:
testing_set <- tibble(
  m = rep(plot_m, length(plor_or)),
  n = rep(plot_n, length(plor_or)),
  or = plor_or
)

testing_set <- testing_set |>
  mutate(method = "zoom")

#### Size plot ----------------------------------------------------------------

# Run size calculations:
size_data <- testing_set |>
  mutate(
    size = pmap(.l = list("m" = m, "n" = n, "alpha" = alpha, "or" = or),
                .f = compute_size_curves)
  ) |>
  unnest(size)

# Pivot and create nicer labels:
size_plot_data <- size_data |>
  pivot_longer(modified:conservative, names_to = "test", values_to = "pi1") |>
  mutate(
    test_label = case_when(
      test == "modified" ~ "Modified",
      test == "woolf" ~ "Woolf",
      test == "sas_freq" ~ "SAS PROC Freq",
      test == "conservative" ~ "Conservative"),
    test_label = factor(test_label, levels = c(
      "Modified", "SAS PROC Freq", "Woolf", "Conservative"
    )),
    or_label = paste0("\\textit{$H_{0} \\, : \\, e^{\\theta_0} = ", or, "\\, (OR\\, = \\,", or, ")$}")
  )

# Plot:

package_logo_plot <- ggplot(data = size_plot_data, aes(x = p0, y = pi1)) +
  geom_line(aes(colour = test_label, linewidth = test_label),
            show.legend = FALSE) +
  #scale_y_continuous(limits = c(-4, 4*alpha*100)) +
  #geom_hline(yintercept = alpha*100, linetype = 2, colour = "grey30",
  #           linewidth = 0.5) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    ylim = c(0, 1.05*alpha*100)
  ) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = c("0", "0.25", "0.5", "0.75", "1")
  ) +
  scale_colour_manual(
    values = c(
      Modified        = "white",
      `SAS PROC Freq` = "grey70",
      Woolf           = "grey70",
      Conservative    = "grey70"
    )
  ) +
  scale_linewidth_manual(
    values = c(
      Modified        = 0.8,
      `SAS PROC Freq` = 0.5,
      Woolf           = 0.5,
      Conservative    = 0.5
    )
  ) +
  theme_void()

# Turn into sticker:

package_logo <- sticker(package_logo_plot,
                        package="modifiedfisher",
                        p_color = "white",
                        p_family = "inter",
                        p_fontface = "bold",
                        p_x = 1, p_y = .5, p_size = 10,
                        s_x = 1, s_y = 1.1, s_width = 1.6, s_height = 1,
                        h_fill = RColorBrewer::brewer.pal(n = 5, "PuBuGn")[[5]],
                        h_color = "#00473B",
                        spotlight = FALSE,
                        filename = "man/figures/logo.png", dpi = 240)
