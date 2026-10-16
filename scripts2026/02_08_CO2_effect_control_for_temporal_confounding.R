# This script is used to explore CO2's effect using spatial variation, because CO2 is strongly correlated to TIME in wet-cold forests. 
# Junna Wang, Sep 17, 2026
# 
rm(list=ls())

library(librarian)
shelf(tidyverse, ggpubr, lme4, dagitty, lmerTest, brms, randomForest, glmmTMB, performance, ranger, caret, Metrics, rsample, terra, car)

trend_driver <- read.csv(file.path("data2026", "potential_trend_driver.csv"))
sites_long <- read.csv(file.path('data2026', 'long_term_ec_sites.csv'))
sites_long <- sites_long %>% mutate(protocal = ifelse(source == "AmeriFlux_BASE", "Eddyproc", "Oneflux"))

trend_driver_site <- trend_driver %>% group_by(site_ID) %>% summarise(NEE_at_sd_min = min(NEE_at_sd, na.rm = T))
trend_driver <- left_join(trend_driver, trend_driver_site, by = "site_ID") %>% mutate(weight = trend_driver$nyear * (NEE_at_sd_min / NEE_at_sd)^2)   # use the same method to add weights.


flux_sites_good_CO2 <- c('CA-Ca1', 'CA-Oas', 'CA-Obs', "CH-Fru", "DE-Obe", "DE-RuR", "DE-Tha", "FI-Hyy",
                         "FI-Let", "FI-Sod", "FR-FBn", "IT-Ro2", "JP-Api", "JP-Tak", "US-GLE", "US-KL3", "US-KL2", "US-KM4",
                         "US-MOz", "US-PFa", "US-Rws", "US-Seg", "US-SRG", "US-SRM", "US-Tw1", "US-UMB", "US-UMd", "US-Whs", 
                         "AU-Whr", "US-Wjs", "ES-LM1", "ES-LM2", "JP-NsC", "JP-NsM", "NL-Vkp", "CH-Oe1", "CH-Oe3", "DE-SfS",
                         "CA-ARB", "CA-ARF")
flux_sites_partial_good_CO2 <- tribble(
  ~site_ID, ~years2remove, 
  "BR-Sa1", c(2020, 2021), 
  "CA-Cbo", 1996, 
  "CA-Gro", 2008, 
  "CA-KP1", 2008, 
  "CA-Man", 2007, 
  "CA-Mer", 1999, 
  "CA-Ojp", c(2002, 2003), 
  "CA-TP4", 2003, 
  "CA-TPD", 2012, 
  "DE-Gri", 2004,
  "DE-Hai", c(2013:2025),
  "ES-LJu", c(2017), 
  "GF-Guy", c(2017:2025),
  "IL-Yat", 2002,
  "JP-Fhk", c(2008, 2009, 2020), 
  "JP-Fjy", 2020, 
  "JP-Spp", c(2002, 2017), 
  "JP-Ta2", 2010,
  "JP-Tef", c(2021, 2022, 2023),  
  "RU-Fyo", 2017,
  "US-ChR", c(2015:2025),
  "US-Ha1", c(2006, 2007),
  "US-Ho1", c(2008, 2009, 2010), 
  "US-ICh", c(2021:2025), 
  "US-ICs", c(2013:2015, 2021), 
  "US-ICt", c(2021:2022), 
  "US-Jo2", 2012, 
  "US-KM3", 2012, 
  "US-Kon", 2007, 
  "US-LL1", c(2012, 2015), 
  "US-Los", 2025, 
  "US-MMS", c(1999, 2002), 
  "US-Mpj", c(2014:2016, 2018:2020), 
  "US-MtB", c(2012, 2013, 2021), 
  "US-NC2", 2009,
  "US-NR1", c(2004, 2005), 
  "US-Oho", c(2004, 2005), 
  "US-Ses", c(2011, 2016:2018), 
  "US-Syv", c(2000:2010, 2019), 
  "US-Ton", c(2004, 2005), 
  "US-Var", 2004,
  "US-Vcm", 2017, 
  "US-Vcp", c(2016, 2017), 
  "US-WCr", c(2021, 2025), 
  "US-Wjs", 2011, 
  "US-Wkg", 2007,
  "AU-Cpr", c(2014, 2018:2020), 
  "AU-Dry", c(2010:2016),
  "AU-Stp", c(2011, 2018), 
  "AU-Tum", c(2005, 2008),
  "AU-Whr", 2024,
  "CH-Rtb", 2022,
  "CN-Dzo", c(2011, 2012, 2017, 2018),
  "CZ-RAJ", c(2014, 2016:2018),
  "DE-Lnf", c(2003:2006, 2019:2024)
)

flux_sites_good_CO2 <- c(flux_sites_good_CO2, flux_sites_partial_good_CO2$site_ID)

trend_driver <- trend_driver[trend_driver$site_ID %in% flux_sites_good_CO2, ]


# adding IGBP, climate, and category into trend driver
if (!'category' %in% names(trend_driver)) {
  trend_driver <- trend_driver %>% left_join(sites_long[, c("site_ID", "IGBP", "Climate.class", "category", "NEEpartition", "protocal", "AI")], by = "site_ID")
}

#####--------
# we only use significant variables
climate_vegs <- c('cold_ENF', 'cold_DBF', 'cold_non_forest', 'dry', 'wet-warm')
year_groups <- list(2001:2005, 2006:2010, 2011:2014, 2015:2018)
#
# df_CO2_size = data.frame()
df_CO2_size <- as.data.frame(matrix(
  NA_real_,
  nrow=length(climate_vegs),
  ncol=2*length(year_groups)
))
rownames(df_CO2_size) <- climate_vegs
names(df_CO2_size) <- as.vector(rbind(
  paste0("estimate_",seq_along(year_groups)),
  paste0("se_",seq_along(year_groups))
))

n_type <- "rows"   # use "rows" or "sites"
min_n_rows <- 11L
min_n_sites <- 4L
min_n <- switch(
  n_type,
  rows=min_n_rows,
  sites=min_n_sites
)

df_CO2_n <- matrix(
  NA_integer_,
  nrow=length(climate_vegs),
  ncol=length(year_groups),
  dimnames=list(climate_vegs,NULL)
)

#
for (iclimate_veg in 1:length(climate_vegs)) {
  
  # iclimate_veg = 5
  
  climate_veg <- climate_vegs[iclimate_veg]
  #
  print(climate_veg)
  #
  if (climate_veg == 'cold_ENF') {
    data_group <- trend_driver %>% 
      # filter(IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold")
      filter(IGBP %in% c('ENF', 'EBF') & category == "cold")
  } else if (climate_veg == 'cold_DBF') {
    data_group <- trend_driver %>% 
      filter(IGBP %in% c('DBF', 'DNF', 'MF') & category == "cold")
  } else if (climate_veg == 'cold_non_forest') {
    data_group <- trend_driver %>% 
      filter(!IGBP %in% c('ENF', 'EBF', 'DBF', 'DNF', 'MF') & category == "cold") 
  } else if (climate_veg == 'dry') {
    data_group <- trend_driver %>% 
      filter(category == "dry") 
  } else if (climate_veg == 'wet-warm') {
    data_group <- trend_driver %>% 
      filter(category == "wet") 
  }
  
  # # calculate VIF between CO2 and TAgsT
  # print(vif(lm(NEE_at ~ CO2 + TA_gs_at, data = data_group)))
  # print(vif(lm(NEE_at ~ CO2 + year_mid, data = data_group)))
  # cor(data_group$CO2,  data_group$year_mid)^2
  
  
    # m_resid <- lm(CO2 ~ year_mid,
    #   data = data_group)
    # data_group$CO2_resid <- residuals(m_resid)

    
    for (j in 1:length(year_groups)) {
      
      year_group = year_groups[[j]]
      if (iclimate_veg %in% c(1, 3:5) && identical(as.integer(year_group), 2001:2005)) {
        next
      }
      
      data_sub = data_group %>% filter(year_mid %in% year_group)
      
      df_CO2_n[iclimate_veg,j] <- switch(
        n_type,
        rows=nrow(data_sub),
        sites=dplyr::n_distinct(data_sub$site_ID,na.rm=TRUE),
        stop("n_type must be 'rows' or 'sites'.")
      )
      
      #
      if (dplyr::n_distinct(data_sub$protocal, na.rm = TRUE) <= 1) {
        if (dplyr::n_distinct(data_sub$IGBP, na.rm = TRUE) <= 1) {
          mod <- lmer(
            NEE_at ~ CO2 + (1 | site_ID) + year_mid, 
            weights = weight, 
            data = data_sub)           
        } else {
          mod <- lmer(
            NEE_at ~ CO2 + (1 | site_ID) + year_mid, #  + (1 | IGBP)
            weights = weight, 
            data = data_sub)          
        }
      } else {
        if (dplyr::n_distinct(data_sub$IGBP, na.rm = TRUE) <= 1) {
          mod <- lmer(
            NEE_at ~ CO2 + (1 | site_ID) + year_mid, #  + (1 | protocal)
            weights = weight, 
            data = data_sub) 
          } else {
          mod <- lmer(
            NEE_at ~ CO2 + (1 | site_ID) + year_mid, #  + (1 | IGBP) + (1 | protocal) 
            weights = weight, 
            data = data_sub)
          }
      }
      print(summary(mod))
      coef_table <- summary(mod)$coefficients
      #
      df_CO2_size[iclimate_veg, 2*j - 1] <- coef_table["CO2", "Estimate"]
      df_CO2_size[iclimate_veg, 2*j] <- coef_table["CO2", "Std. Error"]
    }

    data_group %>% filter(between(year_mid, 2000, 2018)) %>%
      mutate(year_mid_factor = factor(year_mid)) %>%
    ggplot(aes(x = CO2, y = NEE_at, col = year_mid_factor)) +
      geom_point() +
      geom_smooth(method = "lm")
}
rownames(df_CO2_size) <- climate_vegs

########---------------------make a plot to show this result------------------------
# Period labels
period_labels <- vapply(
  year_groups,
  function(x) paste0(min(x), "\u2013", max(x)),
  character(1)
)

# Convert df_CO2_size from wide to long format
plot_data <- dplyr::bind_rows(
  lapply(seq_along(year_groups),function(j) {
    data.frame(
      climate_veg=rownames(df_CO2_size),
      period=period_labels[j],
      estimate=df_CO2_size[[2*j-1]],
      se=df_CO2_size[[2*j]],
      n=df_CO2_n[,j]
    )
  })
) %>%
  mutate(
    ci_lower=estimate-1.96*se,
    ci_upper=estimate+1.96*se,
    period=factor(period,levels=rev(period_labels)),
    climate_veg=factor(climate_veg,levels=climate_vegs),
    n_label=ifelse(is.na(n),"",paste0("n = ",n))
  )

plot_data <- plot_data %>%
  mutate(
    # Reverse effect direction
    estimate_plot = -estimate,
    
    # 95% confidence interval
    ci95_lower = estimate_plot - qnorm(0.975) * se,
    ci95_upper = estimate_plot + qnorm(0.975) * se,
    
    # 90% confidence interval
    ci90_lower = estimate_plot - qnorm(0.95) * se,
    ci90_upper = estimate_plot + qnorm(0.95) * se,
    
    significant_95 = ci95_lower > 0 | ci95_upper < 0,
    significant_90 = ci90_lower > 0 | ci90_upper < 0,
    
    significance = case_when(
      significant_95 ~ "*",
      significant_90 ~ "•",
      TRUE ~ ""
    ),
    
    line_color = if_else(
      significant_95,
      "Significant",
      "Not significant"
    )
  )

#---------------------
facet_order <- c(
  "cold_ENF",
  "cold_DBF",
  "cold_non_forest",
  "wet-warm",
  "dry"
)

facet_labels <- c(
  cold_ENF =
    "Wet-cold climate:\n ENF",
  cold_DBF =
    "Wet-cold climate:\n DBF+DNF+MF",
  cold_non_forest =
    "Wet-cold climate:\nnonforest ecosystems",
  `wet-warm` =
    "Wet-warm climate:\nall ecosystems",
  dry =
    "Dry climate:\nall ecosystems"
)

plot_data <- plot_data %>%
  mutate(
    climate_veg = factor(
      as.character(climate_veg),
      levels = facet_order
    )
  )

plot_data <- plot_data %>%
  mutate(
    show_result=is.finite(n) & n>=min_n,
    n_label=ifelse(show_result,paste0("n = ",n),"")
  )

plot_data_show <- plot_data %>%
  filter(show_result)


#------------------------------
# p_CO2 <- ggplot(plot_data, aes(y = period)) +
p_CO2 <- ggplot(plot_data_show,aes(y=period)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey50"
  ) +
  
  # Display the 95% confidence intervals
  geom_segment(
    aes(
      x = ci95_lower,
      xend = ci95_upper,
      yend = period,
      color = line_color
    ),
    linewidth = 1,
    na.rm = TRUE
  ) +
  
  geom_point(
    aes(x = estimate_plot),
    size = 3,
    shape = 21,
    fill = "white",
    color = "black",
    na.rm = TRUE
  ) +
  
  geom_text(
    aes(
      x=(ci95_lower+ci95_upper)/2,
      label=n_label
    ),
    nudge_y=0.25,
    size=3.2,
    vjust=0,
    na.rm=TRUE
  ) +
  
  # Significance mark above each effect-size point
  geom_text(
    aes(
      x=estimate_plot,
      label=significance
    ),
    nudge_y=-0.22,
    size=5,
    na.rm=TRUE
  ) +
  
  facet_wrap(
    ~ climate_veg,
    nrow = 1,
    scales = "free_x",
    labeller = labeller(climate_veg = as_labeller(facet_labels))
  ) +
  
  scale_color_manual(
    values = c(
      "Significant" = "#2878B5",
      "Not significant" = "grey70"
    ),
    guide = "none"
  ) +
  
  scale_y_discrete(drop = FALSE) +
  
  scale_x_continuous(
    expand=expansion(mult=0.15)
  ) +
  
  labs(
    x = expression(CO[2] ~ "total effect size based on spatial variation of "~CO[2]),
    y = "Period"
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90")
  )

p_CO2

# dir_fig <- '/Users/junnawang/YaleLab/manuscript2/figures'
# ggsave(file.path(dir_fig, "FigS_CO2_effect_spatial_variation.png"), p_CO2, width = 8.5, height = 4.5, dpi = 300)
