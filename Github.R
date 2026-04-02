#Needed libraries
library(tidyverse)
library(data.table)
library(wesanderson)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(minpack.lm)

# Parameters
setwd() #Directory of the file
types = c("static",'different','sealevel','upliftsealevel') #Treatments
file <- "fiji_auto_short_cline.txt" #Name of the file with the data
generations <- list(40,50,60,70,80,90,100,110,120) #list of generations to run through
diff_prob <- 0.002 # Dispersal probability
treatments_names <- list("Static", "Sealevel", "StableUplift", "DifferentUplift") # Names for use

#Functions
return_cline <- function(file, gen, type, diff_prob) {
  data <- read.table(file) %>%
    filter(V2==diff_prob) %>% # Important to put your column name and number correctly
    filter(V5==type) %>%
    filter(V1!=45)
  focal <- pull(data[gen+5]) %>%
    paste(., collapse="\n") %>%
    fread(fill=TRUE)
  focal[focal == -9.9] <- NA
  return(focal)
}

plot_mean_slope <- function(df) {
  treatments_names <- list("Static", "SeaLevel_Change", "StableUplift", "DifferentUplift") # 
  ggplot(df, aes(y = factor(t,treatments_names), x = A_Slope, fill = t)) +
    geom_boxplot() +
    
    facet_wrap(~ generation) +
    theme_bw() +
    labs(
      title = "Mean Slope by Treatment and Generation",
      y = "Treatment",
      x = "Mean Slope"
    ) +
    scale_fill_manual(values = wes_palette("Darjeeling1"),name = "Treatment") + 
    theme(axis.text.y = element_text(face="bold", color="black",  size=8))
}
plot_infl_point <- function(df) {
  df_clean <- df %>% filter(!is.na(i))
  treatments_names <- list("Static", "SeaLevel_Change", "StableUplift", "DifferentUplift") # 
  
  ggplot(df_clean, aes(y = factor(t,treatments_names), x = i, fill = t)) +
    geom_boxplot() +
    facet_wrap(~ generation) +
    theme_bw() +
    labs(
      title = "Location of Inflection Point by Treatment and Generation",
      y = "Treatment",
      x = "Inflection Point Location"
    ) +
    scale_fill_manual(values = wes_palette("Darjeeling1"),name = "Treatment")+
    theme(axis.text.y = element_text(face="bold", color="black",  size=8))+
    geom_vline(xintercept=25.5)
}

pointplot <- function(dataframe,generation_vec,type_vec, prob_vec,
                      max_bin=39,number_of_bins = 13){
  plot_df <- data.frame(name=character(), bins=integer(), Freq=integer(),
                        generation=integer(), type=character(), A_Slope=numeric())
  
  
  bin_size = floor(max_bin/number_of_bins)
  
  for (gen in generation_vec) {
    for (map_type in type_vec){
      
      df <- dataframe %>%
        drop_na() %>%
        filter(abs(generation-gen)<=10) %>%
        filter(generation==gen) %>%
        filter(t==map_type) %>%
        mutate(floor = floor((i)/bin_size))%>%
        filter(A_Slope >= 0)
      
      
      
      
      df_floor_counts <- df %>%
        count(floor, name = "Freq") %>%
        rename(bins = floor) %>%
        complete(bins = 0:number_of_bins, fill = list(Freq = 0))
      
      summ_floor <- df %>% 
        mutate(name = sprintf("%s_%f",map_type,gen),
               generation = gen,
               t = map_type) %>%
        rename(bins = floor)      # do NOT group yet
      
      df_completed <- summ_floor %>%
        ungroup() %>%             # ensure no grouping
        complete(
          bins = 0:number_of_bins,
          name = unique(name)
        )
      
      
      df_completed <- df_completed%>%
        left_join(df_floor_counts,by = 'bins')
      
      
      plot_df <- rbind(plot_df,df_completed)    
    }
    
  }
  treatments_names <- list("Static", "Sealevel", "StableUplift", "DifferentUplift")
  plot <- ggplot(na.omit(filter(plot_df, bins<number_of_bins)),aes(x = bins*3, y = generation)) + 
    geom_point(aes(size = Freq, alpha=Freq, fill=A_Slope), shape=21, stroke=1) +
    scale_fill_viridis_c(direction=-1, option = "viridis") +
    scale_size_continuous(range = c(2, 8)) +
    scale_y_continuous(expand = expansion(add = 10)) +
    theme_bw() + facet_grid(rows = vars(factor(t,treatments_names)), cols=vars(0)) + 
    ylab("Generations into simulation") + xlab("Transect Point") + 
    geom_vline(xintercept=25.5)+
    theme(
      axis.text = element_text(face = "bold", color = "black", size = 10),
      legend.position = "right"
    ) + theme(
      axis.title.x = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black") )
  return(plot)
}
Break_point_stick <- function(file, gen=60, types, diff_prob, treshold=0.0003) {
  
  # Get all the clines to get the slope
  Static <- return_cline(file, gen, types[1], diff_prob)
  StableUplift <- return_cline(file, gen, types[4], diff_prob)
  Sealevel <- return_cline(file, gen, types[3], diff_prob)
  DifferentUplift <- return_cline(file, gen, types[2], diff_prob)
  
  xnew <- seq(1, 39) # Create sequence of x transect values
  
  treatments <- list(Static, StableUplift, DifferentUplift,Sealevel) # List of treatment dataframes
  treatments_names <- list("Static", "StableUplift", "DifferentUplift","Sealevel") # 
  
  y <- 1 # Starting treatment
  
  # Create data frame to store the values after
  df <- data.frame(t=character(), replicate=integer(), i=integer(),
                   A_Slope=double(),
                   stringsAsFactors = FALSE)
  
  # For loop to get all the treatments
  for (treatment in treatments) {
    cur_treat <- treatments_names[[y]]
    y <- y + 1
    x <- 1 # Starting replicate
    
    for (row in seq_len(nrow(treatment))) {
      
      point = xnew
      val = as.numeric(unlist(treatment[row,1:39]))
      daf <- data.frame(point = point, val = val)
      fit <- try(
        nls(val ~ ifelse(point> breakpoint, slope * point + intercept, slope * breakpoint + intercept),
            data = daf,
            start= list(intercept=-50, slope=1, breakpoint=15),
            nls.control(maxiter = 5000, tol = 1e-01, minFactor = 1/1024,
                        printEval = FALSE)),
        silent = TRUE)
      
      if (inherits(fit, "try-error")) {
        break_point <- NA
        A_Slope <- NA
      } else {
        break_point <- summary(fit)$coefficients[3]
        A_Slope <- summary(fit)$coefficients[2]
        
      }
      
      
      
      
      
      
      
      
      # Append to the data frame
      new_row <- data.frame(
        t = cur_treat,
        replicate = x,
        i = break_point,
        A_Slope = A_Slope,
        stringsAsFactors = FALSE
      )
      df <- rbind(df, new_row)
      x <- x + 1
    }
  }
  
  return(df)
}

combine_multiple_conditions <- function(file, generations, types,
                                        diff_prob, thresholds=0.001) {
  # Create new data frame to store values
  combined_data <- data.frame()
  
  #Loop through generations and thresholds and store the values on the combined data frame
  for (gen in generations) {
    for (thresh in thresholds) {
      temp_data <- Break_point_stick(file, gen, types, diff_prob, thresh)
      temp_data$generation <- gen
      temp_data$threshold <- thresh
      combined_data <- rbind(combined_data, temp_data)
    }
  }
  return(combined_data)
}
plot_mean_inflection_over_time <- function(df, treatments_to_compare = c("Static", "StableUplift", "SeaLevel_Change")) {
  
  # Filter for the treatments you want to compare and clean data
  df_clean <- df %>% 
    filter(!is.na(i)) %>%
    filter(t %in% treatments_to_compare)
  
  # Calculate mean and standard error for each treatment and generation
  summary_data <- df_clean %>%
    group_by(generation, t) %>%
    summarise(
      mean_i = mean(i, na.rm = TRUE),
      se_i = sd(i, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Create the plot
  plot <- ggplot(summary_data, aes(x = generation, y = mean_i, color = t, fill = t)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_i - se_i, ymax = mean_i + se_i), 
                alpha = 0.2, color = NA) +
    geom_hline(yintercept = 25.5, linetype = "dashed", color = "black") +
    scale_color_manual(values = wes_palette("Darjeeling1"), name = "Treatment") +
    scale_fill_manual(values = wes_palette("Darjeeling1"), name = "Treatment") +
    theme_bw() +
    labs(
      x = "Generation",
      y = "Mean Inflection Point Location"
    ) +
    theme(
      axis.text = element_text(face = "bold", color = "black", size = 10),
      legend.position = "right"
    ) + theme(legend.position="none")+ theme(
      axis.title.x = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 16, color = "black"),
      axis.text.x  = element_text(size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black")
    )
  
  
  return(plot)
}

# Usage:
plot_mean_inflection_over_time(dataframe)

# Or specify which treatments to compare:
plot_mean_inflection_over_time(dataframe2, 
                               treatments_to_compare = c("Static", 'Sealevel',"StableUplift", "DifferentUplift"))
dataframe <-combine_multiple_conditions(file,generations,types,diff_prob,thresholds=0.01)
plot_mean_slope(dataframe)
plot_infl_point(dataframe)
pointplot(dataframe,generations,treatments_names,diff_prob)
plot_mean_slope_over_time <- function(df, treatments_to_compare = c("Static", "StableUplift", "DifferentUplift")) {
  
  # Filter for the treatments you want to compare and clean data
  df_clean <- df %>% 
    filter(!is.na(A_Slope)) %>%
    filter(t %in% treatments_to_compare) %>%
    filter(A_Slope >= 0)  # Filter positive slopes like in your pointplot function
  
  # Calculate mean and standard error for each treatment and generation
  summary_data <- df_clean %>%
    group_by(generation, t) %>%
    summarise(
      mean_slope = mean(A_Slope, na.rm = TRUE),
      se_slope = sd(A_Slope, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Create the plot
  plot <- ggplot(summary_data, aes(x = generation, y = mean_slope, color = t, fill = t)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_ribbon(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), 
                alpha = 0.2, color = NA) +
    scale_color_manual(values = wes_palette("Darjeeling1"), name = "Treatment") +
    scale_fill_manual(values = wes_palette("Darjeeling1"), name = "Treatment") +
    theme_bw() +
    labs(
      x = "Generation",
      y = "Mean Slope"
    ) +
    theme(
      axis.text = element_text(face = "bold", color = "black", size = 10),
      legend.position = "right"
    ) + theme(legend.position="none")+ theme(
      axis.title.x = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.text.x  = element_text(size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black")
    )
  
  return(plot)
}

# Usage:
plot_mean_slope_over_time(dataframe)

rplot_all_replicates_with_mean <- function(file, gen, type, diff_prob) {
  # Read and filter data
  data <- read.table(file) %>%
    filter(V5 == type) %>%
    filter(V2 == diff_prob) %>%
    filter(V1!=45)
  
  # Get number of replicates
  n_replicates <- nrow(data)
  
  # Initialize list to store all replicate data
  all_data <- list()
  
  # Loop through all replicates
  for (i in 1:n_replicates) {
    focal <- pull(data[i, ], gen + 3) %>%
      fread(text = ., fill = TRUE) %>%
      unlist(use.names = TRUE)
    
    focal[focal == -9.9] <- NA
    
    # Store data with replicate ID
    combined <- data.frame(
      Banded = focal[1:39],
      Transect = seq(1, 39, length.out = 39),
      Replicate = i
    )
    
    all_data[[i]] <- combined
  }
  
  # Combine all replicates into one dataframe
  all_combined <- bind_rows(all_data)
  
  # Calculate mean across replicates for each transect point
  mean_data <- all_combined %>%
    group_by(Transect) %>%
    summarise(Mean_Banded = mean(Banded, na.rm = TRUE), .groups = 'drop')
  
  # Plot
  plot <- ggplot() +
    geom_line(data = all_combined, 
              aes(x = Transect, y = Banded, group = Replicate),
              color = "grey70", alpha = 0.5) +
    geom_line(data = mean_data,
              aes(x = Transect, y = Mean_Banded),
              color = "#F98400FF", size = 1.2) +
    ylim(-0.05, 1.05) +
    xlim(0, 39) +
    labs(y = "Banded", x = "Transect") +
    theme_bw() + theme(
      axis.title.x = element_text(size = 18, color = "black"),
      axis.title.y = element_text(size = 18, color = "black"),
      axis.text.x  = element_text(size = 12, color = "black"),
      axis.text.y  = element_text(size = 12, color = "black")
    )
  
  
  return(plot)
}


# Usage:
rplot_all_replicates_with_mean(file, 100, "different", diff_prob)

