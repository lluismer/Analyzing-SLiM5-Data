# Analyzing-SLiM5-Data

An R pipeline for analyzing and visualizing **hybrid zone clines** from [SLiM 5](https://messerlab.org/slim/) forward-time population genetics simulations.

Built for research on the **Fiji Whistler** (*Pachycephala graeffii*) hybrid zone, testing how geological history — island uplift and sea-level change — shapes where a hybrid zone sits and how sharp it becomes.

---

## Background

A **cline** is the spatial transition in allele (or phenotype) frequency across a landscape. Where two diverging populations meet, that transition forms a hybrid zone, and two properties of the cline carry most of the biological signal:

- **Inflection point** — *where* the transition happens along the transect
- **Slope** — *how sharply* it happens (steeper clines imply stronger selection or weaker dispersal)

This pipeline extracts both from simulation output and compares them across geological scenarios and across time.

## Simulated treatments

Each simulation is run under one of four geological histories:

| Treatment | Scenario |
|---|---|
| `Static` | No geological change — the control |
| `SeaLevel` | Sea level changes over time |
| `StableUplift` | Island uplift at a constant rate |
| `DifferentUplift` | Uplift at differing rates across the landscape |

Clines are sampled across a **39-point spatial transect**, with the geographic midpoint at position **25.5** (drawn as a reference line on every figure — deviation from it means the zone has moved).

## Method

Each replicate cline is fit with a **piecewise "hockey-stick" model** using non-linear least squares (`nls`):

```
frequency = slope × position + intercept,   for position > breakpoint
            slope × breakpoint + intercept, otherwise
```

Fitting returns two parameters per replicate — `breakpoint` (inflection point) and `slope` — which become the response variables for every downstream figure. Replicates that fail to converge are recorded as `NA` and dropped rather than being allowed to bias the summaries.

## Requirements

```r
install.packages(c("tidyverse", "data.table", "ggplot2",
                   "tidyr", "dplyr", "minpack.lm", "wesanderson"))
```

## Usage

Set the working directory and parameters at the top of `Main.R`:

```r
setwd("path/to/your/data")                 # directory containing the data file
file <- "fiji_auto_short_cline.txt"        # SLiM 5 output
generations <- list(40,50,60,70,80,90,100,110,120)
diff_prob <- 0.002                         # dispersal probability to analyze
```

Then build the results table once and reuse it for each figure:

```r
dataframe <- combine_multiple_conditions(file, generations, types, diff_prob)

plot_mean_slope(dataframe)                 # slope by treatment and generation
plot_infl_point(dataframe)                 # inflection point by treatment and generation
plot_mean_slope_over_time(dataframe)       # mean slope across generations, with SE
plot_mean_inflection_over_time(dataframe)  # zone movement across generations, with SE
pointplot(dataframe, generations, treatments_names, diff_prob)

# Raw clines for a single treatment/generation
rplot_all_replicates_with_mean(file, 100, "different", diff_prob)
```

## Functions

| Function | Purpose |
|---|---|
| `return_cline()` | Pulls the cline for one generation / treatment / dispersal rate |
| `Break_point_stick()` | Fits the piecewise model per replicate; returns breakpoint and slope |
| `combine_multiple_conditions()` | Runs the fit across all generations into one tidy data frame |
| `plot_mean_slope()` | Boxplot of cline steepness by treatment, faceted by generation |
| `plot_infl_point()` | Boxplot of zone location by treatment, faceted by generation |
| `plot_mean_slope_over_time()` | Mean slope trajectory with standard-error ribbon |
| `plot_mean_inflection_over_time()` | Mean zone location trajectory with standard-error ribbon |
| `pointplot()` | Bubble plot of binned zone position over time; size = replicate count, fill = slope |
| `rplot_all_replicates_with_mean()` | All replicate clines in grey with the mean overlaid |

## Input format

Whitespace-delimited SLiM 5 output, one row per replicate:

| Column | Meaning |
|---|---|
| `V1` | Replicate / run identifier |
| `V2` | Dispersal probability |
| `V5` | Treatment name |
| `V(gen+5)` | Cline values for that generation, as an embedded delimited string |

Missing data uses the sentinel value `-9.9`, which is converted to `NA` on read.

## License

GPL-3.0
