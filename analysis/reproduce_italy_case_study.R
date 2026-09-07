############################################################
# Reproducibility script for spatconnect: COVID-19 incidence in Italy
# Case study: Italy, NUTS-3 regions, first and second waves
#
# Data:
#   COVID-19 European Regional Tracker
#   - 04_master/italy_data.dta
#   - 04_master/NUTS_POPULATION.dta
#
# Geometry:
#   GISCO / Eurostat NUTS 2016 GeoJSON, downloaded directly
#   without giscoR to avoid httr2 version conflicts.
#
# Main idea:
#   For each wave, define scalar-valued areal data x_i as cumulative COVID-19
#   incidence per 100,000 population over the wave period, and standardise the
#   resulting data vector to z_i.
#
#   For a generic analysis vector y (y = z for superlevel analysis and y = -z
#   for sublevel analysis), A_i is activation, M_i is the local topological
#   merging contribution, and L_i = A_i - M_i is the net local contribution
#   to integrated Betti-0.
#
#   The local decomposition processes positive active-side values from high to
#   low. When an entering area j touches a previously active component C, that
#   fusion contributes one merging unit. One half is assigned to j and the
#   other half is divided equally among the neighbours of j that belong to C.
#   The normalization is therefore component-specific; if j touches q_j active
#   components, its entry generates q_j merging units. This allocation preserves
#   the exact local-to-global activation--merging identity.
#
#   Local inference is formulated as a lower-tail conditional test for L_i.
#   Because A_i is fixed for the focal area under conditional relabelling, this
#   is exactly equivalent to an upper-tail conditional test for M_i.
#
#   Continuous maps display M_i because it has the most direct connectivity
#   interpretation. Activation--merging diagrams retain A_i and M_i, with
#   the diagonal M_i = A_i corresponding to L_i = 0.
#
#   Figure terminology uses "superlevel" and "sublevel" to identify whether
#   the construction is applied to z or to the sign-reversed data vector -z.
#   The mathematical labels A_i, M_i, L_i and p_i are otherwise unchanged.
#
# Required packages:
#   install.packages(c("sf", "spdep", "ggplot2", "RColorBrewer",
#                      "dplyr", "haven"))
#   Install spatconnect from the accompanying GitHub repository.
############################################################

library(sf)
library(spdep)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(haven)
library(spatconnect)

cat("spatconnect version:", as.character(utils::packageVersion("spatconnect")), "\n")

# Avoid spherical geometry issues for contiguity construction.
sf::sf_use_s2(FALSE)

############################################################
# 0. User settings
############################################################

country_code <- "IT"
country_name <- "Italy"
target_nuts_level <- 3
gisco_year <- 2016
population_year <- 2021

# Two illustrative epidemic waves.
# First wave: early national outbreak, spatially concentrated in Northern Italy.
# Second wave: autumn/winter resurgence, broader spread.
waves <- data.frame(
  wave_id = c("first_wave", "second_wave"),
  wave_label = c("First wave", "Second wave"),
  start_date = as.Date(c("2020-02-24", "2020-10-01")),
  end_date   = as.Date(c("2020-05-03", "2020-12-31")),
  stringsAsFactors = FALSE
)

# Global curve thresholds and permutations.
n_grid <- 80
n_perm_global <- 1999

# Local conditional permutation test for integrated L_i.
# For each focal area i, z_i is kept fixed and the remaining values are
# randomly reassigned over the graph. The lower-tail test for L_i is equivalent
# to the upper-tail test for M_i and detects connector roles.
# This is analogous in spirit to the conditional permutation approach used
# for Local Moran / LISA.
#
# The conditional test is substantially slower than full random relabelling:
# computational cost is roughly n_areas * n_perm_local_conditional
# local decompositions per side and wave.
n_perm_local_conditional <- 1999

# Local Moran / LISA permutations. Use the same Monte Carlo resolution as the
# proposed local conditional test.
n_perm_local_moran <- 1999

# Thresholds (in within-wave SD units) used in the descriptive connected-
# component figure. Each wave produces a 2 x 3 display: superlevel and sublevel
# graphs at the same lambda values. These can be changed without affecting any
# inferential result.
component_figure_lambdas <- c(0.50, 0.75, 1.00)

# Local significance level used to define relevant connector areas.
local_alpha <- 0.05

# Use "p" for unadjusted lower-tail p-values for L_i in the descriptive
# typology, while reporting FDR q-values separately. Use "fdr" for a
# stricter typology.
local_significance_rule <- "p"  # allowed: "p", "fdr"

# Empirical tail cut-off used in the descriptive activation--merging typology.
# Superlevel analysis: high activation is defined by the upper eta-tail of z.
# Sublevel analysis: high activation is defined by the upper eta-tail of -z,
# equivalently the lower eta-tail of z.
activation_tail_eta <- 0.15

# Queen contiguity is default because NUTS-3 polygons may meet at points.
queen_contiguity <- TRUE

# Paper-style figure output.
# EPS files are the primary output; PNG previews can be useful for inspection.
save_png_previews <- TRUE
paper_base_size <- 20

# Common dimensions for single-region maps. Keeping these fixed helps make
# exported EPS files visually comparable when included at the same LaTeX width.
single_map_width <- 9
single_map_height <- 7

# Map outlines. Black borders improve legibility when light Brewer colors are used.
map_border_color <- "black"
map_border_linewidth <- 0.18
connector_outline_color <- "black"
connector_outline_linewidth <- 0.85

# Output folder.
output_dir <- paste0(
  "covid_italy_nuts",
  target_nuts_level,
  "_two_waves_topological_outputs"
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
# 1. Basic utilities
############################################################

to_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
}

find_col <- function(dat, exact = character(0), patterns = character(0),
                     required = TRUE, label = "column") {
  nms <- names(dat)
  nms_lower <- tolower(nms)

  for (v in exact) {
    hit <- which(nms_lower == tolower(v))
    if (length(hit) > 0L) return(nms[hit[1]])
  }

  for (p in patterns) {
    hit <- grep(p, nms_lower, perl = TRUE)
    if (length(hit) > 0L) return(nms[hit[1]])
  }

  if (required) {
    stop("Could not identify ", label, ". Available names are:\n",
         paste(nms, collapse = ", "))
  }

  NA_character_
}

parse_date_flexible <- function(x) {
  if (inherits(x, "Date")) return(x)

  if (is.numeric(x)) {
    # Stata daily dates use 1960-01-01 as origin.
    y <- as.Date(x, origin = "1960-01-01")
    if (mean(!is.na(y)) > 0.8) return(y)
  }

  x <- as.character(x)

  formats <- c(
    "%Y-%m-%d",
    "%d/%m/%Y",
    "%m/%d/%Y",
    "%Y/%m/%d",
    "%d-%m-%Y",
    "%d%b%Y",
    "%d%b%y",
    "%d %b %Y",
    "%d %B %Y"
  )

  best <- rep(as.Date(NA), length(x))
  best_n <- -1

  for (f in formats) {
    y <- suppressWarnings(as.Date(x, format = f))
    n_ok <- sum(!is.na(y))
    if (n_ok > best_n) {
      best <- y
      best_n <- n_ok
    }
  }

  y0 <- suppressWarnings(as.Date(x))
  if (sum(!is.na(y0)) > best_n) {
    best <- y0
  }

  best
}

safe_slug <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

############################################################
# 2. Download tracker data and population denominators
############################################################

download_tracker_dta <- function(filename, out_dir) {
  urls <- c(
    paste0(
      "https://raw.githubusercontent.com/asjadnaqvi/",
      "COVID19-European-Regional-Tracker/master/04_master/",
      filename
    ),
    paste0(
      "https://github.com/asjadnaqvi/",
      "COVID19-European-Regional-Tracker/raw/master/04_master/",
      filename
    )
  )

  local_file <- file.path(out_dir, filename)
  last_error <- NULL

  if (!file.exists(local_file)) {
    ok_download <- FALSE

    for (u in urls) {
      message("\nTrying tracker file:")
      message(u)

      ok <- tryCatch({
        utils::download.file(u, destfile = local_file, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e) {
        last_error <<- conditionMessage(e)
        FALSE
      }, warning = function(w) {
        last_error <<- conditionMessage(w)
        FALSE
      })

      if (ok && file.exists(local_file) && file.info(local_file)$size > 0) {
        ok_download <- TRUE
        break
      }
    }

    if (!ok_download) {
      stop("Could not download ", filename, ". Last error: ", last_error)
    }
  }

  as.data.frame(haven::read_dta(local_file))
}

tracker <- download_tracker_dta("italy_data.dta", output_dir)

cat("\nItaly tracker variables:\n")
print(names(tracker))

pop_raw <- download_tracker_dta("NUTS_POPULATION.dta", output_dir)

cat("\nNUTS population variables:\n")
print(names(pop_raw))

############################################################
# 3. Harmonise COVID tracker data
############################################################

nuts_col <- find_col(
  tracker,
  exact = c("nuts_id", "NUTS_ID", "nuts", "NUTS", "id", "geo"),
  patterns = c("^nuts", "nuts.*id", "^geo$"),
  required = TRUE,
  label = "NUTS id column"
)

date_col <- find_col(
  tracker,
  exact = c("date", "Date", "time", "TIME"),
  patterns = c("^date$", "date", "time"),
  required = TRUE,
  label = "date column"
)

daily_cases_col <- find_col(
  tracker,
  exact = c("cases_daily", "daily_cases", "new_cases", "cases_new"),
  patterns = c("case.*daily", "daily.*case", "new.*case"),
  required = FALSE,
  label = "daily cases column"
)

cases_col <- find_col(
  tracker,
  exact = c("cases", "cum_cases", "cumulative_cases", "total_cases"),
  patterns = c("^cases$", "cum.*case", "total.*case"),
  required = FALSE,
  label = "cases column"
)

if (is.na(daily_cases_col) && is.na(cases_col)) {
  stop("Could not find either daily cases or cases column in italy_data.dta.")
}

tracker$NUTS_join <- as.character(tracker[[nuts_col]])
tracker$date_parsed <- parse_date_flexible(tracker[[date_col]])

if (all(is.na(tracker$date_parsed))) {
  stop("Could not parse dates from column: ", date_col)
}

dat_country <- tracker[
  substr(tracker$NUTS_join, 1, 2) == country_code,
]

if (target_nuts_level == 3) {
  dat_country <- dat_country[nchar(dat_country$NUTS_join) == 5, ]
}
if (target_nuts_level == 2) {
  dat_country <- dat_country[nchar(dat_country$NUTS_join) == 4, ]
}

if (nrow(dat_country) == 0L) {
  stop("No rows found for Italy at requested NUTS level.")
}

cat("\nCOVID data rows after country/NUTS filtering:", nrow(dat_country), "\n")
cat("Number of NUTS regions in COVID file:",
    length(unique(dat_country$NUTS_join)), "\n")
cat("Date range in COVID file:\n")
print(range(dat_country$date_parsed, na.rm = TRUE))

############################################################
# 4. Harmonise population denominators
############################################################

get_tracker_population_nuts <- function(
    pop_raw,
    country_code,
    nuts_level = 3,
    population_year = 2021
) {
  pop_nuts_col <- find_col(
    pop_raw,
    exact = c(
      "nuts_id", "NUTS_ID", "nuts", "NUTS", "id", "geo",
      "region", "Region", "code", "Code"
    ),
    patterns = c("^nuts", "nuts.*id", "^geo$", "region", "code"),
    required = TRUE,
    label = "NUTS id column in population file"
  )

  pop_year_col <- find_col(
    pop_raw,
    exact = c("year", "Year", "time", "TIME", "date", "Date"),
    patterns = c("year", "time", "date"),
    required = FALSE,
    label = "year column in population file"
  )

  pop_age_col <- find_col(
    pop_raw,
    exact = c("age", "Age", "AGE"),
    patterns = c("^age$"),
    required = FALSE,
    label = "age column in population file"
  )

  pop_sex_col <- find_col(
    pop_raw,
    exact = c("sex", "Sex", "SEX"),
    patterns = c("^sex$"),
    required = FALSE,
    label = "sex column in population file"
  )

  pop_value_col <- find_col(
    pop_raw,
    exact = c(
      "population", "Population", "POPULATION",
      "pop", "Pop", "POP",
      "value", "Value", "VALUE",
      "obs_value", "OBS_VALUE"
    ),
    patterns = c("population", "^pop$", "pop_", "_pop", "^value$", "obs.*value"),
    required = FALSE,
    label = "population value column in population file"
  )

  if (is.na(pop_value_col)) {
    candidate_cols <- names(pop_raw)
    exclude <- unique(c(
      pop_nuts_col, pop_year_col, pop_age_col, pop_sex_col,
      grep("name|label|level|id|code|nuts|geo|lon|lat|x_|y_",
           candidate_cols, ignore.case = TRUE, value = TRUE)
    ))
    candidate_cols <- setdiff(candidate_cols, exclude)

    numeric_score <- vapply(candidate_cols, function(v) {
      x <- suppressWarnings(as.numeric(pop_raw[[v]]))
      mean(is.finite(x))
    }, numeric(1))

    candidate_cols <- candidate_cols[numeric_score > 0.8]

    if (length(candidate_cols) == 0L) {
      stop("Could not identify a population value column.")
    }

    medians <- vapply(candidate_cols, function(v) {
      median(suppressWarnings(as.numeric(pop_raw[[v]])), na.rm = TRUE)
    }, numeric(1))

    pop_value_col <- candidate_cols[which.max(medians)]
  }

  pop <- pop_raw
  pop$NUTS_join <- as.character(pop[[pop_nuts_col]])
  pop$population_value_tmp <- to_num(pop[[pop_value_col]])

  pop <- pop[substr(pop$NUTS_join, 1, 2) == country_code, ]

  if (nuts_level == 3) pop <- pop[nchar(pop$NUTS_join) == 5, ]
  if (nuts_level == 2) pop <- pop[nchar(pop$NUTS_join) == 4, ]

  if (nrow(pop) == 0L) {
    stop("No population rows matched requested country/NUTS level.")
  }

  selected_year <- NA_integer_

  if (!is.na(pop_year_col)) {
    yy <- suppressWarnings(as.integer(pop[[pop_year_col]]))
    available_years <- sort(unique(yy[is.finite(yy)]))

    if (length(available_years) > 0L) {
      selected_year <- available_years[
        which.min(abs(available_years - population_year))
      ]
      pop <- pop[yy == selected_year, ]
    }
  }

  if (!is.na(pop_sex_col)) {
    sx <- toupper(as.character(pop[[pop_sex_col]]))
    total_codes <- c("T", "TOTAL", "TOTAL_SEX", "BOTH", "ALL")
    if (any(sx %in% total_codes)) {
      pop <- pop[sx %in% total_codes, ]
    }
  }

  age_total_filtered <- FALSE
  if (!is.na(pop_age_col)) {
    ag <- toupper(as.character(pop[[pop_age_col]]))
    total_codes <- c("TOTAL", "T", "ALL", "Y_TOTAL", "Y_GE0", "TOTAL_AGE")
    if (any(ag %in% total_codes)) {
      pop <- pop[ag %in% total_codes, ]
      age_total_filtered <- TRUE
    }
  }

  pop <- pop[
    is.finite(pop$population_value_tmp) &
      pop$population_value_tmp > 0,
  ]

  if (nrow(pop) == 0L) {
    stop("No positive population values after filtering.")
  }

  if (!is.na(pop_age_col) && !age_total_filtered) {
    out <- pop %>%
      group_by(NUTS_join) %>%
      summarise(
        population = sum(population_value_tmp, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    out <- pop %>%
      group_by(NUTS_join) %>%
      summarise(
        population = max(population_value_tmp, na.rm = TRUE),
        .groups = "drop"
      )
  }

  out$population_year <- selected_year

  out <- out[
    is.finite(out$population) &
      out$population > 0,
  ]

  if (nrow(out) == 0L) {
    stop("No usable population denominators after collapsing.")
  }

  out
}

pop_nuts <- get_tracker_population_nuts(
  pop_raw = pop_raw,
  country_code = country_code,
  nuts_level = target_nuts_level,
  population_year = population_year
)

cat("\nPopulation denominators loaded for", nrow(pop_nuts), "regions.\n")

dat_country <- dat_country %>%
  left_join(pop_nuts, by = "NUTS_join")

n_pop_matched <- length(unique(dat_country$NUTS_join[
  is.finite(dat_country$population) & dat_country$population > 0
]))

cat("Population matched for", n_pop_matched,
    "of", length(unique(dat_country$NUTS_join)), "COVID regions.\n")

if (n_pop_matched == 0L) {
  stop("Population did not match COVID NUTS codes.")
}

############################################################
# 5. Construct wave-level cumulative incidence
############################################################

looks_cumulative_by_region <- function(dat, nuts_col, date_col, value_col) {
  tmp <- dat[!is.na(dat[[nuts_col]]) & !is.na(dat[[date_col]]), ]
  tmp$value_num <- to_num(tmp[[value_col]])
  tmp <- tmp[is.finite(tmp$value_num), ]

  if (nrow(tmp) == 0L) return(FALSE)

  diffs <- tmp %>%
    arrange(.data[[nuts_col]], .data[[date_col]]) %>%
    group_by(.data[[nuts_col]]) %>%
    summarise(
      prop_non_decreasing = {
        d <- diff(value_num)
        if (length(d) == 0L) NA_real_ else mean(d >= -1e-8, na.rm = TRUE)
      },
      .groups = "drop"
    )

  median(diffs$prop_non_decreasing, na.rm = TRUE) > 0.95
}

make_wave_incidence <- function(dat, start_date, end_date,
                                nuts_col = "NUTS_join",
                                date_col = "date_parsed",
                                daily_cases_col = NA_character_,
                                cases_col = NA_character_,
                                population_col = "population") {
  dat <- dat[
    !is.na(dat[[date_col]]) &
      dat[[date_col]] <= end_date,
  ]

  if (nrow(dat) == 0L) {
    stop("No data available up to wave end date.")
  }

  if (!is.na(daily_cases_col)) {
    dat$case_value <- pmax(to_num(dat[[daily_cases_col]]), 0)

    win <- dat[
      dat[[date_col]] >= start_date &
        dat[[date_col]] <= end_date,
    ]

    out <- win %>%
      group_by(.data[[nuts_col]]) %>%
      summarise(
        cases_wave = sum(case_value, na.rm = TRUE),
        population = max(.data[[population_col]], na.rm = TRUE),
        n_days = n_distinct(.data[[date_col]]),
        .groups = "drop"
      )

    names(out)[1] <- "NUTS_join"
    out$incidence_per_100k <- 100000 * out$cases_wave / out$population
    out$metric_source <- "sum_daily_cases_over_wave_div_population"
    return(out)
  }

  if (!is.na(cases_col)) {
    dat$case_value <- to_num(dat[[cases_col]])

    is_cumulative <- looks_cumulative_by_region(
      dat = dat,
      nuts_col = nuts_col,
      date_col = date_col,
      value_col = cases_col
    )

    if (is_cumulative) {
      end_dat <- dat %>%
        filter(.data[[date_col]] <= end_date) %>%
        arrange(.data[[nuts_col]], .data[[date_col]]) %>%
        group_by(.data[[nuts_col]]) %>%
        summarise(
          cases_end = dplyr::last(case_value),
          population = dplyr::last(.data[[population_col]]),
          .groups = "drop"
        )
      names(end_dat)[1] <- "NUTS_join"

      before_dat <- dat %>%
        filter(.data[[date_col]] < start_date) %>%
        arrange(.data[[nuts_col]], .data[[date_col]]) %>%
        group_by(.data[[nuts_col]]) %>%
        summarise(
          cases_before = dplyr::last(case_value),
          .groups = "drop"
        )
      names(before_dat)[1] <- "NUTS_join"

      out <- end_dat %>%
        left_join(before_dat, by = "NUTS_join")

      out$cases_before[is.na(out$cases_before)] <- 0
      out$cases_wave <- pmax(out$cases_end - out$cases_before, 0)
      out$n_days <- as.integer(end_date - start_date + 1)
      out$incidence_per_100k <- 100000 * out$cases_wave / out$population
      out$metric_source <- "cumulative_cases_difference_over_wave_div_population"
      return(out)
    }

    win <- dat[
      dat[[date_col]] >= start_date &
        dat[[date_col]] <= end_date,
    ]

    out <- win %>%
      group_by(.data[[nuts_col]]) %>%
      summarise(
        cases_wave = sum(pmax(case_value, 0), na.rm = TRUE),
        population = max(.data[[population_col]], na.rm = TRUE),
        n_days = n_distinct(.data[[date_col]]),
        .groups = "drop"
      )
    names(out)[1] <- "NUTS_join"
    out$incidence_per_100k <- 100000 * out$cases_wave / out$population
    out$metric_source <- "detected_daily_cases_over_wave_div_population"
    return(out)
  }

  stop("Could not construct wave incidence.")
}

wave_incidence_list <- list()

for (k in seq_len(nrow(waves))) {
  inc_k <- make_wave_incidence(
    dat = dat_country,
    start_date = waves$start_date[k],
    end_date = waves$end_date[k],
    daily_cases_col = daily_cases_col,
    cases_col = cases_col,
    population_col = "population"
  )

  inc_k$wave_id <- waves$wave_id[k]
  inc_k$wave_label <- waves$wave_label[k]
  inc_k$wave_start <- waves$start_date[k]
  inc_k$wave_end <- waves$end_date[k]

  inc_k <- inc_k[
    is.finite(inc_k$incidence_per_100k) &
      inc_k$incidence_per_100k >= 0 &
      is.finite(inc_k$population) &
      inc_k$population > 0,
  ]

  wave_incidence_list[[waves$wave_id[k]]] <- as.data.frame(inc_k)

  cat("\n", waves$wave_label[k], " incidence summary:\n", sep = "")
  print(summary(inc_k$incidence_per_100k))
  cat("Metric source:\n")
  print(table(inc_k$metric_source))
}

wave_incidence <- bind_rows(wave_incidence_list)

############################################################
# 6. Download GISCO NUTS geometries without giscoR
############################################################

download_gisco_nuts_geojson <- function(
    country_code,
    nuts_level = 3,
    gisco_year = 2016,
    out_dir = tempdir()
) {
  level_tag <- paste0("LEVL_", nuts_level)

  urls <- c(
    paste0(
      "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/",
      "NUTS_RG_20M_", gisco_year, "_4326_", level_tag, ".geojson"
    ),
    paste0(
      "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/",
      "NUTS_RG_20M_", gisco_year, "_4326.geojson"
    ),
    paste0(
      "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/",
      "NUTS_RG_10M_", gisco_year, "_4326_", level_tag, ".geojson"
    ),
    paste0(
      "https://gisco-services.ec.europa.eu/distribution/v2/nuts/geojson/",
      "NUTS_RG_10M_", gisco_year, "_4326.geojson"
    )
  )

  last_error <- NULL

  for (u in urls) {
    local_file <- file.path(out_dir, basename(u))

    message("\nTrying GISCO NUTS GeoJSON:")
    message(u)

    ok_download <- tryCatch({
      utils::download.file(u, destfile = local_file, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      FALSE
    }, warning = function(w) {
      last_error <<- conditionMessage(w)
      FALSE
    })

    if (!ok_download) next
    if (!file.exists(local_file) || file.info(local_file)$size == 0) next

    geom <- tryCatch({
      sf::st_read(local_file, quiet = TRUE)
    }, error = function(e) {
      last_error <<- conditionMessage(e)
      NULL
    })

    if (is.null(geom) || nrow(geom) == 0L) next

    geom <- sf::st_as_sf(geom)

    if (!("NUTS_ID" %in% names(geom))) {
      id_candidates <- c("nuts_id", "NUTSID", "id", "geo")
      hit <- id_candidates[tolower(id_candidates) %in% tolower(names(geom))]
      if (length(hit) > 0L) {
        names(geom)[tolower(names(geom)) == tolower(hit[1])] <- "NUTS_ID"
      }
    }

    if (!("NUTS_ID" %in% names(geom))) {
      last_error <- "Downloaded geometry has no NUTS_ID-like column."
      next
    }

    if ("CNTR_CODE" %in% names(geom)) {
      geom <- geom[geom$CNTR_CODE == toupper(country_code), ]
    } else {
      geom <- geom[substr(as.character(geom$NUTS_ID), 1, 2) == toupper(country_code), ]
    }

    if ("LEVL_CODE" %in% names(geom)) {
      geom <- geom[as.integer(geom$LEVL_CODE) == as.integer(nuts_level), ]
    } else {
      if (nuts_level == 3) geom <- geom[nchar(as.character(geom$NUTS_ID)) == 5, ]
      if (nuts_level == 2) geom <- geom[nchar(as.character(geom$NUTS_ID)) == 4, ]
    }

    geom <- geom[!sf::st_is_empty(geom), ]

    if (nrow(geom) > 0L) {
      message("Loaded ", nrow(geom), " NUTS geometries.")
      return(geom)
    }

    last_error <- "Geometry downloaded but no rows matched country/level."
  }

  stop(
    "Could not download usable GISCO NUTS geometries. Last error: ",
    last_error
  )
}

nuts_geom <- download_gisco_nuts_geojson(
  country_code = country_code,
  nuts_level = target_nuts_level,
  gisco_year = gisco_year,
  out_dir = output_dir
)

nuts_geom <- sf::st_as_sf(nuts_geom)
nuts_geom <- nuts_geom[!sf::st_is_empty(nuts_geom), ]
nuts_geom <- sf::st_make_valid(nuts_geom)

if (!("NUTS_ID" %in% names(nuts_geom))) {
  stop("GISCO geometry does not contain a NUTS_ID column.")
}

# Base geometry for all waves: keep only regions available in at least one wave.
nuts_ids_needed <- unique(wave_incidence$NUTS_join)

areal_base <- nuts_geom %>%
  inner_join(
    data.frame(NUTS_join = nuts_ids_needed),
    by = c("NUTS_ID" = "NUTS_join")
  )

if (nrow(areal_base) == 0L) {
  stop("No joined geometry rows. Check NUTS version/code compatibility.")
}

areal_base <- sf::st_transform(areal_base, 3035)
areal_base <- sf::st_make_valid(areal_base)
rownames(areal_base) <- NULL

# Use one fixed, slightly padded spatial extent for every map. This prevents
# small differences in panel scaling across waves and variables when EPS files
# are placed side by side in LaTeX.
make_padded_bbox <- function(x, pad_fraction = 0.015) {
  bb <- sf::st_bbox(x)
  x_pad <- as.numeric(bb["xmax"] - bb["xmin"]) * pad_fraction
  y_pad <- as.numeric(bb["ymax"] - bb["ymin"]) * pad_fraction

  c(
    xmin = as.numeric(bb["xmin"] - x_pad),
    xmax = as.numeric(bb["xmax"] + x_pad),
    ymin = as.numeric(bb["ymin"] - y_pad),
    ymax = as.numeric(bb["ymax"] + y_pad)
  )
}

common_map_bbox <- make_padded_bbox(areal_base)

coord_sf_common <- function() {
  ggplot2::coord_sf(
    xlim = common_map_bbox[c("xmin", "xmax")],
    ylim = common_map_bbox[c("ymin", "ymax")],
    datum = NA,
    expand = FALSE,
    clip = "on"
  )
}

cat("\nGeometry units after joining:", nrow(areal_base), "\n")

cat("\nBuilding contiguity graph...\n")

nb <- spdep::poly2nb(
  areal_base,
  queen = queen_contiguity,
  snap = sqrt(.Machine$double.eps)
)

nb_to_edges <- function(nb) {
  edge_list <- list()
  counter <- 1L

  for (i in seq_along(nb)) {
    js <- nb[[i]]

    if (length(js) == 0L) next
    if (length(js) == 1L && js[1] == 0L) next

    js <- js[js > i]

    if (length(js) > 0L) {
      edge_list[[counter]] <- data.frame(i = i, j = js)
      counter <- counter + 1L
    }
  }

  if (length(edge_list) == 0L) {
    stop("No edges found in neighbour list.")
  }

  edges <- do.call(rbind, edge_list)
  edges$i <- as.integer(edges$i)
  edges$j <- as.integer(edges$j)
  rownames(edges) <- NULL
  edges
}

edges <- nb_to_edges(nb)

cat("Number of areal units:", nrow(areal_base), "\n")
cat("Number of undirected edges:", nrow(edges), "\n")
cat("Isolated units:", sum(spdep::card(nb) == 0L), "\n")

############################################################
# 7. Connectivity methodology from spatconnect
############################################################

# The global Betti-0 test and the component-specific local activation--merging
# decomposition are intentionally NOT reimplemented in this script. They are
# called from spatconnect so that the reproducibility analysis and the public
# package use the same computational engine.
#
# Global analysis: spatconnect::global_connectivity()
# Local analysis:  spatconnect::local_connectivity()
#
# The paper case study uses 80 equally spaced non-negative thresholds for the
# global curves, 1999 global relabelings, and 1999 conditional local relabelings
# per focal area.

############################################################
# 8. Plotting and case-study utilities
############################################################


# Brewer-based palettes.
set1 <- RColorBrewer::brewer.pal(9, "Set1")
blues <- RColorBrewer::brewer.pal(9, "Blues")
ylorrd <- RColorBrewer::brewer.pal(9, "YlOrRd")
rdbu <- RColorBrewer::brewer.pal(11, "RdBu")
ylgnbu <- RColorBrewer::brewer.pal(9, "YlGnBu")
brbg <- RColorBrewer::brewer.pal(11, "BrBG")

# Topological relevance colors: keep the same semantic color mapping across
# superlevel and sublevel profiles. Peak and trough connectors use the same
# color because they occupy the same algebraic role in the two tails.
activation_merging_colors_pos <- c(
  "Core connector" = brbg[11],
  "Bridge connector" = brbg[9],
  "Peak connector" = brbg[3],
  "Not significant" = "grey90"
)

activation_merging_colors_neg <- c(
  "Core connector" = brbg[11],
  "Bridge connector" = brbg[9],
  "Trough connector" = brbg[3],
  "Not significant" = "grey90"
)

activation_merging_labels_pos <- c(
  "Core connector" = "Core connector",
  "Bridge connector" = "Bridge connector",
  "Peak connector" = "Peak connector",
  "Not significant" = "Not significant"
)

activation_merging_labels_neg <- c(
  "Core connector" = "Core connector",
  "Bridge connector" = "Bridge connector",
  "Trough connector" = "Trough connector",
  "Not significant" = "Not significant"
)

# Local Moran / LISA colors. High-Low and Low-High are deliberately shown
# using a warm light pink and a light blue. The pink is deliberately kept
# away from violet tones.
lisa_colors <- c(
  "High-High" = set1[1],
  "Low-Low" = set1[2],
  "High-Low" = "#F4A3A8",
  "Low-High" = "#9ECAE1",
  "Not significant" = "grey90"
)

scatter_role_levels_pos <- c(
  "Core connector",
  "Bridge connector",
  "Peak connector",
  "Not significant"
)

scatter_role_levels_neg <- c(
  "Core connector",
  "Bridge connector",
  "Trough connector",
  "Not significant"
)

scatter_role_colors_pos <- c(
  "Core connector" = brbg[11],
  "Bridge connector" = brbg[9],
  "Peak connector" = brbg[3],
  "Not significant" = "grey90"
)

scatter_role_colors_neg <- c(
  "Core connector" = brbg[11],
  "Bridge connector" = brbg[9],
  "Trough connector" = brbg[3],
  "Not significant" = "grey90"
)

# Explicit guides for activation--merging scatter plots. Without setting guide
# orders, ggplot2 can place the size legend and the colour legend in different
# orders depending on the plotted data and on the device. These guides keep the
# connectivity-role legend first and the conditional p-value legend second in
# all waves and tails.
stable_scatter_guides <- function() {
  guides(
    color = guide_legend(
      order = 1,
      override.aes = list(size = 5, alpha = 1)
    ),
    size = guide_legend(order = 2)
  )
}

paper_theme <- function(base_size = paper_base_size) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size = base_size - 2,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      plot.caption = element_text(size = base_size - 3, hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 2),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.key.height = grid::unit(0.52, "cm"),
      legend.key.width = grid::unit(0.52, "cm"),
      strip.text = element_text(face = "bold", size = base_size),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.margin = margin(8, 10, 8, 10)
    )
}

paper_map_theme <- function(base_size = paper_base_size) {
  theme_void(base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = base_size + 2,
        hjust = 0.5,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size = base_size - 2,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      plot.caption = element_text(size = base_size - 3, hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.key.height = grid::unit(0.52, "cm"),
      legend.key.width = grid::unit(0.52, "cm"),
      strip.text = element_text(face = "bold", size = base_size),
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "center",
      plot.margin = margin(8, 10, 8, 10)
    )
}

continuous_fill_guide <- function() {
  guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barheight = grid::unit(3.8, "cm"),
    barwidth = grid::unit(0.35, "cm")
  )
}

force_discrete_fill_legend <- function(values) {
  # For geom_sf(), drop = FALSE can leave empty legend boxes for categories
  # absent from a particular map. Explicitly overriding the legend key fill
  # forces all categories to appear with their assigned colors.
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      byrow = TRUE,
      override.aes = list(
        fill = unname(values),
        colour = rep(map_border_color, length(values)),
        linewidth = rep(0.35, length(values)),
        alpha = rep(1, length(values))
      )
    )
  )
}

add_wave_suffix_to_figure_filename <- function(filename, wave_dir) {
  # Per-wave figures are saved inside first_wave/ and second_wave/.
  # Add a short suffix before the extension so that exported figures have
  # distinct names even when copied to a common folder.
  wave_folder <- basename(normalizePath(wave_dir, mustWork = FALSE))

  wave_suffix <- if (wave_folder == "first_wave") {
    "w1"
  } else if (wave_folder == "second_wave") {
    "w2"
  } else {
    NA_character_
  }

  if (is.na(wave_suffix)) return(filename)

  ext <- tools::file_ext(filename)

  if (nzchar(ext)) {
    stem <- sub(paste0("\\.", ext, "$"), "", filename)
    if (grepl(paste0("_", wave_suffix, "$"), stem)) return(filename)
    return(paste0(stem, "_", wave_suffix, ".", ext))
  }

  if (grepl(paste0("_", wave_suffix, "$"), filename)) return(filename)
  paste0(filename, "_", wave_suffix)
}

save_plot <- function(
    p,
    filename,
    wave_dir = output_dir,
    width = 8,
    height = 6,
    save_png = save_png_previews
) {
  filename <- add_wave_suffix_to_figure_filename(filename, wave_dir)
  path <- file.path(wave_dir, filename)
  ext <- tolower(tools::file_ext(path))

  if (ext == "") {
    path <- paste0(path, ".eps")
    ext <- "eps"
  }

  if (ext == "eps") {
    eps_device <- function(filename, width, height, ...) {
      grDevices::cairo_ps(
        filename = filename,
        width = width,
        height = height,
        onefile = FALSE,
        fallback_resolution = 600,
        ...
      )
    }

    ggsave(
      path,
      p,
      width = width,
      height = height,
      units = "in",
      device = eps_device,
      bg = "white",
      limitsize = FALSE
    )

    if (isTRUE(save_png)) {
      png_path <- sub("\\.eps$", ".png", path, ignore.case = TRUE)
      ggsave(
        png_path,
        p,
        width = width,
        height = height,
        units = "in",
        dpi = 400,
        bg = "white",
        limitsize = FALSE
      )
      message("Saved preview: ", png_path)
    }
  } else {
    ggsave(
      path, p,
      width = width, height = height, units = "in", dpi = 400,
      bg = "white", limitsize = FALSE
    )
  }

  message("Saved: ", path)
  invisible(path)
}

plot_sf_map <- function(
    dat,
    var,
    title,
    subtitle = NULL,
    fill_label = NULL,
    scale_type = c("sequential", "diverging", "pvalue")
) {
  scale_type <- match.arg(scale_type)

  if (is.null(fill_label)) {
    fill_label <- var
  }

  p <- ggplot(dat) +
    geom_sf(
      aes(fill = .data[[var]]),
      color = map_border_color,
      linewidth = map_border_linewidth
    ) +
    coord_sf_common() +
    labs(
      title = title,
      subtitle = subtitle,
      fill = fill_label
    ) +
    paper_map_theme()

  if (scale_type == "sequential") {
    p <- p +
      scale_fill_gradientn(
        colors = ylorrd,
        na.value = "grey90",
        guide = continuous_fill_guide()
      )
  }

  if (scale_type == "diverging") {
    lim <- max(abs(dat[[var]]), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 1

    p <- p +
      scale_fill_gradient2(
        low = rdbu[11],
        mid = rdbu[6],
        high = rdbu[1],
        midpoint = 0,
        limits = c(-lim, lim),
        na.value = "grey90",
        guide = continuous_fill_guide()
      )
  }

  if (scale_type == "pvalue") {
    p <- p +
      scale_fill_gradientn(
        colors = rev(ylgnbu),
        limits = c(0, 1),
        breaks = c(0, 0.05, 0.10, 0.50, 1.00),
        na.value = "grey90",
        guide = continuous_fill_guide()
      )
  }

  p
}

plot_global_beta0 <- function(global_res, title = "Global Betti-0 connectivity curve", subtitle = NULL) {
  df <- data.frame(
    lambda = global_res$curve$lambda,
    beta0_obs = global_res$curve$observed,
    beta0_mean = global_res$curve$relabel_mean,
    beta0_q025 = global_res$curve$relabel_lower,
    beta0_q975 = global_res$curve$relabel_upper
  )

  ggplot(df, aes(x = lambda)) +
    geom_ribbon(
      aes(ymin = beta0_q025, ymax = beta0_q975),
      fill = "grey85",
      color = NA
    ) +
    geom_line(
      aes(y = beta0_mean, linetype = "Random-relabeling mean"),
      linewidth = 0.95,
      color = "grey30"
    ) +
    geom_line(
      aes(y = beta0_obs, linetype = "Observed"),
      linewidth = 1.15,
      color = "black"
    ) +
    scale_linetype_manual(
      values = c("Observed" = "solid", "Random-relabeling mean" = "dashed"),
      breaks = c("Observed", "Random-relabeling mean")
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      caption = NULL,
      x = expression(lambda ~ "(SD units)"),
      y = expression(beta[0](lambda)),
      linetype = NULL
    ) +
    paper_theme()
}


.case_make_dsu <- function(n) {
  parent <- seq_len(n)
  rank <- integer(n)

  find <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }

  union <- function(a, b) {
    ra <- find(a)
    rb <- find(b)
    if (ra == rb) return(FALSE)
    if (rank[ra] < rank[rb]) {
      tmp <- ra; ra <- rb; rb <- tmp
    }
    parent[rb] <<- ra
    if (rank[ra] == rank[rb]) rank[ra] <<- rank[ra] + 1L
    TRUE
  }

  list(find = find, union = union)
}

connected_components_at_threshold <- function(y, edges, lambda) {
  # Connected components of the induced graph on {i: y_i >= lambda}.
  y <- as.numeric(y)
  n <- length(y)
  active <- is.finite(y) & y >= lambda
  component <- rep(NA_integer_, n)

  active_idx <- which(active)
  if (length(active_idx) == 0L) {
    return(list(
      active = active,
      component = component,
      n_active = 0L,
      beta0 = 0L
    ))
  }

  dsu <- .case_make_dsu(n)

  if (nrow(edges) > 0L) {
    for (r in seq_len(nrow(edges))) {
      i <- as.integer(edges$i[r])
      j <- as.integer(edges$j[r])

      if (active[i] && active[j]) {
        dsu$union(i, j)
      }
    }
  }

  roots <- vapply(active_idx, dsu$find, integer(1))
  root_levels <- unique(roots)
  component[active_idx] <- match(roots, root_levels)

  list(
    active = active,
    component = component,
    n_active = length(active_idx),
    beta0 = length(root_levels)
  )
}

make_threshold_component_figure <- function(
    areal,
    z,
    edges,
    lambda_values = component_figure_lambdas,
    wave_label = NULL
) {
  # Descriptive figure showing connected components of the thresholded graphs.
  # To keep the figure visually clean, the rows indicate whether the analysis is
  # superlevel or sublevel, and the columns indicate the common threshold value.
  # The same component palette is reused in every panel, so Component 1, 2, 3,
  # ... always use the same color family across the whole figure.

  lambda_values <- sort(unique(as.numeric(lambda_values)))
  lambda_values <- lambda_values[is.finite(lambda_values) & lambda_values >= 0]

  if (length(lambda_values) == 0L) {
    stop("lambda_values must contain at least one non-negative value.")
  }

  panel_list <- list()
  panel_counter <- 1L

  for (side in c("Superlevel analysis", "Sublevel analysis")) {
    y <- if (side == "Superlevel analysis") z else -z

    for (lambda in lambda_values) {
      cc <- connected_components_at_threshold(y, edges, lambda)

      tmp <- areal
      tmp$analysis_side <- side
      tmp$lambda <- lambda
      tmp$lambda_label <- sprintf("λ = %.2f", lambda)
      tmp$component_id <- cc$component
      tmp$n_active_panel <- cc$n_active
      tmp$beta0_panel <- cc$beta0
      tmp$panel_order <- panel_counter

      panel_list[[panel_counter]] <- tmp
      panel_counter <- panel_counter + 1L
    }
  }

  plot_dat <- do.call(rbind, panel_list)

  side_levels <- c("Superlevel analysis", "Sublevel analysis")
  plot_dat$analysis_side <- factor(plot_dat$analysis_side, levels = side_levels)

  lambda_levels <- sprintf("λ = %.2f", lambda_values)
  plot_dat$lambda_label <- factor(plot_dat$lambda_label, levels = lambda_levels)

  max_components <- suppressWarnings(max(plot_dat$component_id, na.rm = TRUE))
  if (!is.finite(max_components) || max_components <= 0L) {
    max_components <- 1L
  }

  component_palette <- grDevices::hcl.colors(max_components, palette = "Dark 3")
  names(component_palette) <- as.character(seq_len(max_components))

  plot_dat$component_factor <- factor(
    plot_dat$component_id,
    levels = seq_len(max_components)
  )

  p <- ggplot(plot_dat) +
    geom_sf(
      aes(fill = component_factor),
      color = map_border_color,
      linewidth = map_border_linewidth
    ) +
    facet_grid(
      rows = vars(analysis_side),
      cols = vars(lambda_label),
      switch = "y"
    ) +
    coord_sf_common() +
    scale_fill_manual(
      values = component_palette,
      na.value = "grey94",
      guide = "none",
      drop = FALSE
    ) +
    labs(
      title = wave_label,
      subtitle = NULL
    ) +
    paper_map_theme() +
    theme(
      legend.position = "none",
      strip.text.x = element_text(face = "bold", size = paper_base_size - 3),
      strip.text.y.left = element_text(face = "bold", size = paper_base_size - 3, angle = 90),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey96", color = NA),
      plot.title = element_text(face = "bold", size = paper_base_size + 1, hjust = 0.5)
    )

  list(
    plot = p,
    data = plot_dat
  )
}



############################################################
# 10. Run one wave analysis
############################################################

run_wave_topological_analysis <- function(wave_row, seed_offset = 0) {
  wave_id <- wave_row$wave_id
  wave_label <- wave_row$wave_label
  wave_start <- wave_row$start_date
  wave_end <- wave_row$end_date

  wave_dir <- file.path(output_dir, wave_id)
  dir.create(wave_dir, showWarnings = FALSE, recursive = TRUE)


  incidence_k <- wave_incidence %>%
    filter(wave_id == !!wave_id)

  areal <- areal_base %>%
    left_join(
      incidence_k,
      by = c("NUTS_ID" = "NUTS_join")
    )

  areal <- areal[
    is.finite(areal$incidence_per_100k) &
      areal$incidence_per_100k >= 0,
  ]

  rownames(areal) <- NULL

  # Rebuild the graph after dropping any missing incidence regions.
  # Usually this is identical to the base graph, but it is safer.
  nb_k <- spdep::poly2nb(
    areal,
    queen = queen_contiguity,
    snap = sqrt(.Machine$double.eps)
  )
  edges_k <- nb_to_edges(nb_k)

  x_raw <- as.numeric(areal$incidence_per_100k)
  z <- as.numeric(scale(x_raw))

  areal$scalar_raw <- x_raw
  areal$scalar_z <- z

  cat("\n====================================================\n")
  cat("Running wave:", wave_label, "\n")
  cat("Period:", as.character(wave_start), "to", as.character(wave_end), "\n")
  cat("Regions:", nrow(areal), "; edges:", nrow(edges_k), "\n")
  cat("Incidence summary:\n")
  print(summary(areal$scalar_raw))

  # Package analysis. Passing the raw incidence values with reference="mean"
  # and scale="sd" reproduces the within-wave standardization used in the paper.
  # The global analysis uses the same 80-point threshold grid as the final
  # case-study code; the local quantities themselves are integrated exactly.

  cat("\nGlobal test: superlevel graphs for high incidence...\n")
  global_pos <- spatconnect::global_connectivity(
    x = x_raw,
    graph = edges_k,
    reference = "mean",
    scale = "sd",
    side = "superlevel",
    n_perm = n_perm_global,
    seed = 1001 + seed_offset,
    n_grid = n_grid
  )

  cat("\nGlobal test: sublevel graphs for low incidence...\n")
  global_neg <- spatconnect::global_connectivity(
    x = x_raw,
    graph = edges_k,
    reference = "mean",
    scale = "sd",
    side = "sublevel",
    n_perm = n_perm_global,
    seed = 1002 + seed_offset,
    n_grid = n_grid
  )

  significance_for_package <- if (local_significance_rule == "fdr") {
    "adjusted"
  } else {
    "unadjusted"
  }

  cat("\nConditional local test for superlevel merging M_i...\n")
  local_pos <- spatconnect::local_connectivity(
    x = x_raw,
    graph = edges_k,
    reference = "mean",
    scale = "sd",
    side = "superlevel",
    n_perm = n_perm_local_conditional,
    seed = 3001 + seed_offset,
    adjust = "BH",
    alpha = local_alpha,
    activation_tail = activation_tail_eta,
    significance = significance_for_package,
    classify = TRUE,
    progress = TRUE
  )

  cat("\nConditional local test for sublevel merging M_i...\n")
  local_neg <- spatconnect::local_connectivity(
    x = x_raw,
    graph = edges_k,
    reference = "mean",
    scale = "sd",
    side = "sublevel",
    n_perm = n_perm_local_conditional,
    seed = 3002 + seed_offset,
    adjust = "BH",
    alpha = local_alpha,
    activation_tail = activation_tail_eta,
    significance = significance_for_package,
    classify = TRUE,
    progress = TRUE
  )

  activation_cutoff_pos <- local_pos$activation_cutoff
  activation_cutoff_neg <- local_neg$activation_cutoff

  cat("\nEmpirical activation cut-offs for typology:\n")
  cat("  superlevel A_i cut-off:", round(activation_cutoff_pos, 4), "\n")
  cat("  sublevel A_i cut-off:", round(activation_cutoff_neg, 4), "\n")

  pos_res <- local_pos$results
  neg_res <- local_neg$results

  areal$A_pos_index <- pos_res$A
  areal$M_pos_index <- pos_res$M
  areal$L_pos_index <- pos_res$L
  areal$A_neg_index <- neg_res$A
  areal$M_neg_index <- neg_res$M
  areal$L_neg_index <- neg_res$L

  # The package reports the upper-tail conditional p-value for M_i. Because
  # A_i is fixed for the focal area, this is exactly the lower-tail p-value for
  # L_i used in the paper. Existing column names are retained so downstream
  # tables and figure code remain directly comparable with earlier runs.
  areal$L_pos_p_less_cond <- pos_res$p_value
  areal$L_pos_p_less_cond_fdr <- pos_res$p_adjusted
  areal$L_neg_p_less_cond <- neg_res$p_value
  areal$L_neg_p_less_cond_fdr <- neg_res$p_adjusted

  areal$M_pos_excess_cond <- pos_res$M_excess
  areal$M_pos_z_cond <- pos_res$M_z_score
  areal$M_pos_null_mean_cond <- pos_res$M_null_mean
  areal$M_pos_null_sd_cond <- pos_res$M_null_sd
  areal$M_neg_excess_cond <- neg_res$M_excess
  areal$M_neg_z_cond <- neg_res$M_z_score
  areal$M_neg_null_mean_cond <- neg_res$M_null_mean
  areal$M_neg_null_sd_cond <- neg_res$M_null_sd

  areal$L_pos_excess_cond <- -areal$M_pos_excess_cond
  areal$L_pos_z_cond <- -areal$M_pos_z_cond
  areal$L_pos_null_mean_cond <- areal$A_pos_index - areal$M_pos_null_mean_cond
  areal$L_pos_null_sd_cond <- areal$M_pos_null_sd_cond
  areal$L_neg_excess_cond <- -areal$M_neg_excess_cond
  areal$L_neg_z_cond <- -areal$M_neg_z_cond
  areal$L_neg_null_mean_cond <- areal$A_neg_index - areal$M_neg_null_mean_cond
  areal$L_neg_null_sd_cond <- areal$M_neg_null_sd_cond

  areal$L_pos_relevant <- if (local_significance_rule == "fdr") {
    areal$L_pos_p_less_cond_fdr < local_alpha
  } else {
    areal$L_pos_p_less_cond < local_alpha
  }
  areal$L_neg_relevant <- if (local_significance_rule == "fdr") {
    areal$L_neg_p_less_cond_fdr < local_alpha
  } else {
    areal$L_neg_p_less_cond < local_alpha
  }
  areal$L_pos_relevant[is.na(areal$L_pos_relevant)] <- FALSE
  areal$L_neg_relevant[is.na(areal$L_neg_relevant)] <- FALSE

  areal$AM_class_pos <- pos_res$profile
  areal$AM_class_neg <- neg_res$profile

  # LISA benchmark
  cat("\nLISA benchmark...\n")

  lw_k <- spdep::nb2listw(
    nb_k,
    style = "W",
    zero.policy = TRUE
  )

  set.seed(6001 + seed_offset)

  lisa <- spdep::localmoran_perm(
    x = areal$scalar_z,
    listw = lw_k,
    nsim = n_perm_local_moran,
    zero.policy = TRUE,
    alternative = "two.sided"
  )

  areal$local_moran_I <- lisa[, "Ii"]

  p_candidates <- c(
    "Pr(folded) Sim",
    "Pr(z != E(Ii)) Sim",
    "Pr(z != E(Ii))"
  )

  p_col <- p_candidates[p_candidates %in% colnames(lisa)][1]

  if (is.na(p_col)) {
    p_col <- grep("Pr", colnames(lisa), value = TRUE)[1]
  }

  alpha <- 0.05

  areal$local_moran_p <- lisa[, p_col]
  areal$local_moran_p_fdr <- p.adjust(areal$local_moran_p, method = "BH")

  areal$lag_z <- spdep::lag.listw(
    lw_k,
    areal$scalar_z,
    zero.policy = TRUE
  )

  areal$lisa_cluster <- "Not significant"

  areal$lisa_cluster[
    areal$scalar_z > 0 &
      areal$lag_z > 0 &
      areal$local_moran_p < alpha
  ] <- "High-High"

  areal$lisa_cluster[
    areal$scalar_z < 0 &
      areal$lag_z < 0 &
      areal$local_moran_p < alpha
  ] <- "Low-Low"

  areal$lisa_cluster[
    areal$scalar_z > 0 &
      areal$lag_z < 0 &
      areal$local_moran_p < alpha
  ] <- "High-Low"

  areal$lisa_cluster[
    areal$scalar_z < 0 &
      areal$lag_z > 0 &
      areal$local_moran_p < alpha
  ] <- "Low-High"

  areal$lisa_cluster <- factor(
    areal$lisa_cluster,
    levels = c(
      "High-High",
      "Low-Low",
      "High-Low",
      "Low-High",
      "Not significant"
    )
  )

  areal$wave_id <- wave_id
  areal$wave_label <- wave_label
  areal$wave_start <- wave_start
  areal$wave_end <- wave_end

  # Tables
  tab_AM_pos_lisa <- table(
    AM_class_pos = areal$AM_class_pos,
    LISA = areal$lisa_cluster
  )
  tab_AM_neg_lisa <- table(
    AM_class_neg = areal$AM_class_neg,
    LISA = areal$lisa_cluster
  )

  summary_table <- data.frame(
    wave_id = wave_id,
    wave_label = wave_label,
    wave_start = as.character(wave_start),
    wave_end = as.character(wave_end),
    n_areas = nrow(areal),
    n_edges = nrow(edges_k),
    n_isolates = sum(spdep::card(nb_k) == 0L),
    n_perm_global = n_perm_global,
    n_perm_local_conditional = n_perm_local_conditional,
    n_perm_local_moran = n_perm_local_moran,
    spatconnect_version = as.character(utils::packageVersion("spatconnect")),
    global_integration = global_pos$integration_method,
    local_allocation_rule = local_pos$allocation_rule,
    local_alpha = local_alpha,
    local_significance_rule = local_significance_rule,
    activation_tail_eta = activation_tail_eta,
    activation_cutoff_pos = activation_cutoff_pos,
    activation_cutoff_neg = activation_cutoff_neg,
    high_D_abs = global_pos$discrepancy,
    high_p_abs = global_pos$p_value,
    low_D_abs = global_neg$discrepancy,
    low_p_abs = global_neg$p_value,
    max_A_pos = max(areal$A_pos_index, na.rm = TRUE),
    max_M_pos = max(areal$M_pos_index, na.rm = TRUE),
    max_M_pos_excess_cond = max(areal$M_pos_excess_cond, na.rm = TRUE),
    max_A_neg = max(areal$A_neg_index, na.rm = TRUE),
    max_M_neg = max(areal$M_neg_index, na.rm = TRUE),
    max_M_neg_excess_cond = max(areal$M_neg_excess_cond, na.rm = TRUE),
    n_L_pos_cond_p_lt_alpha = sum(areal$L_pos_p_less_cond < local_alpha, na.rm = TRUE),
    n_L_pos_cond_fdr_lt_alpha = sum(areal$L_pos_p_less_cond_fdr < local_alpha, na.rm = TRUE),
    n_L_neg_cond_p_lt_alpha = sum(areal$L_neg_p_less_cond < local_alpha, na.rm = TRUE),
    n_L_neg_cond_fdr_lt_alpha = sum(areal$L_neg_p_less_cond_fdr < local_alpha, na.rm = TRUE),
    n_lisa_HH = sum(areal$lisa_cluster == "High-High", na.rm = TRUE),
    n_lisa_LL = sum(areal$lisa_cluster == "Low-Low", na.rm = TRUE),
    n_pos_core_connector = sum(areal$AM_class_pos == "Core connector", na.rm = TRUE),
    n_pos_bridge_connector = sum(areal$AM_class_pos == "Bridge connector", na.rm = TRUE),
    n_pos_peak_connector = sum(areal$AM_class_pos == "Peak connector", na.rm = TRUE),
    n_pos_not_significant = sum(areal$AM_class_pos == "Not significant", na.rm = TRUE),
    n_neg_core_connector = sum(areal$AM_class_neg == "Core connector", na.rm = TRUE),
    n_neg_bridge_connector = sum(areal$AM_class_neg == "Bridge connector", na.rm = TRUE),
    n_neg_trough_connector = sum(areal$AM_class_neg == "Trough connector", na.rm = TRUE),
    n_neg_not_significant = sum(areal$AM_class_neg == "Not significant", na.rm = TRUE)
  )

  print(summary_table)

  write.csv(
    summary_table,
    file.path(wave_dir, paste0(wave_id, "_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    as.data.frame.matrix(tab_AM_pos_lisa),
    file.path(wave_dir, paste0(wave_id, "_tab_AM_positive_vs_lisa.csv"))
  )

  write.csv(
    as.data.frame.matrix(tab_AM_neg_lisa),
    file.path(wave_dir, paste0(wave_id, "_tab_AM_negative_vs_lisa.csv"))
  )

  # Figures
  p_raw <- plot_sf_map(
    areal,
    var = "scalar_raw",
    title = "Cumulative COVID-19 incidence",
    subtitle = wave_label,
    fill_label = "Incidence per 100,000",
    scale_type = "sequential"
  )
  print(p_raw)
  save_plot(p_raw, "01_raw_cumulative_incidence_per100k.eps", wave_dir, single_map_width, single_map_height)

  p_z <- plot_sf_map(
    areal,
    var = "scalar_z",
    title = "Standardized incidence",
    subtitle = wave_label,
    fill_label = expression(z[i]),
    scale_type = "diverging"
  )
  print(p_z)
  save_plot(p_z, "02_standardized_incidence_z.eps", wave_dir, single_map_width, single_map_height)

  # Descriptive thresholded-graph figure: connected components for the same
  # lambda values in the superlevel and sublevel analyses.
  component_fig <- make_threshold_component_figure(
    areal = areal,
    z = z,
    edges = edges_k,
    lambda_values = component_figure_lambdas,
    wave_label = wave_label
  )
  print(component_fig$plot)
  save_plot(
    component_fig$plot,
    "02b_connected_components_thresholds.eps",
    wave_dir,
    width = 13,
    height = 8.4
  )

  write.csv(
    sf::st_drop_geometry(component_fig$data)[
      ,
      c(
        "NUTS_ID", "analysis_side", "lambda", "component_id",
        "n_active_panel", "beta0_panel"
      )
    ],
    file.path(wave_dir, paste0(wave_id, "_threshold_component_membership.csv")),
    row.names = FALSE
  )

  p_global_pos <- plot_global_beta0(
    global_pos,
    title = "Superlevel Betti-0 connectivity",
    subtitle = wave_label
  )
  print(p_global_pos)
  save_plot(p_global_pos, "03_global_curve_high_incidence.eps", wave_dir, 8.5, 6.2)

  p_global_neg <- plot_global_beta0(
    global_neg,
    title = "Sublevel Betti-0 connectivity",
    subtitle = wave_label
  )
  print(p_global_neg)
  save_plot(p_global_neg, "04_global_curve_low_incidence.eps", wave_dir, 8.5, 6.2)

  p_A_pos <- plot_sf_map(
    areal,
    var = "A_pos_index",
    title = "Superlevel activation",
    subtitle = wave_label,
    fill_label = expression(A[i]),
    scale_type = "sequential"
  )
  print(p_A_pos)
  save_plot(p_A_pos, "05_A_pos_activation.eps", wave_dir, single_map_width, single_map_height)

  p_M_pos <- plot_sf_map(
    areal,
    var = "M_pos_index",
    title = "Superlevel merging contribution",
    subtitle = wave_label,
    fill_label = expression(M[i]),
    scale_type = "sequential"
  )
  print(p_M_pos)
  save_plot(p_M_pos, "06_M_pos_merging.eps", wave_dir, single_map_width, single_map_height)

  p_M_pos_excess <- plot_sf_map(
    areal,
    var = "M_pos_excess_cond",
    title = "Conditional excess superlevel merging",
    subtitle = wave_label,
    fill_label = expression(M[i] - E[0](M[i])),
    scale_type = "diverging"
  )
  print(p_M_pos_excess)
  save_plot(p_M_pos_excess, "07_M_pos_excess_conditional.eps", wave_dir, single_map_width, single_map_height)

  p_pos_pvalues <- plot_sf_map(
    areal,
    var = "L_pos_p_less_cond",
    title = "Superlevel connector p-values",
    subtitle = wave_label,
    fill_label = expression(p[i]),
    scale_type = "pvalue"
  )
  print(p_pos_pvalues)
  save_plot(p_pos_pvalues, "08_L_pos_pvalues_conditional.eps", wave_dir, single_map_width, single_map_height)

  p_AM_pos <- ggplot(areal) +
    geom_sf(aes(fill = AM_class_pos), color = map_border_color, linewidth = map_border_linewidth, show.legend = TRUE, key_glyph = draw_key_rect) +
    coord_sf_common() +
    scale_fill_manual(
      values = activation_merging_colors_pos,
      breaks = names(activation_merging_colors_pos),
      limits = names(activation_merging_colors_pos),
      labels = activation_merging_labels_pos,
      drop = FALSE
    ) +
    force_discrete_fill_legend(activation_merging_colors_pos) +
    labs(
      title = "Superlevel connectivity profiles",
      subtitle = wave_label,
      fill = "Connectivity role"
    ) +
    paper_map_theme()
  print(p_AM_pos)
  save_plot(p_AM_pos, "09_significant_upper_local_profiles.eps", wave_dir, single_map_width, single_map_height)

  am_scatter_pos_df <- st_drop_geometry(areal) %>%
    mutate(
      p_cond_plot = pmax(L_pos_p_less_cond, 1 / (n_perm_local_conditional + 1)),
      p_strength = -log10(p_cond_plot),
      scatter_role = dplyr::case_when(
        AM_class_pos == "Core connector" ~ "Core connector",
        AM_class_pos == "Bridge connector" ~ "Bridge connector",
        AM_class_pos == "Peak connector" ~ "Peak connector",
        TRUE ~ "Not significant"
      ),
      scatter_role = factor(scatter_role, levels = scatter_role_levels_pos)
    )

  p_AM_scatter_pos <- ggplot(
    am_scatter_pos_df,
    aes(x = A_pos_index, y = M_pos_index, color = scatter_role, size = p_strength)
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dotted",
      color = "grey35",
      linewidth = 0.7
    ) +
    geom_vline(
      xintercept = activation_cutoff_pos,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.75
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey70",
      linewidth = 0.35
    ) +
    geom_point(alpha = 0.85) +
    scale_color_manual(
      values = scatter_role_colors_pos,
      breaks = scatter_role_levels_pos,
      limits = scatter_role_levels_pos,
      drop = FALSE
    ) +
    scale_size_continuous(
      range = c(2.2, 7.0),
      limits = c(0, -log10(1 / (n_perm_local_conditional + 1))),
      breaks = -log10(c(0.20, 0.05, 0.01)),
      labels = c("0.20", "0.05", "0.01"),
      name = expression("Conditional"~p[i])
    ) +
    stable_scatter_guides() +
    labs(
      title = "Superlevel activation–merging diagram",
      subtitle = wave_label,
      x = expression(A[i]),
      y = expression(M[i]),
      color = "Connectivity role"
    ) +
    paper_theme() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top"
    )
  print(p_AM_scatter_pos)
  save_plot(p_AM_scatter_pos, "10_AM_scatter_high_incidence_conditional.eps", wave_dir, 8.8, 6.5)

  p_A_neg <- plot_sf_map(
    areal,
    var = "A_neg_index",
    title = "Sublevel activation",
    subtitle = wave_label,
    fill_label = expression(A[i]),
    scale_type = "sequential"
  )
  print(p_A_neg)
  save_plot(p_A_neg, "11_A_neg_activation.eps", wave_dir, single_map_width, single_map_height)

  p_M_neg <- plot_sf_map(
    areal,
    var = "M_neg_index",
    title = "Sublevel merging contribution",
    subtitle = wave_label,
    fill_label = expression(M[i]),
    scale_type = "sequential"
  )
  print(p_M_neg)
  save_plot(p_M_neg, "12_M_neg_merging.eps", wave_dir, single_map_width, single_map_height)

  p_M_neg_excess <- plot_sf_map(
    areal,
    var = "M_neg_excess_cond",
    title = "Conditional excess sublevel merging",
    subtitle = wave_label,
    fill_label = expression(M[i] - E[0](M[i])),
    scale_type = "diverging"
  )
  print(p_M_neg_excess)
  save_plot(p_M_neg_excess, "13_M_neg_excess_conditional.eps", wave_dir, single_map_width, single_map_height)

  p_neg_pvalues <- plot_sf_map(
    areal,
    var = "L_neg_p_less_cond",
    title = "Sublevel connector p-values",
    subtitle = wave_label,
    fill_label = expression(p[i]),
    scale_type = "pvalue"
  )
  print(p_neg_pvalues)
  save_plot(p_neg_pvalues, "14_L_neg_pvalues_conditional.eps", wave_dir, single_map_width, single_map_height)

  p_AM_neg <- ggplot(areal) +
    geom_sf(aes(fill = AM_class_neg), color = map_border_color, linewidth = map_border_linewidth, show.legend = TRUE, key_glyph = draw_key_rect) +
    coord_sf_common() +
    scale_fill_manual(
      values = activation_merging_colors_neg,
      breaks = names(activation_merging_colors_neg),
      limits = names(activation_merging_colors_neg),
      labels = activation_merging_labels_neg,
      drop = FALSE
    ) +
    force_discrete_fill_legend(activation_merging_colors_neg) +
    labs(
      title = "Sublevel connectivity profiles",
      subtitle = wave_label,
      fill = "Connectivity role"
    ) +
    paper_map_theme()
  print(p_AM_neg)
  save_plot(p_AM_neg, "15_significant_lower_local_profiles.eps", wave_dir, single_map_width, single_map_height)

  am_scatter_neg_df <- st_drop_geometry(areal) %>%
    mutate(
      p_cond_plot = pmax(L_neg_p_less_cond, 1 / (n_perm_local_conditional + 1)),
      p_strength = -log10(p_cond_plot),
      scatter_role = dplyr::case_when(
        AM_class_neg == "Core connector" ~ "Core connector",
        AM_class_neg == "Bridge connector" ~ "Bridge connector",
        AM_class_neg == "Trough connector" ~ "Trough connector",
        TRUE ~ "Not significant"
      ),
      scatter_role = factor(scatter_role, levels = scatter_role_levels_neg)
    )

  p_AM_scatter_neg <- ggplot(
    am_scatter_neg_df,
    aes(x = A_neg_index, y = M_neg_index, color = scatter_role, size = p_strength)
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dotted",
      color = "grey35",
      linewidth = 0.7
    ) +
    geom_vline(
      xintercept = activation_cutoff_neg,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.75
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey70",
      linewidth = 0.35
    ) +
    geom_point(alpha = 0.85) +
    scale_color_manual(
      values = scatter_role_colors_neg,
      breaks = scatter_role_levels_neg,
      limits = scatter_role_levels_neg,
      drop = FALSE
    ) +
    scale_size_continuous(
      range = c(2.2, 7.0),
      limits = c(0, -log10(1 / (n_perm_local_conditional + 1))),
      breaks = -log10(c(0.20, 0.05, 0.01)),
      labels = c("0.20", "0.05", "0.01"),
      name = expression("Conditional"~p[i])
    ) +
    stable_scatter_guides() +
    labs(
      title = "Sublevel activation–merging diagram",
      subtitle = wave_label,
      x = expression(A[i]),
      y = expression(M[i]),
      color = "Connectivity role"
    ) +
    paper_theme() +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.box.just = "top"
    )
  print(p_AM_scatter_neg)
  save_plot(p_AM_scatter_neg, "15b_AM_scatter_low_incidence_conditional.eps", wave_dir, 8.8, 6.5)

  p_lisa <- ggplot(areal) +
    geom_sf(aes(fill = lisa_cluster), color = map_border_color, linewidth = map_border_linewidth, show.legend = TRUE, key_glyph = draw_key_rect) +
    coord_sf_common() +
    scale_fill_manual(values = lisa_colors, breaks = names(lisa_colors), limits = names(lisa_colors), drop = FALSE) +
    force_discrete_fill_legend(lisa_colors) +
    labs(
      title = "LISA cluster map",
      subtitle = wave_label,
      fill = "LISA cluster"
    ) +
    paper_map_theme()
  print(p_lisa)
  save_plot(p_lisa, "16_lisa_cluster_map.eps", wave_dir, single_map_width, single_map_height)

  areal$AM_pos_connector <- areal$AM_class_pos %in% c(
    "Core connector",
    "Bridge connector",
    "Peak connector"
  )

  p_lisa_AM_pos <- ggplot(areal) +
    geom_sf(aes(fill = lisa_cluster), color = map_border_color, linewidth = map_border_linewidth, show.legend = TRUE, key_glyph = draw_key_rect) +
    geom_sf(
      data = areal[areal$AM_pos_connector, ],
      fill = NA,
      color = connector_outline_color,
      linewidth = connector_outline_linewidth
    ) +
    coord_sf_common() +
    scale_fill_manual(values = lisa_colors, breaks = names(lisa_colors), limits = names(lisa_colors), drop = FALSE) +
    force_discrete_fill_legend(lisa_colors) +
    labs(
      title = "LISA with superlevel connectivity profiles",
      subtitle = wave_label,
      fill = "LISA cluster"
    ) +
    paper_map_theme()
  print(p_lisa_AM_pos)
  save_plot(p_lisa_AM_pos, "17_lisa_with_significant_upper_profiles.eps", wave_dir, single_map_width, single_map_height)

  # Save enriched geometry.
  out_gpkg <- file.path(wave_dir, paste0(wave_id, "_covid_topological_results.gpkg"))
  sf::st_write(areal, out_gpkg, delete_dsn = TRUE, quiet = TRUE)
  cat("Saved enriched spatial object:", out_gpkg, "\n")

  list(
    wave_id = wave_id,
    wave_label = wave_label,
    areal = areal,
    summary = summary_table,
    global_pos = global_pos,
    global_neg = global_neg,
    local_pos = local_pos,
    local_neg = local_neg,
    local_L_pos_cond = local_pos,
    local_L_neg_cond = local_neg,
    component_figure = component_fig,
    tab_AM_pos_lisa = tab_AM_pos_lisa,
    tab_AM_neg_lisa = tab_AM_neg_lisa
  )
}

############################################################
# 11. Run all wave analyses and save combined outputs
############################################################

results <- list()

for (k in seq_len(nrow(waves))) {
  results[[waves$wave_id[k]]] <- run_wave_topological_analysis(
    wave_row = waves[k, ],
    seed_offset = 10000 * k
  )
}

combined_summary <- bind_rows(lapply(results, function(x) x$summary))

write.csv(
  combined_summary,
  file.path(output_dir, "combined_two_wave_summary.csv"),
  row.names = FALSE
)

print(combined_summary)

# Combined map dataset.
combined_areal <- bind_rows(lapply(results, function(x) x$areal))

sf::st_write(
  combined_areal,
  file.path(output_dir, "combined_two_wave_topological_results.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

# Faceted comparison maps.
p_compare_raw <- ggplot(combined_areal) +
  geom_sf(aes(fill = scalar_raw), color = map_border_color, linewidth = map_border_linewidth) +
  facet_wrap(~ wave_label) +
  coord_sf_common() +
  scale_fill_gradientn(colors = ylorrd, na.value = "grey90") +
  labs(
    title = "Italy COVID-19 cumulative incidence by wave",
    subtitle = NULL,
    fill = "Incidence"
  ) +
  paper_map_theme()

print(p_compare_raw)
save_plot(p_compare_raw, "combined_01_raw_incidence_by_wave.eps", output_dir, 12, 7)

p_compare_M <- ggplot(combined_areal) +
  geom_sf(aes(fill = M_pos_index), color = map_border_color, linewidth = map_border_linewidth) +
  facet_wrap(~ wave_label) +
  coord_sf_common() +
  scale_fill_gradientn(colors = ylorrd, na.value = "grey90") +
  labs(
    title = "Superlevel merging contribution by wave",
    subtitle = NULL,
    fill = expression(M[i])
  ) +
  paper_map_theme()

print(p_compare_M)
save_plot(p_compare_M, "combined_02_M_pos_by_wave.eps", output_dir, 12, 7)

p_compare_M_excess <- ggplot(combined_areal) +
  geom_sf(aes(fill = M_pos_excess_cond), color = map_border_color, linewidth = map_border_linewidth) +
  facet_wrap(~ wave_label) +
  coord_sf_common() +
  scale_fill_gradient2(
    low = rdbu[11],
    mid = rdbu[6],
    high = rdbu[1],
    midpoint = 0,
    na.value = "grey90"
  ) +
  labs(
    title = "Conditional excess superlevel merging by wave",
    subtitle = NULL,
    fill = expression(M[i] - E[0](M[i]))
  ) +
  paper_map_theme()

print(p_compare_M_excess)
save_plot(p_compare_M_excess, "combined_03_M_pos_excess_conditional_by_wave.eps", output_dir, 12, 7)

p_compare_AM <- ggplot(combined_areal) +
  geom_sf(aes(fill = AM_class_pos), color = map_border_color, linewidth = map_border_linewidth, show.legend = TRUE, key_glyph = draw_key_rect) +
  facet_wrap(~ wave_label) +
  coord_sf_common() +
  scale_fill_manual(
      values = activation_merging_colors_pos,
      breaks = names(activation_merging_colors_pos),
      limits = names(activation_merging_colors_pos),
      labels = activation_merging_labels_pos,
      drop = FALSE
    ) +
    force_discrete_fill_legend(activation_merging_colors_pos) +
  labs(
    title = "Superlevel connectivity profiles",
    subtitle = NULL,
    fill = "Connectivity role"
  ) +
  paper_map_theme()

print(p_compare_AM)
save_plot(p_compare_AM, "combined_04_significant_upper_local_profiles_by_wave.eps", output_dir, 12, 7)

cat("\nFinished Italy two-wave COVID analysis.\n")
cat("Main output folder:\n")
cat(output_dir, "\n")
