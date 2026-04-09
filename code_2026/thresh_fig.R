
########################################################
#FUNCTION DEFINITIONS
########################################################
library(magrittr)

calc_det_prob <- function(data, selected_taxon_level = "scientificName", selected_taxon_id = "All", pool_primers = FALSE) {

  oop <- options("dplyr.summarise.inform")
  options(dplyr.summarise.inform = FALSE)
  # reset option on exit
  on.exit(options(dplyr.summarise.inform = oop))



  data <- data %>%
    dplyr::group_by(
      protocol_ID,
      primer,
      month,
      year,
      .data[[selected_taxon_level]],
      samp_name
    ) %>%
    dplyr::summarise(
      detected = as.integer(any(detected == 1)),
      .groups = "drop"
    )

  bad <- !data$detected %in% c(0, 1)

  if (any(bad)) {
    stop(
      "calc_det_prob(): detected must be 0/1 after collapsing.\n",
      "Found ", sum(bad), " invalid values.\n",
      "Unique invalid values: ",
      paste(unique(data$detected[bad]), collapse = ", ")
    )
  }

  if (pool_primers) {
    data %<>%
      dplyr::mutate(.,
                    id = paste0(protocol_ID, ";", .data[[selected_taxon_level]]),
                    id.yr = paste0(protocol_ID, ";", .data[[selected_taxon_level]], ";ALLPRIMERS;", year)
      )
  } else {
    data %<>%
      dplyr::mutate(.,
                    id = paste0(protocol_ID, ";", .data[[selected_taxon_level]], ";", primer),
                    id.yr = paste0(protocol_ID, ";", .data[[selected_taxon_level]], ";", primer, ";", year)
      )
  }

  # Create a variable so detection probability is calculated separately for each
  # protocol ID, version, selected_taxon_level, and primer

  # create new list variables to store outputs
  lnd <- length(unique(data$id))
  newP <- vector("list", lnd)
  names(newP) <- unique(data$id)
  SUM <- COM <- comps <- vector("list", length(unique(data$id)))
  names(COM) <- unique(data$id)

  # calculate detection probabilities - year aggregated
  for (occurrence in unique(data$id)) {
    SDF <- data %>%
      dplyr::filter(id == occurrence) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(
        n = dplyr::n(),
        nd = sum(detected),
        p = nd/n,
        s = sqrt(p * (1 - p)/n)
      ) %>%
      as.data.frame()

    if (any(SDF$n > 1)) {
      newP[[occurrence]] <- SDF
    }
  }
  newP_agg <- newP[lengths(newP) != 0]


  # calculate monthly detection probability for each year
  lny <- length(unique(data$id.yr))
  newP <- vector("list", lny)
  names(newP) <- unique(data$id.yr)
  SUM <- COM <- comps <- vector("list", length(unique(data$id.yr)))
  names(COM) <- unique(data$id.yr)

  for (occurrence in unique(data$id.yr)) {
    SDF <- data %>%
      dplyr::filter(id.yr == occurrence) %>%
      dplyr::group_by(month) %>%
      dplyr::summarise(
        n = dplyr::n(),
        nd = sum(detected),
        p = nd/n,
        s = sqrt(p * (1 - p)/n)
      ) %>%
      as.data.frame()

    if (any(SDF$n > 1)) {
      newP[[occurrence]] <- SDF
    }
  }
  newP <- newP[lengths(newP) != 0]
  newP_yr <- Map(cbind, newP, id.yr = names(newP))

  newP_yr <- lapply(newP_yr, function(x) {
    dplyr::mutate(x,
                  protocol_ID = stringr::word(x$id.yr, 1, sep = stringr::fixed(";")),
                  year = stringr::word(x$id.yr, -1, sep = stringr::fixed(";")),
                  yr.mo = paste0(year, ";", month),
                  id.yr = NULL
    )
  }) # to obtain variation among years

  list(newP_agg = newP_agg, newP_yr = newP_yr)
}

drop_all_zero_taxa <- function(df, taxon_col) {
  df %>%
    group_by(
      station,
      protocol_ID,
      .data[[taxon_col]]
    ) %>%
    filter(any(detected == 1)) %>%
    ungroup()
}




filter_nondetections_specifics <- function(df_subset, distance = 500) {
  sf::sf_use_s2(TRUE)

  if (nrow(df_subset) == 0) return(df_subset)

  # 1) Coordinate-level detection flag (keep ALL rows after join)
  coord_detection_summary <- df_subset %>%
    group_by(Longitude, Latitude) %>%
    summarise(
      coord_has_detection = any(detected == 1),
      .groups = "drop"
    )

  df_subset <- df_subset %>%
    left_join(coord_detection_summary,
              by = c("Longitude", "Latitude"))

  # 2) Unique coordinates with coord_id
  coords <- df_subset %>%
    distinct(Longitude, Latitude, coord_has_detection) %>%
    mutate(coord_id = row_number())

  df_subset <- df_subset %>%
    left_join(coords, by = c("Longitude", "Latitude", "coord_has_detection"))

  # If only one coordinate, the only way to be "within distance of a detection"
  # is if that coordinate itself has a detection.
  if (nrow(coords) == 1) {
    return(df_subset %>% filter(coord_has_detection))
  }

  # 3) sf points in lon/lat (global-safe)
  coords_sf <- sf::st_as_sf(
    coords,
    coords = c("Longitude", "Latitude"),
    crs = 4326
  )

  # 4) Neighbor list within distance (meters)
  nbrs <- sf::st_is_within_distance(coords_sf, coords_sf, dist = distance)

  # 5) Propagate coord_has_detection through neighbor lists
  coord_detected_within_dist <- vapply(
    seq_along(nbrs),
    function(i) any(coords$coord_has_detection[nbrs[[i]]]),
    logical(1)
  )

  # 6) Attach to all rows + final filter rule:
  # keep rows if:
  #   - taxon was ever detected at that coordinate (coord_has_detection)
  #   OR
  #   - taxon detected within distance of that coordinate (detected_within_dist)
  df_subset %>%
    mutate(detected_within_dist = coord_detected_within_dist[coord_id]) %>%
    filter(coord_has_detection | detected_within_dist)
}


filter_nondetections_all <- function(df,
                                     distance = 500,
                                     selected_taxon_level = "scientificName",
                                     selected_taxon_id = "All") {

  stopifnot(selected_taxon_level %in% names(df))
  stopifnot(all(c("protocol_ID", "Longitude", "Latitude", "detected") %in% names(df)))

  # # Optional filter to one taxon
  # if (!identical(selected_taxon_id, "All")) {
  #   df <- df %>% filter(.data[[selected_taxon_level]] == selected_taxon_id)
  # }
  print(nrow(df))
  if (nrow(df) == 0) return(df)

  df %>%
    group_by(protocol_ID, primer, .data[[selected_taxon_level]]) %>%
    group_modify(~ filter_nondetections_specifics(.x, distance = distance)) %>%
    ungroup()
}

scale_newprob <- function(data, newprob, selected_taxon_level = "scientificName") {
  CPscaled <- lapply(newprob, function(species_list) {

    # 1. Extract the maximum p across all elements for this species
    all_p <- unlist(lapply(species_list, function(df) df$p))
    max_p <- max(all_p, na.rm = TRUE)

    # 2. Apply minimal change: scale by max_p instead of scale_prop
    lapply(species_list, function(y) {
      data.frame(y) |>
        dplyr::mutate(
          scaleP = dplyr::case_when(
            p == 1 ~ 1,
            p == 0 ~ 0,
            TRUE   ~ p / max_p
          )
        )
    })
  })
  DFmo <- lapply(CPscaled$newP_agg, function(x) {
    out <- data.frame(
      month = 1:12,
      detect = NA_integer_,
      nondetect = NA_integer_,
      scaleP = NA_real_
    )
    out$detect[x$month] <- x$nd
    out$nondetect[x$month] <- x$n - x$nd
    out$scaleP[x$month] <- x$scaleP
    out
  }) |>
    do.call(what = rbind) |>
    dplyr::mutate(
      id = rep(names(CPscaled$newP_agg), each = 12)
    ) |>
    dplyr::select(id, month, detect, nondetect, scaleP) |>
    dplyr::tibble()
  row.names(DFmo) <- NULL

  DFmo[c("protocol_ID", selected_taxon_level, "primer")] <- stringr::str_split_fixed(DFmo$id, ";", 3)

  taxa_columns <- c("All", "kingdom", "phylum", "class", "order", "family", "genus", "scientificName")
  cols_to_keep <- taxa_columns[1:which(taxa_columns == selected_taxon_level)]

  DFmo <- DFmo |>
    dplyr::left_join(
      unique(data[, cols_to_keep]),
      by = selected_taxon_level,
      multiple = "first"
    )

  # Interpolate missing months
  DFmo$det_int <- NA
  DFmo$nd_int <- NA
  DFmo$fill <- NA # add column

  for (taxon in unique(DFmo$id)) {
    DF1 <- DFmo[DFmo$id == taxon, ]

    # then add code for interpolation that starts with DF2 = .....
    # add dataframe above and below to help will fills for jan and dec. Needed to have 4 copyies because of max function used below
    DF2 <- rbind(
      cbind(DF1, data.frame(G = 1)),
      cbind(DF1, data.frame(G = 2)),
      cbind(DF1, data.frame(G = 3)),
      cbind(DF1, data.frame(G = 4))
    )

    DF2$det_int <- DF2$detect
    DF2$nd_int <- DF2$nondetect
    DF2$fill <- DF2$scaleP

    # which months are NA and define groups with sequential NAs
    month_na_id <- which(is.na(DF2$detect))
    nagroups <- cumsum(c(1, abs(month_na_id[-length(month_na_id)] - month_na_id[-1]) > 1))

    # identify which NA groups are in G = 2 or 3 (ignore 1 and 2)
    nagroupsG <- list()
    for (i in unique(nagroups)) {
      nagroupsG[[i]] <- max(DF2$G[month_na_id[which(nagroups == i)]])
    }
    nagroupsGv <- unlist(nagroupsG)
    nagroupsf <- which(nagroupsGv %in% 2:3)

    # loop over final NA groups and fill in using average
    for (i in unique(nagroupsf)) {
      DF2$det_int[month_na_id[which(nagroups == i)]] <- (DF2$detect[min(
        month_na_id[which(nagroups == i)]) - 1] + DF2$detect[max(
          month_na_id[which(nagroups == i)]) + 1])/2

      DF2$nd_int[month_na_id[which(nagroups == i)]] <- (DF2$nondetect[min(
        month_na_id[which(nagroups == i)]) - 1] + DF2$nondetect[max(
          month_na_id[which(nagroups == i)]) + 1])/2

      DF2$fill[month_na_id[which(nagroups == i)]] <- (DF2$scaleP[min(
        month_na_id[which(nagroups == i)]) - 1] + DF2$scaleP[max(
          month_na_id[which(nagroups == i)]) + 1])/2
    }
    # then put values from DF3 back into DF$sp.pr. This assumes that the months are all in the correct order (jan to dec) in DF3 and test_interp
    # DF3 is final DF with fills
    DF3 <- DF2[DF2$G == 2, ]

    DFmo$det_int[DFmo$id == taxon] <- DF3$det_int
    DFmo$nd_int[DFmo$id == taxon] <- DF3$nd_int
    DFmo$fill[DFmo$id == taxon] <- DF3$fill
  }

  # Pscaled_month <- DFmo %>%
  #  dplyr::ungroup()

  # scale and interpolate each month separately
  DFyr <- lapply(CPscaled$newP_yr, function(x) {
    out <- data.frame(
      month = 1:12,
      detect = NA_integer_,
      nondetect = NA_integer_,
      scaleP = NA_real_
    )
    out$detect[x$month] <- x$nd
    out$nondetect[x$month] <- x$n - x$nd
    out$scaleP[x$month] <- x$scaleP
    out
  }) |>
    do.call(what = rbind) |>
    dplyr::mutate(
      id = rep(names(CPscaled$newP_yr), each = 12)
    ) |>
    dplyr::select(id, month, detect, nondetect, scaleP) |>
    dplyr::tibble()
  row.names(DFyr) <- NULL

  DFyr[c("protocol_ID", selected_taxon_level, "primer", "year")] <- stringr::str_split_fixed(DFyr$id, ";", 4)

  taxa_columns <- c("All", "kingdom", "phylum", "class", "order", "family", "genus", "scientificName")
  cols_to_keep <- taxa_columns[1:which(taxa_columns == selected_taxon_level)]

  DFyr <- DFyr |>
    dplyr::left_join(
      unique(data[, cols_to_keep]),
      by = selected_taxon_level,
      multiple = "first"
    )

  # Interpolate missing months
  DFyr$fill <- NA # add column

  for (taxon in unique(DFyr$id)) {
    DF1 <- DFyr[DFyr$id == taxon, ]

    # then add code for interpolation that starts with DF2 = .....
    # add dataframe above and below to help will fills for jan and dec. Needed to have 4 copies because of max function used below
    DF2 <- rbind(
      cbind(DF1, data.frame(G = 1)),
      cbind(DF1, data.frame(G = 2)),
      cbind(DF1, data.frame(G = 3)),
      cbind(DF1, data.frame(G = 4))
    )

    DF2$fill <- DF2$scaleP

    # which months are NA and define groups with sequential NAs
    month_na_id <- which(is.na(DF2$scaleP))
    nagroups <- cumsum(c(1, abs(month_na_id[-length(month_na_id)] - month_na_id[-1]) > 1))

    # identify which NA groups are in G = 2 or 3 (ignore 1 and 2)
    nagroupsG <- list()
    for (i in unique(nagroups)) {
      nagroupsG[[i]] <- max(DF2$G[month_na_id[which(nagroups == i)]])
    }
    nagroupsGv <- unlist(nagroupsG)
    nagroupsf <- which(nagroupsGv %in% 2:3)

    # loop over final NA groups and fill in using average
    for (i in unique(nagroupsf)) {
      DF2$fill[month_na_id[which(nagroups == i)]] <- (DF2$scaleP[min(month_na_id[which(nagroups == i)]) - 1] + DF2$scaleP[max(month_na_id[which(nagroups == i)]) + 1]) / 2
    }
    # then put values from DF3 back into DF$id.sp.pr.yr. This assumes that the months are all in the correct order (jan to dec) in DF3 and test_interp
    # DF3 is final DF with fills
    DF3 <- DF2[DF2$G == 2, ]
    DFyr$fill[DFyr$id == taxon] <- DF3$fill
  }

  DFmo$year = NA
  scaledprobs = dplyr::bind_rows(DFmo, DFyr)
  #  list(Pscaled_month = DFmo, Pscaled_year = DFyr)
  return(scaledprobs)
}

thresh_fig <- function(
    threshold,
    scaledprobs,
    species_name = NULL) {

  thresh_slc <- as.character(seq(50, 95, 5))
  threshold <- match.arg(arg = threshold, choices = thresh_slc)
  thresh <- data.frame(
    values = paste0("thresh", thresh_slc),
    labels = thresh_slc
  )

  thresh.value <- switch(threshold,
                         thresh$values[thresh$labels == threshold]
  )

  df <- scaledprobs %>%
    dplyr::filter(is.na(year)) %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      nd = sum(detect, na.rm = TRUE),
      n = sum(detect, nondetect, na.rm = TRUE),
      fill = mean(fill, na.rm = TRUE),
      scaleP = mean(scaleP, na.rm = TRUE)
    ) %>%
    dplyr::mutate(
      thresh95 = (fill >= 0.94999) * 1,
      thresh90 = (fill >= 0.89999) * 1,
      thresh85 = (fill >= 0.84999) * 1,
      thresh80 = (fill >= 0.79999) * 1,
      thresh75 = (fill >= 0.74999) * 1,
      thresh70 = (fill >= 0.69999) * 1,
      thresh65 = (fill >= 0.64999) * 1,
      thresh60 = (fill >= 0.59999) * 1,
      thresh55 = (fill >= 0.54999) * 1,
      thresh50 = (fill >= 0.49999) * 1
    )

  ggplot2::ggplot(df) +
    ggplot2::geom_hline(
      mapping = ggplot2::aes(yintercept = y),
      data.frame(y = c(0:4) / 4),
      color = "darkgrey"
    ) +
    ggplot2::geom_vline(
      mapping = ggplot2::aes(xintercept = x),
      data.frame(x = 1:12),
      color = "lightgrey"
    ) +
    ggplot2::geom_col(
      dplyr::filter(df, !!dplyr::ensym(thresh.value) %in% "1"),
      mapping = ggplot2::aes(x = month, y = fill), fill = viridis::viridis(1), position = "dodge2",
      width = 0.9, show.legend = FALSE, alpha = .9
    ) +
    ggplot2::geom_col(
      dplyr::filter(df, !!dplyr::ensym(thresh.value) %in% "0"),
      mapping = ggplot2::aes(x = month, y = fill), fill = "darkgrey", position = "dodge2",
      width = 0.9, show.legend = FALSE, alpha = .9
    ) +
    # To make the interpolated data stand out
    ggpattern::geom_col_pattern(
      dplyr::filter(df, is.nan(scaleP)),
      mapping = ggplot2::aes(x = month, y = fill),
      position = "dodge2",
      width = 0.9, show.legend = FALSE, fill = NA,
      pattern_color = "white", pattern_density = 0.05, pattern_spacing = 0.015, pattern_key_scale_factor = 0.6
    ) +
    ggplot2::coord_polar() +
    ggplot2::scale_x_continuous(
      limits = c(0.5, 12.5),
      breaks = 1:12,
      labels = month.abb
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      caption = species_name
    ) +
    ggplot2::labs(
      x = NULL, y = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.caption = element_text(hjust = 0.5, size = 32),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    ) +
    theme_circle
}


# custom theme for all figures
theme_circle <- ggplot2::theme(
  # plotting components

  ## drop minor gridlines
  panel.grid = ggplot2::element_blank(),
  # change grid lines to gray
  #  panel.grid.major =  element_line(color = "#d0d0d0"),
  # fill the plot and panel spaces with grey and remove border
  #  panel.background = element_blank(),
  # plot.background = element_blank(),
  panel.border = ggplot2::element_blank(),
  # adjust the margins of plots and remove axis ticks
  # plot.margin = ggplot2::margin(b = 25, l = 25,
  #                              unit = "pt"),
  axis.ticks = ggplot2::element_blank(),
  # change text family, size, and adjust position of titles
  text = ggplot2::element_text(
    family = "sans", size = 24),
  axis.text = ggplot2::element_text(
    colour = "#939598", size = 20),
  axis.text.y = ggplot2::element_blank(),
  axis.title = ggplot2::element_text(colour = "#5A5A5A",
                                     size = 24),
  axis.line = ggplot2::element_blank(),
  plot.title = ggplot2::element_text(
    face = "bold",
    size = 30,
    hjust = 0,
    colour = "#5A5A5A",
    margin = ggplot2::margin(b = 5, unit = "pt")),
  plot.title.position = "plot",
  plot.subtitle = ggplot2::element_text(size = 24,
                                        margin = ggplot2::margin(b = 25, unit = "pt"),
                                        colour = "#5A5A5A",
                                        hjust = 0),
  legend.title.align = 1,
  legend.text = ggplot2::element_text(size = 20,
                                      colour = "#939598"),
  legend.position = "right",
  legend.box.just = "right",
  legend.key.spacing.y = ggplot2::unit(20, "pt"),
  legend.spacing.y = ggplot2::unit(20, "pt"),
  legend.title = ggplot2::element_text(colour = "#5A5A5A",
                                       margin = ggplot2::margin(b = 20))

)


smooth_fig <- function(scaledprobs, species_name = NULL) {
  # Use the month-level scaled probability directly
  data <- scaledprobs %>%
    dplyr::filter(!is.na(scaleP), !is.na(year)) %>%
    dplyr::mutate(month = as.numeric(month),
                  year = as.factor(year))  # ensure year is factor for coloring

  Dsummary24 <- Dsummary12 <- data
  Dsummary12$month <- data$month + 12
  Dsummary24$month <- data$month + 24

  Dsummary_comb <- rbind(data, Dsummary12, Dsummary24)

  loessmod <- loess(scaleP ~ month, Dsummary_comb, span = 3 / 12)

  NEW <- data.frame(month = seq(1, 35, 0.1))
  NEW$PRED <- predict(loessmod, newdata = NEW$month)

  NEW2 <- NEW[NEW$month > 12 & NEW$month <= 24, ]
  NEW2$month <- NEW2$month - 12

  ggplot2::ggplot() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = y), data.frame(y = c(0:4)/4), color = "lightgrey") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = x), data.frame(x = 0:12), color = "lightgrey") +
    ggplot2::geom_path(data = NEW2, ggplot2::aes(x = month, y = PRED), colour = "blue") +
    ggplot2::geom_point(data = data, ggplot2::aes(x = month, y = scaleP, col = year), alpha = 0.9, size = 5) +
    ggplot2::coord_polar(clip = "off") +
    ggplot2::labs(col = "Year", x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_manual(values = palette("Alphabet")) +
    ggplot2::scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12, labels = month.abb) +
    ggplot2::guides(colour = ggplot2::guide_legend(label.position = "left", label.hjust = 1)) +
    ggplot2::scale_y_continuous(limits = c(-0.1, 1.01), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      caption = species_name
    ) +
    theme_circle
}



calc_window <- function(threshold, scaledprobs) {
  oop <- options("dplyr.summarise.inform")
  options(dplyr.summarise.inform = FALSE)
  # reset option on exit
  on.exit(options(dplyr.summarise.inform = oop))

  thresh_slc <- seq(50, 95, 5) %>% as.character()
  threshold <- match.arg(threshold, choices = thresh_slc)
  thresh <- data.frame(
    values = 0.49999 + seq(0, 0.45, 0.05),
    labels = thresh_slc
  )

  thresh.value <- thresh$values[thresh$labels == threshold]

  df <- scaledprobs %>%
    dplyr::filter(is.na(year)) %>%

    dplyr::group_by(month) %>%
    dplyr::summarise(
      nd = sum(det_int, na.rm = TRUE),
      n = sum(det_int, nd_int, na.rm = TRUE),
      fill = mean(fill, na.rm = TRUE),
      scaleP = mean(scaleP, na.rm = TRUE)
    )

  df_thresh <- df %>%
    dplyr::filter(
      fill >= thresh.value
    )

  # selects only a single month >= threshold
  multmonth <- df_thresh %>%
    dplyr::group_by(
      consec = cumsum(c(1, diff(month) != 1)))


  # selects species with consecutive months being December-January >= threshold
  if (length(unique(multmonth$consec)) == 2) {
    if (any(multmonth$month == 12) & any(multmonth$month == 1)) {
      consecmonth <- multmonth[order(multmonth$consec, decreasing = TRUE),]
      # all of the species with consecutive windows AND single month window
    } else {
      print("returning NULL at this place.")
      consecmonth <- NULL
    }

    consec.det <- consecmonth

  } else {

    if (length(unique(multmonth$consec)) == 1) {
      # selects species with consecutive months >= threshold
      consecmonth1 <- multmonth %>%
        dplyr::group_by(
          consec) %>%
        dplyr::ungroup()

      consec.det <- consecmonth1
    } else {
      print("there is not one window...returning NULL.")
      consec.det <- NULL
    }
  }

  # create df where species have no discernible window (for now, this means more than one non-consecutive period)
  if(!is.null(consec.det)){
    optwin <- consec.det %>%
      dplyr::mutate(window = "inwindow") %>%
      dplyr::ungroup()

    inwindow <- optwin %>%
      dplyr::group_by(
        window) %>%
      dplyr::summarise(
        detect = sum(nd, na.rm = TRUE),
        nondetect = sum(n, na.rm = TRUE)
      ) %>%
      dplyr::ungroup()

    # filter out all of the species that have a window from the dataset
    nowin <- dplyr::anti_join(df, optwin, by = "month") %>%
      dplyr::mutate(window = "outsidewindow")

    outsidewindow <- nowin %>%
      dplyr::group_by(
        window) %>%
      dplyr::summarise(
        detect = sum(nd, na.rm = TRUE),
        nondetect = sum(n, na.rm = TRUE)
      )

    window_sum <- dplyr::bind_rows(inwindow, outsidewindow)

    # Fisher's exact test to compare detection probability within/outside the window for each species and primer
    if (nrow(window_sum) == 2) {

      opt_sampling <- optwin %>%
        dplyr::summarise(
          len = length(month),
          thresh = unique(threshold),
          period = unique(dplyr::case_when(
            len == 1 ~ paste(month.abb[month]),
            len != 1 ~ paste0(month.abb[dplyr::first(month)], "-",
                              month.abb[dplyr::last(month)]))
          ))

      f <- fisher.test(matrix(
        c(
          window_sum$detect[window_sum$window == "inwindow"],
          window_sum$detect[window_sum$window == "outsidewindow"],
          window_sum$nondetect[window_sum$window == "inwindow"],
          window_sum$nondetect[window_sum$window == "outsidewindow"]
        ),
        nrow = 2
      ))

      fshTest <- data.frame(
        "odds ratio" = f$estimate,
        "p value" = scales::pvalue(f$p.value),
        "Lower CI" = f$conf.int[1],
        "Upper CI" = f$conf.int[2],
        check.names = FALSE
      ) %>%
        dplyr::mutate(
          confidence = factor(
            dplyr::case_when(
              `p value` == "<0.001" ~ "Very high",
              `p value` < "0.01" ~ "High",
              `p value` < "0.05" ~ "Medium",
              `p value` >= "0.05" ~ "Low",
              is.na(`p value`) ~ "No optimal period"),
            levels = c("Very high", "High", "Medium", "Low", "No optimal period"),
            ordered = TRUE))

      return(
        list(
          opt_sampling = opt_sampling,
          fshTest = fshTest))

    }


  } else {
    warning("No optimal detection window")
    return(NULL)

  }
}

calc_window_just_duals <- function(threshold, scaledprobs) {
  oop <- options("dplyr.summarise.inform")
  options(dplyr.summarise.inform = FALSE)
  on.exit(options(dplyr.summarise.inform = oop))

  # Threshold mapping
  thresh_slc <- seq(50, 95, 5) %>% as.character()
  threshold <- match.arg(threshold, choices = thresh_slc)
  thresh <- data.frame(
    values = 0.49999 + seq(0, 0.45, 0.05),
    labels = thresh_slc
  )
  thresh.value <- thresh$values[thresh$labels == threshold]

  # Aggregate per month
  df <- scaledprobs %>%
    dplyr::filter(is.na(year)) %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      nd = sum(det_int, na.rm = TRUE),
      n = sum(det_int, nd_int, na.rm = TRUE),
      fill = mean(fill, na.rm = TRUE),
      scaleP = mean(scaleP, na.rm = TRUE),
      .groups = "drop"
    )

  # Select months exceeding threshold
  df_thresh <- df %>% dplyr::filter(fill >= thresh.value)
  if (nrow(df_thresh) == 0) {
    warning("No months exceed threshold")
    return(NULL)
  }

  # Identify consecutive month blocks
  df_thresh <- df_thresh %>%
    dplyr::arrange(month) %>%
    dplyr::mutate(window_id = cumsum(c(1, diff(month) != 1)))

  # Only keep species with exactly 2 windows
  n_windows <- length(unique(df_thresh$window_id))
  if (n_windows != 2) {
    message("Skipping: not exactly two windows (found ", n_windows, ")")
    return(NULL)
  }

  # Summarise counts per window
  window_counts <- df_thresh %>%
    dplyr::group_by(window_id) %>%
    dplyr::summarise(
      detect = sum(nd, na.rm = TRUE),
      nondetect = sum(n, na.rm = TRUE),
      months = paste0(month.abb[month], collapse = "-"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(window_label = paste0("window_", row_number()))

  # Fisher's exact test between the two windows
  f <- fisher.test(matrix(
    c(
      window_counts$detect[1],
      window_counts$detect[2],
      window_counts$nondetect[1],
      window_counts$nondetect[2]
    ),
    nrow = 2
  ))

  fshTest <- data.frame(
    "odds ratio" = f$estimate,
    "p value" = scales::pvalue(f$p.value),
    "Lower CI" = f$conf.int[1],
    "Upper CI" = f$conf.int[2],
    check.names = FALSE
  )

  return(list(
    window_counts = window_counts,
    fshTest = fshTest,
    threshold = threshold
  ))
}



calc_window_dual <- function(threshold, scaledprobs) {
  oop <- options("dplyr.summarise.inform")
  options(dplyr.summarise.inform = FALSE)
  on.exit(options(dplyr.summarise.inform = oop))

  thresh_slc <- seq(50, 95, 5) %>% as.character()
  threshold <- match.arg(threshold, choices = thresh_slc)
  thresh <- data.frame(
    values = 0.49999 + seq(0, 0.45, 0.05),
    labels = thresh_slc
  )
  thresh.value <- thresh$values[thresh$labels == threshold]

  df <- scaledprobs %>%
    dplyr::filter(is.na(year)) %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      nd = sum(det_int, na.rm = TRUE),
      n = sum(det_int, nd_int, na.rm = TRUE),
      fill = mean(fill, na.rm = TRUE),
      scaleP = mean(scaleP, na.rm = TRUE),
      .groups = "drop"
    )

  df_thresh <- df %>% dplyr::filter(fill >= thresh.value)

  if (nrow(df_thresh) == 0) {
    warning("No months exceed the threshold")
    return(NULL)
  }

  # Identify consecutive month blocks
  df_thresh <- df_thresh %>%
    dplyr::arrange(month) %>%
    dplyr::mutate(consec = cumsum(c(1, diff(month) != 1)))

  # Combine all months in any window into a single "inwindow"
  inwindow_months <- df_thresh$month
  optwin <- df %>%
    dplyr::mutate(
      window = ifelse(month %in% inwindow_months, "inwindow", "outsidewindow")
    )

  # Summarise for Fisher test
  window_sum <- optwin %>%
    dplyr::group_by(window) %>%
    dplyr::summarise(
      detect = sum(nd, na.rm = TRUE),
      nondetect = sum(n, na.rm = TRUE),
      .groups = "drop"
    )

  # Only perform Fisher test if there are exactly two categories
  if (nrow(window_sum) == 2) {

    opt_sampling <- optwin %>%
      dplyr::filter(window == "inwindow") %>%
      dplyr::summarise(
        len = length(month),
        thresh = unique(threshold),
        period = if (len == 1) month.abb[month] else paste0(month.abb[first(month)], "-", month.abb[last(month)])
      )

    f <- fisher.test(matrix(
      c(
        window_sum$detect[window_sum$window == "inwindow"],
        window_sum$detect[window_sum$window == "outsidewindow"],
        window_sum$nondetect[window_sum$window == "inwindow"],
        window_sum$nondetect[window_sum$window == "outsidewindow"]
      ),
      nrow = 2
    ))

    fshTest <- data.frame(
      "odds ratio" = f$estimate,
      "p value" = scales::pvalue(f$p.value),
      "Lower CI" = f$conf.int[1],
      "Upper CI" = f$conf.int[2],
      check.names = FALSE
    ) %>%
      dplyr::mutate(
        confidence = factor(
          dplyr::case_when(
            `p value` == "<0.001" ~ "Very high",
            `p value` < "0.01" ~ "High",
            `p value` < "0.05" ~ "Medium",
            `p value` >= "0.05" ~ "Low",
            is.na(`p value`) ~ "No optimal period"
          ),
          levels = c("Very high", "High", "Medium", "Low", "No optimal period"),
          ordered = TRUE
        )
      )

    return(
      list(
        opt_sampling = opt_sampling,
        fshTest = fshTest
      )
    )
  } else {
    warning("No optimal detection window or more than one non-consecutive period")
    return(NULL)
  }
}



jaccard_test <- function(scaledprobs, threshold) {

  thresh_slc <- seq(50, 95, 5) %>% as.character()
  threshold <- match.arg(threshold, choices = thresh_slc)
  thresh <- data.frame(values = 0.49999 + seq(0, 0.45, 0.05),
                       labels = thresh_slc)

  thresh.value <- thresh$values[thresh$labels == threshold]

  data <- scaledprobs %>%
    dplyr::filter(!is.na(year)) %>%
    dplyr::group_by(year, month) %>%
    dplyr::summarise(
      detect = sum(detect, na.rm = TRUE),
      n = sum(detect, nondetect, na.rm = TRUE),
      fill = mean(fill, na.rm = TRUE)
    )

  df <- data %>%
    dplyr::mutate(
      fill_bin = case_when(
        fill < thresh.value ~ NA,
        fill >= thresh.value ~ 1
      )) %>%
    tidyr::pivot_wider(id_cols = year,
                       names_from = month,
                       values_from = fill_bin) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var="year") %>%
    as.matrix()

  n_year = data %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(n = sum(n, na.rm=TRUE))


  # df = df[purrr::map(df, nrow) > 1]

  n_mths <- df %>%
    as.data.frame() %>%
    #    mutate(year = rownames(.)) %>%
    mutate(
      n_mths = select(., 1:12) %>%
        rowSums(na.rm = TRUE)) %>%
    select(#year,
      n_mths)

  if (nrow(df) > 1) {

    compare = t(combn(nrow(df), 2, FUN=function(x) df[x[1],] == df[x[2],]))

  } else {
    compare = NA
  }

  if (all(is.na(compare)) == FALSE ){

    n_union = data.frame(
      n_union = combn(nrow(n_mths), 2, FUN=function(x)
        max(sum(df[x[1],], na.rm = TRUE),
            sum(df[x[2],], na.rm = TRUE))),
      n_wt = combn(nrow(n_year), 2, FUN=function(x)
        sum(n_year[x[1],"n"], n_year[x[2],"n"]))
    )

    rownames(compare) = combn(nrow(df), 2, FUN=function(x) paste0("year",x[1],"_year",x[2]))

    jacc_df <- compare %>%
      as.data.frame() %>%
      mutate(group = rownames(.)) %>%
      rowwise() %>%
      mutate(
        n_int = #select(., V1:V12) %>%
          rowSums(pick(V1:V12), na.rm = TRUE)) %>%
      select(group, n_int) %>%
      cbind(n_union) %>%
      mutate(
        J_index = n_int/n_union,
        J_wt = J_index*n_wt
      )


    jacc_wt_mean <- jacc_df %>%
      dplyr::summarise(
        wt_mean = sum(J_wt)/sum(n_wt) *100
      ) %>%
      dplyr::mutate(
        wt_text = dplyr::case_when(
          between(wt_mean, 0, 9.99) == TRUE ~ "Very low",
          between(wt_mean, 10, 29.99) == TRUE ~ "Low",
          between(wt_mean, 30, 69.99) == TRUE ~ "Medium",
          between(wt_mean, 70, 89.99) == TRUE ~ "High",
          between(wt_mean, 90, 100) == TRUE ~ "Very high"
        )
      )

    return(c(jacc_wt_mean$wt_mean, jacc_wt_mean$wt_text))
  } else {
    return(paste("Multi-year data not available"))
  }
}



#############################################
source('./AIS_eDNA_data_prep.R')
df <- dfRawClean %>% mutate(protocol_ID = 1, samp_name = materialSampleID)
df$All <- "ALL AIS SPECIES"


all_species <- unique(df$scientificName)
selected_species <- all_species[1]
all_regions <- unique(df$region)
selected_region <- all_regions[5]

df_filtered <- df %>% filter(scientificName == selected_species, region == selected_region)




new_probs <- calc_det_prob(df_filtered, selected_taxon_level = "scientificName", selected_taxon_id = selected_species, pool_primers = TRUE)
scaledprobs <- scale_newprob(df_filtered, new_probs, selected_taxon_level = "scientificName")

thresh_fig("75", scaledprobs, selected_species)
smooth_fig(scaledprobs, selected_species)


win <- calc_window("75", scaledprobs)
jaccard_results <- jaccard_test(scaledprobs, "75")


##############################
#ALL WHEELS IN ONE FIGURE
##############################

for (sp in all_species) {
  for (reg in all_regions) {

    df_sub <- df %>%
      filter(scientificName == sp, region == reg)

    # Skip completely empty datasets (optional)
    if (nrow(df_sub) == 0) next

    new_probs <- calc_det_prob(
      df_sub,
      selected_taxon_level = "scientificName",
      selected_taxon_id = sp,
      pool_primers = TRUE
    )

    # Skip if no usable output
    if (length(new_probs$newP_agg) == 0) next

    scaledprobs <- scale_newprob(
      df_sub,
      new_probs,
      selected_taxon_level = "scientificName"
    )

    p <- thresh_fig("75", scaledprobs, paste(sp, reg))

    # Clean filename (important!)
    file_name <- paste0(
      "thresh_plots/",
      gsub(" ", "_", sp), "_",
      reg,
      ".png"
    )

    ggsave(
      filename = file_name,
      plot = p,
      width = 6,
      height = 6,
      units = "in",
      dpi = 300
    )
  }
}
