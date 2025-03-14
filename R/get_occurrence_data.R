#' Get Occurrence Data
#'
#' This function imports and processes biodiversity occurrence data from various sources.
#' It supports local CSV files, data frames already in R, or GBIF downloads (in ZIP format).
#'
#' @param data Input data. Can be:
#'   - File path to a local CSV file (when source_type = "local_csv").
#'   - A data frame already loaded in R (when source_type = "data_frame").
#'   - NULL for GBIF download (when source_type = "gbif").
#' @param source_type Source type of the data. One of:
#'   - "local_csv": Load data from a CSV file.
#'   - "data_frame": Use an existing data frame in R.
#'   - "gbif": Download and process data from a GBIF ZIP URL.
#' @param gbif_zip_url URL of the GBIF ZIP file (required when source_type = "gbif").
#' @param download_dir Directory to save the downloaded GBIF data (default: temporary directory).
#' @param sep Separator for reading CSV files (default: ",").
#'
#' @return A processed data frame in long or wide format. Returns:
#'   - "long format" if sp_name and either pa or abundance columns are detected.
#'   - "wide format" if species columns (sp_*) are present.
#'
#' @export
#'
#' @examples
#' # Load data from a local CSV file:
#' # local_data <- get_occurrence_data(data = "path/to/local.csv", source_type = "local_csv")
#'
#' # Use an existing data frame:
#' # df_data <- get_occurrence_data(data = local_data, source_type = "data_frame")
#'
#' # Download GBIF data:
#' # gbif_data <- get_occurrence_data(source_type = "gbif",
#' #                                  gbif_zip_url = "https://api.gbif.org/v1/occurrence/download/request/0038969-240906103802322.zip",
#' #                                  download_dir = "path/to/download/dir")
get_occurrence_data <- function(data = NULL,
                                source_type = c("local_csv", "data_frame", "gbif"),
                                gbif_zip_url = NULL,
                                download_dir = tempdir(),
                                sep = ",") {
  # Ensure required packages are loaded
  required_packages <- c("dplyr", "tidyr", "httr", "data.table")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required but not installed.")
  }

  source_type <- match.arg(source_type)

  # Helper function to detect columns
  detect_columns <- function(data, possible_names) {
    matched <- intersect(tolower(names(data)), tolower(possible_names))
    if (length(matched) > 1) {
      preferred_order <- match(tolower(possible_names), tolower(matched))
      matched <- matched[order(preferred_order)]
    }
    if (length(matched) > 1) matched <- matched[1]
    if (length(matched) == 1)
      return(names(data)[tolower(names(data)) %in% matched])
    NULL
  }

  # Read data based on source_type
  if (source_type == "local_csv") {
    if (!is.character(data) || !file.exists(data))
      stop("Invalid file path provided for 'local_csv'.")
    data <- utils::read.csv(data, sep = sep, stringsAsFactors = FALSE,
                            row.names = NULL, check.names = FALSE)
  } else if (source_type == "data_frame") {
    if (!is.data.frame(data))
      stop("For 'data_frame', 'data' must be a valid data frame.")
  } else if (source_type == "gbif") {
    if (is.null(gbif_zip_url))
      stop("For 'gbif', 'gbif_zip_url' is required.")
    if (!dir.exists(download_dir))
      dir.create(download_dir, recursive = TRUE)

    zip_path <- file.path(download_dir, basename(gbif_zip_url))
    httr::GET(gbif_zip_url, httr::write_disk(zip_path, overwrite = TRUE))

    csv_files <- utils::unzip(zip_path, list = TRUE)
    occurrence_file <- csv_files$Name[grepl("(occurrence\\.txt|\\.csv)$",
                                            csv_files$Name, ignore.case = TRUE)]
    if (length(occurrence_file) == 0)
      stop("No valid occurrence file found in the ZIP.")

    utils::unzip(zip_path, files = occurrence_file, exdir = download_dir)
    data <- data.table::fread(file.path(download_dir, occurrence_file),
                              sep = "\t", data.table = FALSE)
  }

  # Detect key columns
  column_sets <- list(
    site_id = c("site_id", "site", "sample", "id", "plot"),
    x       = c("x", "lon", "long", "longitude", "x_coord", "decimalLongitude"),
    y       = c("y", "lat", "latitude", "y_coord", "decimalLatitude"),
    sp_name = c("sp_name", "name", "species", "scientific", "spp", "verbatimScientificName"),
    pa      = c("pa", "presence", "obs"),
    abund   = c("abund", "abundance", "count", "total")
  )
  detected_columns <- lapply(column_sets, detect_columns, data = data)

  # Validate coordinates
  if (is.null(detected_columns$x) || is.null(detected_columns$y))
    stop("Latitude/Longitude columns not detected.")
  if (any(data[[detected_columns$x]] < -180 | data[[detected_columns$x]] > 180,
          na.rm = TRUE)) {
    stop("Longitude values should be between -180 and 180.")
  }
  if (any(data[[detected_columns$y]] < -90 | data[[detected_columns$y]] > 90,
          na.rm = TRUE)) {
    stop("Latitude values should be between -90 and 90.")
  }

  # Create site_id if missing
  if (is.null(detected_columns$site_id)) {
    unique_coords <- unique(data[, c(detected_columns$x, detected_columns$y), drop = FALSE])
    unique_coords$site_id <- seq_len(nrow(unique_coords))
    data <- dplyr::left_join(data, unique_coords, by = c(detected_columns$x, detected_columns$y))
    detected_columns$site_id <- "site_id"
  }

  # If no value column specified, create a presence-absence column
  if (is.null(detected_columns$pa) && is.null(detected_columns$abund)) {
    data <- data %>% dplyr::mutate(pa = 1)
    message("No value column specified. A presence-absence column ('pa') was created and filled with 1.")
    detected_columns$pa <- "pa"
  }

  # Determine format and process
  sp_columns <- grep("^sp_", names(data), value = TRUE)
  is_long <- !is.null(detected_columns$sp_name) &&
    (!is.null(detected_columns$pa) || !is.null(detected_columns$abund))
  is_wide <- length(sp_columns) > 0

  if (is_long) {
    data <- data %>%
      dplyr::rename(site_id = detected_columns$site_id,
                    x = detected_columns$x,
                    y = detected_columns$y,
                    sp_name = detected_columns$sp_name)
    if (!is.null(detected_columns$pa)) {
      data <- data %>%
        dplyr::mutate(pa = ifelse(is.na(data[[detected_columns$pa]]) |
                                    data[[detected_columns$pa]] == "",
                                  1, suppressWarnings(as.numeric(data[[detected_columns$pa]]))))
    } else {
      data <- data %>%
        dplyr::mutate(abund = suppressWarnings(as.numeric(data[[detected_columns$abund]])))
    }
    return(data)
  } else if (is_wide) {
    data <- data %>%
      tidyr::pivot_longer(cols = sp_columns, names_to = "sp_name", values_to = "value")
    return(data)
  } else {
    stop("Unable to determine data format.")
  }
}
