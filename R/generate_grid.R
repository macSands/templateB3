#' Generate Spatial Grid
#'
#' This function generates a spatial grid over a geographic extent and assigns grid IDs to points.
#' It also computes cell centroids and, for geographic grids (in degrees or minutes), a mapsheet code.
#' Finally, it summarizes specified columns within each grid cell while preserving additional metadata
#' provided via `extra_cols`.
#'
#' @param data A data frame containing point data with x and y coordinates.
#' @param x_col Character. Column name for x-coordinates (default: "x").
#' @param y_col Character. Column name for y-coordinates (default: "y").
#' @param grid_size Numeric. Size of the grid cells. For geographic data (EPSG:4326) this is in degrees.
#' @param sum_col_range Numeric vector. Range of columns to summarize within each grid cell.
#' @param extra_cols Character vector of additional columns to retain in the output (optional).
#' @param crs_epsg Numeric. EPSG code for the coordinate reference system (default: 4326).
#' @param unit Character. One of "deg", "min", "sec", or "m". For geographic data use "deg" (default).
#'
#' @return A list containing:
#'   - `grid`: terra raster object of the generated grid with unique grid IDs.
#'   - `grid_sf`: sf object of the grid as polygons, with added centroid coordinates and (if applicable) mapsheet codes.
#'   - `block_sp`: Data frame summarizing grid cell contents (including extra_cols if provided).
#'
#' @import sf
#' @import terra
#' @import dplyr
#' @export
#'
#' @examples
#' # Simulated Example for generate_grid
#' set.seed(123)
#' data <- data.frame(
#'   x = runif(100, -10, 10),
#'   y = runif(100, -10, 10),
#'   species1 = rpois(100, 5),
#'   species2 = rpois(100, 3),
#'   recordedBy = sample(LETTERS, 100, replace = TRUE)
#' )
#'
#' # Generate grid with grid cells of size 1 degree,
#' # summarizing species counts (columns 3 and 4)
#' # and retaining the 'recordedBy' metadata.
#' grid_result <- generate_grid(
#'   data,
#'   x_col = "x",
#'   y_col = "y",
#'   grid_size = 1,
#'   sum_col_range = 3:4,
#'   extra_cols = "recordedBy",
#'   unit = "deg"
#' )
#'
#' # Inspect the summarized block spatial data
#' print(grid_result$block_sp)
#'
#' # Plot the grid polygons with grid IDs
#' plot(grid_result$grid_sf["grid_id"], main = "Grid Polygons with IDs")
#'
#' # Optionally, plot the terra raster grid
#' terra::plot(grid_result$grid, main = "Terra Raster Grid")
generate_grid <- function(data,
                          x_col = "x",       # Column name for x-coordinates
                          y_col = "y",       # Column name for y-coordinates
                          grid_size = 0.5,   # Grid cell size (in degrees for EPSG:4326)
                          sum_col_range = NULL,  # Columns (by index) to summarize
                          extra_cols = NULL, # Additional columns to retain (if available)
                          crs_epsg = 4326,   # Coordinate reference system
                          unit = c("deg", "min", "sec", "m")
) {
  # Ensure required packages are available
  required_packages <- c("sf", "terra", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed.")
    }
  }

  unit <- match.arg(unit)

  # Validate inputs
  if (!x_col %in% names(data) || !y_col %in% names(data)) {
    stop("Specified x or y columns do not exist in the input data frame.")
  }
  if (is.null(sum_col_range)) {
    stop("Please specify `sum_col_range` for columns to sum.")
  }

  # Convert data to an sf object using the specified CRS
  points_sf <- sf::st_as_sf(data, coords = c(x_col, y_col), crs = crs_epsg)

  # Determine the bounding box of the points and extend it by 2 grid cells in each direction.
  bounds <- sf::st_bbox(points_sf)
  bounds["xmin"] <- floor(bounds["xmin"] / grid_size) * grid_size - 2 * grid_size
  bounds["ymin"] <- floor(bounds["ymin"] / grid_size) * grid_size - 2 * grid_size
  bounds["xmax"] <- ceiling(bounds["xmax"] / grid_size) * grid_size + 2 * grid_size
  bounds["ymax"] <- ceiling(bounds["ymax"] / grid_size) * grid_size + 2 * grid_size
  bb <- bounds  # bb now is a proper bbox object

  # --- Create grid as an sf object (integrating functionality from create_grid) ---
  grid_polygons <- sf::st_make_grid(sf::st_as_sfc(bb), cellsize = grid_size, what = "polygons")
  grid_sf <- sf::st_sf(geometry = grid_polygons, crs = crs_epsg)

  # Compute cell centroids and attach centroid coordinates
  centroids <- sf::st_centroid(grid_sf)
  coords <- sf::st_coordinates(centroids)
  grid_sf$centroid_lon <- coords[, 1]
  grid_sf$centroid_lat <- coords[, 2]

  # If using a graticule (degrees or minutes), add a mapsheet code.
  if (unit %in% c("deg", "min")) {
    lon_int <- floor(grid_sf$centroid_lon)
    lat_int <- floor(grid_sf$centroid_lat)
    lon_dir <- ifelse(grid_sf$centroid_lon >= 0, "E", "W")
    lat_dir <- ifelse(grid_sf$centroid_lat >= 0, "N", "S")
    sub_lon <- "B"  # Placeholder sub-cell indicator
    sub_lat <- "B"
    grid_sf$mapsheet <- sprintf("%s%03d%s%02d%s%s",
                                lon_dir, abs(lon_int),
                                lat_dir, abs(lat_int),
                                sub_lon, sub_lat)
  }

  # Assign unique grid IDs to the sf grid
  grid_sf$grid_id <- as.character(seq_len(nrow(grid_sf)))

  # --- Create a terra raster grid for additional operations (e.g., visualization) ---
  grid_rast <- terra::rast(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"],
    resolution = grid_size, crs = paste0("EPSG:", crs_epsg)
  )
  grid_rast[] <- seq_len(terra::ncell(grid_rast))
  names(grid_rast) <- "grid_id"

  # --- Assign grid IDs to points via spatial join ---
  # Points falling inside a grid cell inherit that cell's grid_id.
  points_join <- sf::st_join(points_sf, grid_sf[, c("grid_id")], join = sf::st_within)
  data$grid_id <- as.character(points_join$grid_id)

  # --- Summarize data within each grid cell ---
  sum_cols <- names(data)[sum_col_range]
  # If extra_cols are provided, include them as grouping variables
  if (!is.null(extra_cols)) {
    block_sp <- data %>%
      dplyr::group_by(grid_id, dplyr::across(dplyr::any_of(extra_cols))) %>%
      dplyr::summarize(across(all_of(sum_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  } else {
    block_sp <- data %>%
      dplyr::group_by(grid_id) %>%
      dplyr::summarize(across(all_of(sum_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  }

  # Merge summarized data with grid cell attributes (centroids and mapsheet code)
  grid_info <- grid_sf %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::select(grid_id, centroid_lon, centroid_lat, dplyr::any_of("mapsheet"))
  block_sp <- dplyr::left_join(block_sp, grid_info, by = "grid_id") %>%
    dplyr::select(grid_id, centroid_lon, centroid_lat, dplyr::any_of("mapsheet"), dplyr::everything())

  # row.names(block_sp) <- block_sp$grid_id

  return(list(grid = grid_rast,
              grid_sf = grid_sf,
              block_sp = as.data.frame(block_sp)))
}

