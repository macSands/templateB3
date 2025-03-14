#' Get Environmental Data
#'
#' This function retrieves environmental data for spatial points, either by downloading data from
#' specified sources (e.g., WorldClim, SoilGrids) or using local raster files. Data is extracted for
#' a buffer region around the input points, and missing values are interpolated using surrounding values.
#'
#' @param data A data frame containing spatial points with x and y coordinates.
#' @param buffer_km Numeric. Buffer distance around the input points in kilometers (default: 10).
#' @param var Character. Variable to retrieve (e.g., "bio" for bioclimatic variables, "elev" for elevation).
#' @param res Numeric. Resolution of the environmental data in arc-minutes (default: 2.5).
#' @param path Character. Directory path to store downloaded data or to load local rasters (default: "data/").
#' @param source Character. Data source: "geodata" for downloading or "local" for local rasters (default: "geodata").
#' @param year Numeric. Year for temporal data (e.g., for population or footprint datasets; required for specific datasets).
#' @param model Character. Climate model for CMIP6 projections (e.g., "ACCESS-CM2").
#' @param ssp Numeric. Shared Socioeconomic Pathway for CMIP6 (e.g., 245 for SSP2-4.5).
#' @param time Character. Time period for CMIP6 projections (e.g., "2041-2060").
#' @param depth Character. Depth range for soil data (e.g., "0-5cm").
#' @param stat Character. Statistic to summarize soil data (e.g., "mean").
#'
#' @return A list containing:
#'   - `env_rast`: Raster of the environmental data for the area of interest.
#'   - `sites_sf`: Spatial points as an sf object.
#'   - `env_df`: Data frame combining input data and extracted environmental variables.
#'
#' @import sf
#' @import terra
#' @import dplyr
#' @import geodata
#' @import zoo
#' @import tidyselect
#' @export
#'
#' @examples
#' # Example usage:
#' data = data.frame(site_id = 1:5, x = c(10, 12, 14, 16, 18), y = c(20, 22, 24, 26, 28))
#' env_data = get_enviro_data(data, buffer_km = 5, var = "bio", res = 2.5, path = "data/")
#' terra::plot(env_data$env_rast[[1]])
#' points(env_data$sites_sf)

get_enviro_data = function(data,
                           buffer_km = 10,
                           source = "geodata", # Options: 'local'
                           var = "bio", # Options: c('bio','elev','footprint','population','soil_world')
                           res = 2.5, # Resolution needed for 'worldclim_global','elev', 'population'
                           path = "data/",
                           year = NULL, # Needed for 'footprint' and 'population' data
                           depth = NULL, # Needed for 'soil_world'
                           stat = "mean", # Needed for 'soil_world'
                           model = NULL, # Needed for climate CMIP6 projections
                           ssp = NULL, # Needed for climate CMIP6 projections
                           time = NULL # Needed for climate CMIP6 projections
) {
  # Load required packages
  required_packages = c("terra", "sf", "dplyr", "geodata", "zoo", "tidyselect")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is required but not installed.")
  })

  # Ensure the output directory exists
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created directory: ", path)
  }

  # Validate input data
  if (!is.data.frame(data)) stop("Input must be a data frame.")

  # Identify coordinate columns
  x_cols = c("x", "lon", "longitude", "x_coord", "decimalLongitude")
  y_cols = c("y", "lat", "latitude", "y_coord", "decimalLatitude")
  site_cols = c("site_id", "grid_id")

  find_col = function(cols) {
    match = intersect(tolower(names(data)), tolower(cols))
    if (length(match) == 0) return(NULL)
    names(data)[tolower(names(data)) %in% match]
  }

  site_col = find_col(site_cols)
  x = find_col(x_cols)
  y = find_col(y_cols)
  if (is.null(x) || is.null(y)) stop("Coordinate columns (x, y) not found in the input data.")

  # Generate AOI
  message("Generating Area of Interest (AOI)...")
  data_xy = data %>% dplyr::select(all_of(c(site_col, x, y))) %>% dplyr::distinct()
  sites_sf = sf::st_as_sf(data_xy, coords = c(x, y), crs = 4326)
  aoi_extent = sf::st_buffer(sf::st_convex_hull(sf::st_union(sites_sf)), dist = buffer_km * 1000)

  # Download or load environmental data
  env_rast = switch(
    source,
    geodata = {
      message("Downloading data using geodata...")
      switch(
        var,
        bio = geodata::worldclim_global(var, res, path),
        elev = geodata::worldclim_global(var, res, path),
        footprint = geodata::footprint(year, path),
        population = geodata::population(year, res, path),
        soil_world = geodata::soil_world(var, depth, stat, path),
        stop("Unsupported variable for geodata source.")
      )
    },
    local = {
      message("Using local raster data...")
      if (inherits(var, "SpatRaster")) {
        var
      } else if (is.character(var) && file.exists(var)) {
        terra::rast(var)
      } else {
        stop("Invalid input for `var`. Must be a SpatRaster object or valid file path.")
      }
    },
    stop("Unsupported data source. Use 'geodata' or 'local'.")
  )

  if (!is.null(env_rast)) {
    env_rast = terra::crop(env_rast, terra::vect(aoi_extent))
    names(env_rast) = sub(".*_", "", names(env_rast))  # Simplify names
  }

  if (is.null(env_rast)) stop("Failed to load environmental data.")

  # Extract data for input points
  message("Extracting environmental data for input points...")
  env_data = terra::extract(env_rast, terra::vect(sites_sf))
  env_df = dplyr::bind_cols(data_xy, env_data)
  # Interpolate missing values
  env_df = env_df %>% dplyr::mutate(across(where(is.numeric), ~ zoo::na.approx(.x, na.rm = FALSE, rule = 2)))

  # # Extract raster values at site locations
  # env_data = terra::extract(env_rast, terra::vect(sites_sf))
  #
  # # Combine the site coordinates with the extracted data
  # env_df = dplyr::bind_cols(data_xy, env_data)
  # # print(str(env_df))
  #
  # # Interpolate missing values using nearby points
  # interpolate_missing_values = function(df) {
  #   # Select coordinates and numeric columns
  #   coords = df %>% dplyr::select(x, y)
  #   # print(head(coords))
  #   # numeric_cols = df %>% dplyr::select(where(is.numeric))
  #   numeric_cols <- df %>% dplyr::select(where(is.numeric), -ID) %>% names()
  #   # print(numeric_cols)
  #
  #   # Use spatial interpolation for each numeric column
  #   for (col in colnames(numeric_cols)) {
  #     missing_idx = which(is.na(df[[col]]))
  #     if (length(missing_idx) > 0) {
  #       # Create a spatial object for interpolation
  #       known_points = terra::vect(coords[!is.na(df[[col]]), ], crs = crs(env_rast))
  #       known_values = numeric_cols[[col]][!is.na(df[[col]])]
  #
  #       # Convert missing points to spatial object
  #       missing_points = terra::vect(coords[missing_idx, ], crs = crs(env_rast))
  #
  #       # Perform interpolation (bilinear or nearest neighbor)
  #       interpolated_values = terra::interpolate(known_points, known_values, missing_points, method = "bilinear")
  #
  #       # Assign interpolated values back to missing indices
  #       df[[col]][missing_idx] = interpolated_values
  #     }
  #   }
  #   return(df)
  # }
  #
  # # Apply interpolation function to fill missing values
  # env_df = interpolate_missing_values(env_df)

  # Save results in the global environment
  assign("env_rast", env_rast, envir = .GlobalEnv)
  assign("sites_sf", sites_sf, envir = .GlobalEnv)
  assign("env_df", env_df, envir = .GlobalEnv)

  return(list(env_rast = env_rast, sites_sf = sites_sf, env_df = env_df))
}
