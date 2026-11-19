# Map vector layers: flow as particle-animation data and as arrow spokes.
# Both take a flow's east (u) and north (v) components as rasters.

format_velocity_data <- function(u, v) {

      ext <- terra::ext(u)
      dims <- dim(u)

      # Extract values row-major, top to bottom
      u_vals <- as.vector(t(terra::as.matrix(u, wide = TRUE)))
      v_vals <- as.vector(t(terra::as.matrix(v, wide = TRUE)))

      # Replace NA with 0
      u_vals[is.na(u_vals)] <- 0
      v_vals[is.na(v_vals)] <- 0

      list(
            list(
                  header = list(
                        parameterCategory = 2,
                        parameterNumber = 2,
                        lo1 = ext$xmin,
                        la1 = ext$ymax,
                        lo2 = ext$xmax,
                        la2 = ext$ymin,
                        dx = terra::res(u)[1],
                        dy = terra::res(u)[2],
                        nx = dims[2],
                        ny = dims[1]
                  ),
                  data = u_vals
            ),
            list(
                  header = list(
                        parameterCategory = 2,
                        parameterNumber = 3,
                        lo1 = ext$xmin,
                        la1 = ext$ymax,
                        lo2 = ext$xmax,
                        la2 = ext$ymin,
                        dx = terra::res(u)[1],
                        dy = terra::res(u)[2],
                        nx = dims[2],
                        ny = dims[1]
                  ),
                  data = v_vals
            )
      ) |>
            jsonlite::toJSON(auto_unbox = TRUE, digits = 6)
}

# Components are averaged over blocks of `subsample` cells, so opposing flows
# within a block cancel, as they do in the flow itself
make_spokes <- function(u, v, subsample = 10, arrow_scale = 0.5) {
      u_sub <- terra::aggregate(u, fact = subsample, fun = "mean", na.rm = TRUE)
      v_sub <- terra::aggregate(v, fact = subsample, fun = "mean", na.rm = TRUE)

      coords <- terra::xyFromCell(u_sub, 1:ncell(u_sub))
      du <- values(u_sub, mat = FALSE)
      dv <- values(v_sub, mat = FALSE)

      valid <- !is.na(du) & !is.na(dv)

      spokes <- data.frame(
            x0 = coords[valid, 1],
            y0 = coords[valid, 2],
            x1 = coords[valid, 1] + du[valid] * arrow_scale,
            y1 = coords[valid, 2] + dv[valid] * arrow_scale,
            mag = sqrt(du[valid]^2 + dv[valid]^2)
      )

      lines <- lapply(seq_len(nrow(spokes)), function(i) {
            st_linestring(matrix(c(spokes$x0[i], spokes$y0[i],
                                   spokes$x1[i], spokes$y1[i]),
                                 ncol = 2, byrow = TRUE))
      })

      points <- lapply(seq_len(nrow(spokes)), function(i) {
            st_point(c(spokes$x0[i], spokes$y0[i]))
      })

      list(
            lines = st_sf(mag = spokes$mag, geometry = st_sfc(lines, crs = 4326)),
            points = st_sf(mag = spokes$mag, geometry = st_sfc(points, crs = 4326))
      )
}
