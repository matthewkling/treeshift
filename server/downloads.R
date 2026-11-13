# DOWNLOADS --------------------------------------
output$download_raster <- downloadHandler(
      filename = function() {
            paste(input$download_filename, ".tif", sep="")
      },

      content = function(file) {
            log_event("download", what = "raster", query = log_query())
            terra::writeRaster(display_layer(), file)
      }
)
output$download_points <- downloadHandler(
      filename = function() {
            paste(input$download_filename, "_points.csv", sep="")
      },

      content = function(file) {
            # isolate() because this runs outside a render context, where a
            # req() reaching up through the reactive chain surfaces as a 500
            # (which the browser saves as download_points.html) rather than
            # as a quiet no-op. The realistic path here is "focal site
            # analogs" selected with no site clicked yet; a header-only file
            # is the truthful answer to that, and always a valid CSV.
            pd <- tryCatch(isolate(export_points()),
                           error = function(e) export_points_empty())
            log_event("download", what = "points", rows = nrow(pd), query = log_query())

            write_csv(pd, file)
      }
)
# The settings record (helpers/provenance.R)
output$download_record <- downloadHandler(
      filename = function() {
            paste0(input$download_filename, "_settings.txt")
      },

      content = function(file) {
            log_event("download", what = "settings record", query = log_query())
            writeLines(provenance_text(current_state(), url_ctx, link_url(), app_info), file)
      }
)
