# DOWNLOADS --------------------------------------
output$download_raster <- downloadHandler(
      filename = function() {
            paste(input$download_filename, ".tif", sep="")
      },

      content = function(file) {
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

            write_csv(pd, file)
      }
)
output$download_metadata <- downloadHandler(
      filename = function() {
            paste(input$download_filename, ".txt", sep="")
      },

      content = function(file) {
            params <- names(input)
            out <- c()
            for(p in params){
                  out <- c(out, paste(p, input[[p]], sep = ": "))
            }
            if(input$stat != "connectivity"){
                  out <- out[!grepl("cs_", params)]
            }
            writeLines(out, file)
      }
)
